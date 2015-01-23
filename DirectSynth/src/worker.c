/*
 * worker.c
 *
 *  Created on: Dec 11, 2014
 *      Author: scott
 */

#include "include.h"
#include "structure.h"
#include "functions.h"
#include "defs.h"


/*if my_id>N:
	Obtain process/point mapping
	while true:
		Request next available task
		if no more tasks:
			break
		For each rupture variation in task
			generate genslip
		Get list of SGTs needed
		Construct requests of (process, SGTs) (if request is too big, split it into pieces)
		for each request:
			Request data from process (<N)
			Perform convolution on data and SRF
		Calculate PSA, RotD
		Send output data to process N to aggregate and write
	exit
*/
void do_work(struct sgtmaster* sgtmast, struct sgtindex* sgtindx, struct sgtfileparams* sgtfilepar, int* proc_points, int num_sgt_handlers, char stat[64], float slat, float slon, int run_id, float det_max_freq, float stoch_max_freq, int run_PSA, int run_rotd, int my_id);
void get_task(worker_msg* w_msg, MPI_Datatype* worker_msg_type, int task_manager_id, int my_id);


int worker(int num_sgt_handlers, struct sgtfileparams* sgtfilepar, char stat[64], float slat, float slon, int run_id, float det_max_freq, float stoch_max_freq, int run_PSA, int run_rotd, int my_id) {
	//Construct MPI datatype
	MPI_Datatype sgtmast_type, sgtindx_type;
	construct_sgtmast_datatype(&sgtmast_type);
	construct_sgtindx_datatype(&sgtindx_type);

	//Get sgtmast, sgtindx info
	struct sgtmaster sgtmast;
	struct sgtindex sgtindx;
	if (debug) write_log("Receiving sgtmast, sgtindx from master.");
	check_bcast(&sgtmast, 1, sgtmast_type, 0, MPI_COMM_WORLD, "Error receiving sgtmast, aborting.", my_id);
	check_bcast(&sgtindx, sgtmast.globnp, sgtindx_type, 0, MPI_COMM_WORLD, "Error receiving sgtindx, aborting.", my_id);

	if (debug) close_log();
        MPI_Finalize();
        exit(0);

	//Get point-to-process mapping
	int* proc_points = check_malloc(sizeof(int)*(num_sgt_handlers+1));
	if (debug) write_log("Receiving SGT point-to-process mapping from master.");
	check_bcast(proc_points, num_sgt_handlers+1, MPI_INT, 0, MPI_COMM_WORLD, "Error receiving SGT point-to-process mapping, aborting.", my_id);

	//Now start doing tasks
	do_work(&sgtmast, &sgtindx, sgtfilepar, proc_points, num_sgt_handlers, stat, slat, slon, run_id, det_max_freq, stoch_max_freq, run_PSA, run_rotd, my_id);

	free(proc_points);
	return 0;
}


void do_work(struct sgtmaster* sgtmast, struct sgtindex* sgtindx, struct sgtfileparams* sgtfilepar, int* proc_points, int num_sgt_handlers, char stat[64],
		float slat, float slon, int run_id, float det_max_freq, float stoch_max_freq, int run_PSA, int run_rotd, int my_id) {
	int task_manager_id = num_sgt_handlers;

	worker_msg msg;
	MPI_Datatype worker_msg_type;
	task_info t_info;
	t_info.sgtmast = sgtmast;
	t_info.sgtindx = sgtindx;
	t_info.sgtfilepar = sgtfilepar;
	construct_worker_message_datatype(&msg);
	while (1) {
		if (debug) write_log("Requesting work from task manager.");
		get_task(&msg, &worker_msg_type, task_manager_id, my_id);
		if (msg.msg_type==WORK_RESPONSE) {
			t_info.task = &(msg.task);
			if (debug) {
				char buf[256];
				sprintf(buf, "Received task from task manager: source %d, rupture %d, starting rv %d, ending rv %d.", t_info.task->source_id, t_info.task->rupture_id, t_info.task->starting_rv_id, t_info.task->ending_rv_id);
				write_log(buf);
			}
			run_synth(&t_info, proc_points, num_sgt_handlers, stat, slat, slon, det_max_freq, stoch_max_freq, run_PSA, run_rotd, my_id);
		} else if (msg.msg_type==WORK_COMPLETE) {
			if (debug) write_log("No more work to do, shutting down.");
			break;
		} else {
			fprintf(stderr, "%d) SGT handler received message of unknown type %d from task manager, aborting.", my_id, msg.msg_type);
			if (debug) close_log();
			MPI_Finalize();
			exit(4);
		}
	}
}


void get_task(worker_msg* w_msg, MPI_Datatype* worker_msg_type, int task_manager_id, int my_id) {
	manager_msg msg;
	msg.msg_src = my_id;
	msg.msg_type = WORK_REQUEST;
	check_send(&msg, 2, MPI_INT, task_manager_id, WORK_REQUEST_TAG, MPI_COMM_WORLD, "Error requesting work from task manager, aborting.", my_id);
	check_recv(w_msg, 1, *worker_msg_type, task_manager_id, MPI_ANY_TAG, MPI_COMM_WORLD, "Error receiving work from task manager, aborting.", my_id);
}
