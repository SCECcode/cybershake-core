/*
 * master.c
 *
 *  Created on: Dec 11, 2014
 *      Author: scott
 */

#include "include.h"
#include "structure.h"
#include "defs.h"
#include "functions.h"

/*
 * 	Read in SGT header info
	Broadcast sgtmast to processes less than N
	Broadcast sgtindx to processes less than N
	Determine which processes are responsible for which points
	Broadcast process/point mapping to all processes
	Send sgtheader data to appropriate process < N
	Listen for messages:
		if output:
			write to disk
		if work complete from process N:
			exit
 */

int master(struct sgtfileparams* sgtfilepar, MPI_Comm* sgt_handler_comm, int num_sgt_readers) {
	//Read in SGT header info
	struct sgtmaster sgtmast;
	struct sgtindex sgtindx;
	get_sgtpars(&sgtfilepar,&sgtmast,&sgtindx);

	//Construct MPI datatype
	MPI_Datatype sgtmast_type, sgtindx_type;
	construct_sgtmast_datatype(&sgtmast_type);
	construct_sgtindx_datatype(&sgtindx_type);

	//Broadcast sgtmast, sgtindx to everyone - I think the workers need it too
	if (debug) write_log("Sending sgtmast, sgtindx to all.");
	check_bcast(sgtmast, 1, sgtmast_type, 0, MPI_COMM_WORLD, "Error broadcasting sgtmast, aborting.", 0);
	check_bcast(sgtindx, sgtmast.globnp, sgtindx_type, 0, MPI_COMM_WORLD, "Error broadcasting sgtindx, aborting.", 0);

	//Assign points and send out header info
	distribute_sgt_data(sgtfilepar, &sgtmast, &sgtindx, sgt_handler_comm, num_sgt_readers);

	//Listen for messages and respond accordingly
	master_listen();

	//Everything is done; at least, it better be...
	return 0;
}

void master_listen() {
	//If a process is sending output, write to disk
	//If process N has finished work, exit

	//Create file pointer cache, most recently used at top
	fp_cache_entry* fp_cache = NULL;

	master_msg msg;
	while (1) {
		if (debug) write_log("Listening for messages.");
		check_recv(&msg, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, "Error listening for message, aborting.", 0);
		if (msg.msg_type==DATA_FILE) {
			if (debug) {
				char buf[256];
				sprintf(buf, "Received request from %d to write data file.", msg.msg_src);
				write_log(buf);
			}
			//The worker is now going to send us a data file
			int data_size = msg.msg_data;
			int data_src = msg.msg_src;
			data_file df;
			df.data = check_malloc(sizeof(char)*data_size)
			if (debug) write_log("Receiving data file.");
			check_recv(&df, 256+data_size, MPI_BYTE, data_src, DATA_TAG, MPI_COMM_WORLD, "Error receiving data file, aborting.", 0);
			//Write data to file
			write_to_file(fp_cache, data_size, df);
			free(df.data);
		} else if (msg.msg_type==WORKERS_COMPLETE) {
			//close down file pointers and exit
			if (debug) write_log("Received workers complete message, shutting down.");
			remove_fp_cache(&fp_cache);
			break;
		} else {
			fprintf(stderr, "Master received a message with unknown type %d from id %d, aborting.", msg.msg_type, msg.msg_src);
			if (debug) close_log();
			MPI_Finalize()
			exit(4);
		}
	}

}

void write_to_file(fp_cache_entry* fp_cache, int data_size, data_file df) {
	//Use file pointer cache to avoid opening and closing files as much
	if (fp_cache==NULL) {
		create_fp_cache(&fp_cache);
	}
	//This will put our item in position 0
	find_and_use_fp(fp_cache, df.filename);
	if (debug) {
		char buf[256];
		sprintf(buf, "Writing to data file %s.", df.filename);
		write_log(buf);
	}
	frite(fp_cache[0].fp, 1, data_size, df.data);
}

void distribute_sgt_data(struct sgtfileparams* sgtfilepar, struct sgtmaster* sgtmast, struct sgtindex* sgtindx, MPI_Comm* sgt_handler_comm, int num_sgt_readers) {
	int* proc_points = check_malloc(sizeof(int)*(num_sgt_readers+1));

	//Assign points to processes
	assign_sgt_points(sgtmast, sgtindx, sgt_handler_comm, num_sgt_readers, proc_points);
	//Read and distribute header info to everyone
	send_header_data(sgtfilepar, sgtmast, sgt_handler_comm, num_sgt_readers, proc_points);

	free(proc_points);
}

void send_header_data(struct sgtfileparams* sgtfilepar, struct sgtmaster* sgtmast, MPI_Comm* sgt_handler_comm, int num_sgt_readers, int* proc_points) {
    MPI_Datatype sgtheader_type;
    construct_sgtheader_datatype(&sgtheader_type);

	struct sgtheader* buffer = NULL;
	//So we can iterate
	char* header_files[3];
	header_files[0] = sgtfilepar->xfile_header;
	header_files[1] = sgtfilepar->yfile_header;
	header_files[2] = sgtfilepar->zfile_header;
	FILE* header_fps[3];
	header_fps[0] = sgtfilepar->x_head_fp;
	header_fps[1] = sgtfilepar->y_head_fp;
	header_fps[2] = sgtfilepar->z_head_fp;
	int i, j;
	if (debug) write_log("Sending header data to SGT handlers.");
	for (i=0; i<3; i++) {
		if (header_files[i][0] != '\0') {
			fseek(header_fps[i], sizeof(struct sgtmaster) + sgtmast.globnp*sizeof(struct sgtindex), SEEK_SET);
			for (j=1; j<num_sgt_readers; j++) {
				buffer = check_realloc(buffer, sizeof(struct sgtheader)*(proc_points[j+1] - proc_points[j]));
				freed(header_fps[j], buffer, sizeof(struct sgtheader), proc_points[j+1] - proc_points[j]);
				check_send(buffer, proc_points[j+1] - proc_points[j], sgtheader_type, j, SGT_HEADER_TAG, sgt_handler_comm, "Error sending SGT header data, aborting.", 0);
			}
			fclose(header_fps[i]);
		}
	}
	free(buffer);
}

void assign_sgt_points(struct sgtmaster* sgtmast, struct sgtindex* sgtindx, MPI_Comm* sgt_handler_comm, int num_sgt_readers, int* proc_points) {
	//We will exclude processor 0 from managing points so it can just manage I/O
	//The indices in the proc points array correspond directly to the processor IDs, we just ignore 0
	int i;
	//We are not including 0 in the list of processes with points
	proc_points[0] = -1;
	proc_points[1] = 0;
	proc_points[num_sgt_readers] = sgtmast.globnp;
	float avg_points_per_proc = ((float)sgtmast.globnp)/((float)(num_sgt_readers-1));
	int cur_point_loc = 0;
	for (i=1; i<num_sgt_readers; i++) {
		if (i==num_sgt_readers-1) {
			proc_points[i+1] = sgtmast.globnp;
		} else {
			proc_points[i+1] = (int)(avg_points_per_proc * (i+1));
		}
	}
	//Just broadcast the points, since we want the workers to have the points-to-process mapping too
	if (debug) write_log("Broadcasting point-to-process mapping.");
	check_bcast(proc_points, num_sgt_readers+1, MPI_INT, 0, MPI_COMM_WORLD, "Error broadcasting SGT point-to-process mapping to all, aborting.", 0);
}
