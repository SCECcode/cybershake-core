/*
 * master.c
 *
 *  Created on: Dec 11, 2014
 *      Author: scott
 */

#include "include.h"
#include "defs.h"
#include "structure.h"
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
void master_listen(int* task_tuples, int num_ruptures, char* site, int run_id, int run_PSA, int run_rotd, int run_duration, int run_period_duration, int run_vert_rsp);
void write_to_file(int data_size, data_file_metadata* df, char* data);
void distribute_sgt_data(struct sgtfileparams* sgtfilepar, struct sgtmaster* sgtmast, struct sgtindex* sgtindx, MPI_Comm* sgt_handler_comm, int num_sgt_readers);
void send_header_data(struct sgtfileparams* sgtfilepar, struct sgtmaster* sgtmast, MPI_Comm* sgt_handler_comm, int num_sgt_readers, int* proc_points);
void assign_sgt_points(struct sgtmaster* sgtmast, struct sgtindex* sgtindx, MPI_Comm* sgt_handler_comm, int num_sgt_readers, int* proc_points);


int master(struct sgtfileparams* sgtfilepar, MPI_Comm* sgt_handler_comm, int num_sgt_readers, char* stat, int run_id, int run_PSA, int run_rotd, int run_duration, int run_period_duration, int run_vert_rsp) {
	//Read in SGT header info
	struct sgtmaster sgtmast;
	struct sgtindex* sgtindx;
	if (debug) write_log("Master reading sgtpars.");
	get_sgtpars(sgtfilepar,&sgtmast,&sgtindx);

	//Construct MPI datatype
	MPI_Datatype sgtmast_type, sgtindx_type;
	construct_sgtmast_datatype(&sgtmast_type);
	construct_sgtindx_datatype(&sgtindx_type);

	//Broadcast sgtmast, sgtindx to everyone - I think the workers need it too
	if (debug) write_log("Sending sgtmast to all.");
	check_bcast(&sgtmast, 1, sgtmast_type, 0, MPI_COMM_WORLD, "Error broadcasting sgtmast, aborting.", 0);
	if (debug) write_log("Sending sgtindx to all.");
	check_bcast(sgtindx, sgtmast.globnp, sgtindx_type, 0, MPI_COMM_WORLD, "Error broadcasting sgtindx, aborting.", 0);
	
	//Assign points and send out header info
	distribute_sgt_data(sgtfilepar, &sgtmast, sgtindx, sgt_handler_comm, num_sgt_readers);

	//Get list of (source id, rupture id, #rup vars) from task manager so we can checkpoint
	int num_ruptures;
	check_recv(&num_ruptures, 1, MPI_INT, num_sgt_readers, NUM_RUPTURES_TAG, MPI_COMM_WORLD, "Error receiving num ruptures from task manager, aborting.", 0);
	int* task_tuples = check_malloc(sizeof(int)*3*num_ruptures);
	check_recv(task_tuples, num_ruptures*3, MPI_INT, num_sgt_readers, VARIATION_INFO_TAG, MPI_COMM_WORLD, "Error receiving variation info tuples from task manager, aborting.", 0);

	//Listen for messages and respond accordingly
	master_listen(task_tuples, num_ruptures, stat, run_id, run_PSA, run_rotd, run_duration, run_period_duration, run_vert_rsp);

	free(task_tuples);
	//Everything is done; at least, it better be...
	return 0;
}

void master_listen(int* task_tuples, int num_ruptures, char* site, int run_id, int run_PSA, int run_rotd, int run_duration, int run_period_duration, int run_vert_rsp) {
	//Construct easy-access task tuples for monitoring status for restart
	short src_rup_table[MAX_SOURCE_ID][MAX_RUPTURE_ID];
	int i;
	if (debug) {
		char buf[256];
		sprintf(buf, "%d ruptures to be processed.", num_ruptures);
		write_log(buf);
	}
	for (i=0; i<num_ruptures; i++) {
		int src = task_tuples[3*i];
		int rup = task_tuples[3*i+1];
		int num_files = (1+run_PSA+run_rotd+run_duration+run_period_duration+run_vert_rsp)*(task_tuples[3*i+2]);
		if (debug) {
			char buf[256];
			sprintf(buf, "Source %d, rupture %d has %d files.", src, rup, num_files);
			write_log(buf);
		}
		if (src>=MAX_SOURCE_ID || rup>=MAX_RUPTURE_ID) {
			fprintf(stderr, "Source ID %d, rupture ID %d is in the input file, but MAX_SOURCE_ID=%d and MAX_RUPTURE_ID=%d.  Change these values in defs.h.\n", src, rup, MAX_SOURCE_ID, MAX_RUPTURE_ID);
                        if (debug) close_log();
                        MPI_Finalize();
                        exit(5);
		}

		src_rup_table[src][rup] = (short)(num_files);
		//If the output files already exist, clear them
		char filename[256];
		sprintf(filename, "Seismogram_%s_%d_%d_%d.grm", site, run_id, src, rup);
		unlink(filename);
		if (run_PSA) {
			sprintf(filename, "PeakVals_%s_%d_%d_%d.bsa", site, run_id, src, rup);
			unlink(filename);
		}
		if (run_rotd) {
			sprintf(filename, "RotD_%s_%d_%d_%d.rotd", site, run_id, src, rup);
			unlink(filename);
		}
		if (run_duration) {
			sprintf(filename, "Duration_%s_%d_%d_%d.dur", site, run_id, src, rup);
			unlink(filename);
		}
		if (run_period_duration) {
			sprintf(filename, "PeriodDuration_%s_%d_%d_%d.dur", site, run_id, src, rup);
			unlink(filename);
		}
		if (run_vert_rsp) {
			sprintf(filename, "VerticalRSP_%s_%d_%d_%d.rsp", site, run_id, src, rup);
			unlink(filename);
		}
	}
	//Open checkpoint file
	FILE* checkpoint_fp = fopfile(CHECKPOINT_FILE, "a");

	//If a process is sending output, write to disk
	//If process N has finished work, exit

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
			data_file_metadata df;
			if (debug) write_log("Receiving data file metadata.");
			check_recv(&df, 4*sizeof(int) + 256, MPI_BYTE, data_src, DATA_FILENAME_TAG, MPI_COMM_WORLD, "Error receiving data file metadata, aborting.", 0);
			char* data = check_malloc(data_size);
			if (debug) write_log("Receiving data contents.");
			check_recv(data, data_size, MPI_BYTE, data_src, DATA_TAG, MPI_COMM_WORLD, "Error receiving data contents, aborting.", 0);
			//Write data to file
			write_to_file(data_size, &df, data);
			free(data);
			//Update src_rup_table
			src_rup_table[df.src_id][df.rup_id] -= df.ending_rv - df.starting_rv;
			if (src_rup_table[df.src_id][df.rup_id]==0) {
				//This rupture is done; sync the files and write to checkpoint log
				if (debug) {
					char buf[256];
					sprintf(buf, "Source %d, rupture %d finished.  fsyncing files and writing to checkpoint log.", df.src_id, df.rup_id);
					write_log(buf);
				}
				char filename[256];
				sprintf(filename, "Seismogram_%s_%d_%d_%d.grm", site, run_id, df.src_id, df.rup_id);
				fsync_and_close(filename);
				if (run_PSA) {
					sprintf(filename, "PeakVals_%s_%d_%d_%d.bsa", site, run_id, df.src_id, df.rup_id);
					fsync_and_close(filename);
				}
				if (run_rotd) {
					sprintf(filename, "RotD_%s_%d_%d_%d.rotd", site, run_id, df.src_id, df.rup_id);
					fsync_and_close(filename);
				}
				if (run_duration) {
					sprintf(filename, "Duration_%s_%d_%d_%d.dur", site, run_id, df.src_id, df.rup_id);
					fsync_and_close(filename);
				}
				if (run_period_duration) {
					sprintf(filename, "PeriodDuration_%s_%d_%d_%d.dur", site, run_id, df.src_id, df.rup_id);
                    fsync_and_close(filename);
				}
				fprintf(checkpoint_fp, "%d %d\n", df.src_id, df.rup_id);
				fflush(checkpoint_fp);
			}
		} else if (msg.msg_type==WORKERS_COMPLETE) {
			//close down file pointers and exit
			if (debug) write_log("Received workers complete message, shutting down.");
			remove_fp_cache();
			break;
		} else {
			fprintf(stderr, "Master received a message with unknown type %d from id %d, aborting.", msg.msg_type, msg.msg_src);
			if (debug) close_log();
			MPI_Finalize();
			exit(4);
		}
	}

}

void write_to_file(int data_size, data_file_metadata* df, char* data) {
	//Use file pointer cache to avoid opening and closing files as much
	//This will put our item in position 0
	FILE* out_fp = find_and_use_fp(df->filename);
	if (debug) {
		char buf[256];
		sprintf(buf, "Writing %d bytes to data file %s.", data_size, df->filename);
		write_log(buf);
	}
	frite(out_fp, data, 1, data_size);
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
			fseek(header_fps[i], sizeof(struct sgtmaster) + sgtmast->globnp*sizeof(struct sgtindex), SEEK_SET);
			for (j=1; j<num_sgt_readers; j++) {
				buffer = check_realloc(buffer, sizeof(struct sgtheader)*(proc_points[j+1] - proc_points[j]));
				freed(header_fps[i], buffer, sizeof(struct sgtheader), proc_points[j+1] - proc_points[j]);
				if (debug) {
					char buf[256];
					sprintf(buf, "Master sending header info to handler %d.", j);
					write_log(buf);
				}
				check_send(buffer, proc_points[j+1] - proc_points[j], sgtheader_type, j, SGT_HEADER_TAG, *sgt_handler_comm, "Error sending SGT header data, aborting.", 0);
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
	proc_points[num_sgt_readers] = sgtmast->globnp;
	//We have the option to use an input file which specifies the # of points each handler is responsible for
	char point_mapping_file[256];
	point_mapping_file[0]='\0';
	getpar("pt_mapping_file","s",point_mapping_file);
	if (point_mapping_file[0]!='\0') {
		if (debug) write_log("Using point mapping file.");
		FILE* fp_in;
		char line[256];
		int id, num_pts, hits;
		fopfile_ro(point_mapping_file, &fp_in);
		//This file has a header, then
		//Handler_ID #SGT_pts #Hits
		fgets(line, 256, fp_in);
		for(i=1; i<num_sgt_readers; i++) {
			fscanf(fp_in, "%d %d %d\n", &id, &num_pts, &hits);
			proc_points[i+1] = proc_points[i]+num_pts;
		}
		fclose(fp_in);
	} else {
		float avg_points_per_proc = ((float)sgtmast->globnp)/((float)(num_sgt_readers-1));
		int cur_point_loc = 0;
		for (i=1; i<num_sgt_readers; i++) {
			if (i==num_sgt_readers-1) {
				proc_points[i+1] = sgtmast->globnp;
			} else {
				proc_points[i+1] = (int)(avg_points_per_proc * (i));
			}
		}
	}
	//Just broadcast the points, since we want the workers to have the points-to-process mapping too
	if (debug) write_log("Broadcasting point-to-process mapping.");
	check_bcast(proc_points, num_sgt_readers+1, MPI_INT, 0, MPI_COMM_WORLD, "Error broadcasting SGT point-to-process mapping to all, aborting.", 0);
}
