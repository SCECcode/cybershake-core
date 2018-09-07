/*
 * sgt_handler.c
 *
 *  Created on: Dec 11, 2014
 *      Author: scott
 */

#include "include.h"
#include "defs.h"
#include "structure.h"
#include "functions.h"

#include "cfuhash.h"

/*
 * if my_id>0 and my_id<N:
	Obtain sgtmast
	Obtain sgtindx
	Obtain process/point mapping
	Obtain appropriate sgtheader data
	Read SGT data with MPI I/O
	Listen for messages:
		if request for data:
			service request for SGT data
		if work complete from process N:
			exit
 */

void handler_listen(struct sgtmaster* sgtmast, struct sgtindex* sgtindx, struct sgtheader* sgthead[3], char** sgtdata, int num_my_points, int my_offset, int num_comps, int my_id);
void handle_SGT_request(long long* sgts_requested, int num_sgts_requested, cfuhash_table_t* hash_table, struct sgtmaster* sgtmast, struct sgtheader* sgthead[3], char** sgtdata, int num_my_points, int my_id, int num_comps, int dest, MPI_Datatype sgtheader_type, MPI_Datatype sgtdata_type);
void get_SGT_data(struct sgtfileparams* sgtfilepar, struct sgtmaster* sgtmast, int num_comps, int* proc_points, struct sgtheader* (*sgthead)[3], char*** sgtdata, int my_id, MPI_Comm* sgt_handler_comm, MPI_Comm* sgt_readers_comm);
int readSGT_MPI(char *sgt_fname, int nPoints, int nMyPoints, int nt, int indexMyPoint, char *sgt_buf, int separate_header_flag, MPI_Comm* sgt_handler_comm, int my_id);

int sgt_handler(struct sgtfileparams* sgtfilepar, int num_comps, MPI_Comm* sgt_handler_comm, int num_sgt_handlers, MPI_Comm* sgt_readers_comm, int my_id) {
	struct sgtmaster sgtmast;
	struct sgtindex* sgtindx;

	//Construct MPI datatype
	MPI_Datatype sgtmast_type, sgtindx_type;
	construct_sgtmast_datatype(&sgtmast_type);
	construct_sgtindx_datatype(&sgtindx_type);

	//Get sgtmast, sgtindx info
	if (debug) write_log("Receiving sgtmast from master.");
	check_bcast(&sgtmast, 1, sgtmast_type, 0, MPI_COMM_WORLD, "Error receiving sgtmast, aborting.", my_id);
	sgtindx = check_malloc(sizeof(struct sgtindex) * sgtmast.globnp);
        if (debug) write_log("Receiving sgtindx from master.");
	check_bcast(sgtindx, sgtmast.globnp, sgtindx_type, 0, MPI_COMM_WORLD, "Error receiving sgtindx, aborting.", my_id);

	//Get point-to-process mapping
	int* proc_points = check_malloc(sizeof(int)*(num_sgt_handlers+1));
	if (debug) write_log("Receiving point-to-process mapping from master.");
	check_bcast(proc_points, num_sgt_handlers+1, MPI_INT, 0, MPI_COMM_WORLD, "Error receiving SGT point-to-process mapping, aborting.", my_id);

	//Get SGT data
	struct sgtheader *sgthead[3];
	char** sgtdata;
	int num_my_points = proc_points[my_id+1] - proc_points[my_id];
	int my_offset = proc_points[my_id];
	get_SGT_data(sgtfilepar, &sgtmast, num_comps, proc_points, &sgthead, &sgtdata, my_id, sgt_handler_comm, sgt_readers_comm);

	//Listen for messages
	handler_listen(&sgtmast, sgtindx, sgthead, sgtdata, num_my_points, my_offset, num_comps, my_id);

	int i;
	for (i=0; i<3; i++) {
		free(sgtdata[i]);
		free(sgthead[i]);
	}
	free(sgtdata);
	free(sgtindx);
	free(proc_points);
	MPI_Type_free(&sgtmast_type);
	MPI_Type_free(&sgtindx_type);
	return 0;
}


void handler_listen(struct sgtmaster* sgtmast, struct sgtindex* sgtindx, struct sgtheader* sgthead[3], char** sgtdata, int num_my_points, int my_offset, int num_comps, int my_id) {
	/* Listen for messages:
			if request for data:
				service request for SGT data
			if work complete from process N:
				exit
	*/
	int i;

	//Use hashtable to keep mapping of SGT xyz to sgtdata index for quick access
	//Shouldn't have more than about 14k SGT points per process
	if (debug) write_log("Creating hash table of SGT points.");
	cfuhash_table_t* hash_table = cfuhash_new_with_initial_size(20000);
	cfuhash_set_flag(hash_table, CFUHASH_NO_LOCKING);
	//cfuhash doesn't copy the value, so we need dedicated mem for this
	int* value_array = check_malloc(sizeof(int)*num_my_points);
	int key_len = sizeof(long long);
	int val_len = sizeof(int);
	for (i=0; i<num_my_points; i++) {
		value_array[i] = i;
		cfuhash_put_data(hash_table, &(sgtindx[i+my_offset].indx), key_len, &(value_array[i]), val_len, NULL);
	}

	//Get SGT header data
    MPI_Datatype sgtheader_type;
    construct_sgtheader_datatype(&sgtheader_type);
    MPI_Datatype sgtdata_type;
    construct_sgtdata_datatype(&sgtdata_type, sgtmast->nt);

	handler_msg msg;
	long long* sgts_requested = NULL;
	while (1) {
		if (debug) write_log("Listening for messages.");
		check_recv(&msg, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, "Error listening for message, aborting.", my_id);
		if (msg.msg_type==SGT_REQUEST) {
			if (debug) {
				char buf[256];
				sprintf(buf, "Received request for %d SGTs from process %d.", msg.num_sgts_requested, msg.msg_src);
				write_log(buf);
			}
			//The worker is about to request a series of SGTs
			int num_sgts_requested = msg.num_sgts_requested;
			int data_src = msg.msg_src;
			//SGTs requested via series of long longs with the index
			sgts_requested = check_realloc(sgts_requested, sizeof(long long)*num_sgts_requested);
			if (debug) write_log("Receiving list of SGTs.");
			check_recv(sgts_requested, num_sgts_requested, MPI_LONG_LONG, data_src, SGT_REQUEST_TAG, MPI_COMM_WORLD, "Error receiving SGT request, aborting.", my_id);
			handle_SGT_request(sgts_requested, num_sgts_requested, hash_table, sgtmast, sgthead, sgtdata, num_my_points, my_id, num_comps, data_src, sgtheader_type, sgtdata_type);
		} else if (msg.msg_type==WORKERS_COMPLETE) {
			if (debug) write_log("Received workers complete message, shutting down.");
			break;
		} else {
			fprintf(stderr, "%d) SGT handler received message of unknown type %d from %d, aborting.", my_id, msg.msg_type, msg.msg_src);
			if (debug) close_log();
			MPI_Finalize();
			exit(4);
		}
	}
	//Clean up
	cfuhash_destroy(hash_table);
	free(value_array);
	free(sgts_requested);
	MPI_Type_free(&sgtheader_type);
	MPI_Type_free(&sgtdata_type);
}

void handle_SGT_request(long long* sgts_requested, int num_sgts_requested, cfuhash_table_t* hash_table, struct sgtmaster* sgtmast, struct sgtheader* sgthead[3],
		char** sgtdata, int num_my_points, int my_id, int num_comps, int dest, MPI_Datatype sgtheader_type, MPI_Datatype sgtdata_type) {
	//Send only the components which exist
	//Send all of X, then all of Y, then all of Z
	int i, j;
	//Send 1 message per component, reusing sending buffer to minimize size of buffer
	//Keep track of indices so we don't have to find for second component
	char* sending_buffer = check_malloc(sizeof(float)*(long long)sgtmast->nt*N_SGTvars*num_sgts_requested);
	int* sgts_to_send_index_list = check_malloc(sizeof(int)*num_sgts_requested);
	int key_len = sizeof(long long);
	size_t val_len;
	if (debug) write_log("Looking up SGTs in hash table.");
	for (i=0; i<num_sgts_requested; i++) {
		//Find this SGT in my data
		int* tmp;
		int found = cfuhash_get_data(hash_table, &(sgts_requested[i]), key_len, &tmp, &val_len);
		if (!found) {
			fprintf(stderr, "%d) Could not find SGT %ld in my hash table.  I am troubled.  Aborting since this is an error somewhere.", my_id, sgts_requested[i]);
			if (debug) close_log();
			MPI_Finalize();
			exit(4);
		}
		sgts_to_send_index_list[i] = *tmp;
	}
	//Send header info
        for (j=0; j<num_sgts_requested; j++) {
	        int sgt_index = sgts_to_send_index_list[j];
                memcpy(sending_buffer+j*sizeof(struct sgtheader), &(sgthead[0][sgt_index]), sizeof(struct sgtheader));
	}
        if (debug) {
        	char buf[256];
                sprintf(buf, "Sending SGT header data to process %d.", dest);
                write_log(buf);
	}
        check_send(sending_buffer, num_sgts_requested, sgtheader_type, dest, SGT_HEADER_DATA_TAG, MPI_COMM_WORLD, "Error sending SGT header data, aborting.", my_id);
	//Prepare buffer for sending
	//Don't send 0s for missing components; instead, let worker copy things appropriately
	for (i=0; i<num_comps; i++) {
		//Send raw data info
		for (j=0; j<num_sgts_requested; j++) {
			int sgt_index = sgts_to_send_index_list[j];
			memcpy(sending_buffer+(long long)j*sizeof(float)*N_SGTvars*sgtmast->nt, sgtdata[i]+(long long)sgt_index*N_SGTvars*sgtmast->nt*sizeof(float), sizeof(float)*N_SGTvars*sgtmast->nt);
		}
		if (debug) {
			char buf[256];
			sprintf(buf, "Sending raw SGT data, component %d, to process %d.", i, dest);
			write_log(buf);
		}
		check_send(sending_buffer, num_sgts_requested, sgtdata_type, dest, SGT_RAW_DATA_TAG, MPI_COMM_WORLD, "Error sending raw SGT data, aborting.", my_id);
	}
	if (debug) write_log("SGT request complete.");

	free(sending_buffer);
	free(sgts_to_send_index_list);
}


void get_SGT_data(struct sgtfileparams* sgtfilepar, struct sgtmaster* sgtmast, int num_comps, int* proc_points, struct sgtheader* (*sgthead)[3], char*** sgtdata,
		int my_id, MPI_Comm* sgt_handler_comm, MPI_Comm* sgt_readers_comm) {
	//Get SGT header data
    MPI_Datatype sgtheader_type;
    construct_sgtheader_datatype(&sgtheader_type);

	char* header_files[3];
	header_files[0] = sgtfilepar->xfile_header;
	header_files[1] = sgtfilepar->yfile_header;
	header_files[2] = sgtfilepar->zfile_header;

	int num_my_points = proc_points[my_id+1] - proc_points[my_id];

	if (debug) write_log("Receiving header info from master.");
	int i;
	for (i=0; i<3; i++) {
		(*sgthead)[i] = NULL;
		if (header_files[i][0] != '\0') {
			//printf("%d) receiving component %d header info.\n", my_id, i);
			(*sgthead)[i] = check_malloc(num_my_points*sizeof(struct sgtheader));
			//printf("%d) receiving %d points of header info into address %ld\n", my_id, num_my_points, (*sgthead)[i]);
			check_recv((*sgthead)[i], num_my_points, sgtheader_type, 0, SGT_HEADER_TAG, *sgt_handler_comm, "Error receiving SGT header, aborting.", my_id);
		}
	}
	//We are assuming that all the moments are the same; check to make sure that's a correct assumption, if not, abort
	if (num_comps>=2) {
		if ((*sgthead)[0][0].xmom!=(*sgthead)[1][0].ymom) {
			fprintf(stderr, "X-moment and Y-moment are not the same.  We assume this is true so we don't have to send the Y header separately, so we can't handle this and must abort.\n");
			if (debug) close_log();
			MPI_Finalize();
			exit(7);
		}	
	} 
	if (num_comps==3) {
        	if ((*sgthead)[0][0].xmom!=(*sgthead)[2][0].zmom) {
                        fprintf(stderr, "X-moment and Z-moment are not the same.  We assume all moments are equal so we don't have to send each header separately, so we can't handle this and must abort.\n"); 
                        if (debug) close_log();
                        MPI_Finalize();
                        exit(8);
                }
	}
	printf("%d) In get SGT data, sgthead[0] = %ld, sgthead[1] = %ld\n", my_id, (*sgthead)[0], (*sgthead)[1]);
	

		//Read SGT data
	char* sgt_files[3];
	sgt_files[0] = sgtfilepar->xfile;
	sgt_files[1] = sgtfilepar->yfile;
	sgt_files[2] = sgtfilepar->zfile;
	//Allocate memory for SGT storage
	*sgtdata = check_malloc(3*sizeof(char*));
	int separate_header_flag[3] = {0, 0, 0};
	for (i=0; i<3; i++) {
		(*sgtdata)[i] = NULL;
		if (header_files[i][0] != '\0') {
			separate_header_flag[i] = 1;
			(*sgtdata)[i] = check_malloc((long long)num_my_points*sizeof(float)*N_SGTvars*sgtmast->nt);
		} else {
			(*sgtdata)[i] = check_malloc((long long)num_my_points*(sizeof(struct sgtheader)+sizeof(float)*N_SGTvars*sgtmast->nt));
		}
	}


	for (i=0; i<3; i++) {
		if (sgt_files[i][0] != '\0') {
			if (debug) {
				char buf[256];
				sprintf(buf, "Reading from SGT file %s.", sgt_files[i]);
				write_log(buf);
			}
			readSGT_MPI(sgt_files[i], sgtmast->globnp, num_my_points, sgtmast->nt, proc_points[my_id], (*sgtdata)[i], separate_header_flag[0], sgt_readers_comm, my_id);
		}
	}
}

int readSGT_MPI(char *sgt_fname, int nPoints, int nMyPoints, int nt, int indexMyPoint, char *sgt_buf, int separate_header_flag, MPI_Comm* sgt_readers_comm, int my_id) {
  MPI_File fh;
  int err;
  int size;
  long sgt_pointDataSize;
  MPI_Datatype filetype;
  MPI_Status stat;
  MPI_Offset disp = 0;

  sgt_pointDataSize = sizeof(struct sgtheader) + nt*sizeof(float)*N_SGTvars;

  if (!separate_header_flag) {
	disp += sizeof(struct sgtmaster) + nPoints*(sizeof(struct sgtindex));
	sgt_pointDataSize = sizeof(struct sgtheader) + nt*sizeof(float)*N_SGTvars;
  } else {
	sgt_pointDataSize = nt*sizeof(float)*N_SGTvars;
  }
  disp += indexMyPoint*sgt_pointDataSize;

  //printf("%d) nMyPoints: %d, sgt_pointDataSize: %ld\n", my_id, nMyPoints, sgt_pointDataSize);
  //Check to make sure first argument to contiguous is less than INT_MAX
  if ((long)nMyPoints*(long)sgt_pointDataSize > INT_MAX) {
	fprintf(stderr, "%d) nMyPoints (%d) times sgt_pointDataSize (%d) is greater than INT_MAX, and will crash if passed to MPI_Type_contiguous.  Increase the number of sgt_handlers by at least %.1fx\n", my_id, nMyPoints, sgt_pointDataSize, ((float)nMyPoints*(float)sgt_pointDataSize)/INT_MAX);
	exit(2);
  }
  //MPI_Type_contiguous(nMyPoints*sgt_pointDataSize, MPI_BYTE, &filetype);
  MPI_Type_contiguous(nMyPoints*sgt_pointDataSize/sizeof(float), MPI_FLOAT, &filetype);
  MPI_Type_commit(&filetype);
  MPI_Type_size(filetype, &size);
  if(my_id==1) {
		if (debug) {
			char buf[256];
			sprintf(buf, "filetype size (supposedly nMyPoints*sgt_pointDataSize=%ld)=%d\n", nMyPoints*sgt_pointDataSize, size);
			write_log(buf);
		}
  }
  printf("%d) disp to read MPI=%ld\n", my_id, (long)disp);

  err = MPI_File_open(*sgt_readers_comm, sgt_fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  //err = MPI_File_set_view(fh, disp, MPI_BYTE, filetype, "native", MPI_INFO_NULL);
  err = MPI_File_set_view(fh, disp, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
  //printf("%d) Ready to read file\n", rank);
  //fflush(stdout);
  //err = MPI_File_read_at_all(fh, 0, sgt_buf, 1, filetype, &stat);
  err = MPI_File_read_all(fh, sgt_buf, 1, filetype, &stat);

  err = MPI_File_close(&fh);

  MPI_Type_free(&filetype);

  return 0;
}
