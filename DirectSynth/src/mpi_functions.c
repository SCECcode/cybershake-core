/*
 * mpi_functions.c
 *
 *  Created on: Dec 12, 2014
 *      Author: scott
 */
#include "include.h"
#include "structure.h"
#include "defs.h"
#include "functions.h"

void getMPIInfo(int *my_id, int* num_procs) {
	int rc = MPI_Comm_rank(MPI_COMM_WORLD, my_id);
	if (rc!=MPI_SUCCESS) {
		fprintf(stderr, "Error computing rank.");
		MPI_Finalize();
		exit(1);
	}
	rc = MPI_Comm_size(MPI_COMM_WORLD, num_procs);
	if (rc!=MPI_SUCCESS) {
		fprintf(stderr, "%d) Error computing size.", *my_id);
		MPI_Finalize();
		exit(2);
	}
}

void constructSGTHandlerComm(int sgt_handlers, MPI_Comm* sgt_handler_comm) {
	//Create a SGT handler communicator
	int* sgt_handler_ranks = check_malloc(sizeof(int)*sgt_handlers);
	int i;
	for (i=0; i<sgt_handlers; i++) {
		sgt_handler_ranks[i] = i;
	}
	MPI_Group world_group, sgt_handler_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	MPI_Group_incl(world_group, sgt_handlers, sgt_handler_ranks, &sgt_handler_group);
	MPI_Comm_create(MPI_COMM_WORLD, sgt_handler_group, sgt_handler_comm);
	free(sgt_handler_ranks);
}

void constructSGTReaderComm(int sgt_handlers, MPI_Comm* sgt_readers_comm) {
	//Constructs a communicator to do the SGT reading (ranks 1 - sgt_handlers-1)
	int* sgt_reader_ranks = check_malloc(sizeof(int)*(sgt_handlers-1));
	int i;
	for (i=0; i<sgt_handlers-1; i++) {
		sgt_reader_ranks[i] = i+1;
	}
	MPI_Group world_group, reader_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	MPI_Group_incl(world_group, sgt_handlers-1, sgt_reader_ranks, &reader_group);
	MPI_Comm_create(MPI_COMM_WORLD, reader_group, sgt_readers_comm);
	free(sgt_reader_ranks);
}


void check_bcast(void* buf, int num_items, MPI_Datatype type, int root, MPI_Comm comm, char* error_msg, int my_id) {
    /*if (my_id==root) {
	printf("%d) sending %d items from location %ld.\n", my_id, num_items, buf);
    } else {
        printf("%d) receiving %d items to location %ld from %d.\n", my_id, num_items, buf, root);
    }
    fflush(stdout);*/
    if (my_id==root) {
	if (debug) {
		char buf[256];
		sprintf(buf, "Broadcasting from process %d.", my_id);
		write_log(buf);
	}
    }
    int error = MPI_Bcast(buf, num_items, type, root, comm);
    if (error!=MPI_SUCCESS) {
    	fprintf(stderr, "MPI bcast error.\n");
	fflush(stderr);
	char string[256];
	int err_len;
	MPI_Error_string(error, string, &err_len);
        fprintf(stderr, "Process %d: %s\n", my_id, error_msg);
        fprintf(stderr, "Error message: %s\n", string);
        MPI_Finalize();
        exit(1);
    }
}

void check_send(void* buf, int num_items, MPI_Datatype type, int dest, int tag, MPI_Comm comm, char* error_msg, int my_id) {
	if (debug) {
		char buf[256];
		sprintf(buf, "Sending message to processor %d.", dest);
		write_log(buf);
	}
	int error = MPI_Send(buf, num_items, type, dest, tag, comm);
    if (error!=MPI_SUCCESS) {
    	fprintf(stderr, "MPI send error.\n");
		fflush(stderr);
		char string[256];
		int err_len;
	    MPI_Error_string(error, string, &err_len);
        fprintf(stderr, "From process %d to %d: %s\n", my_id, dest, error_msg);
        fprintf(stderr, "Error message: %s\n", string);
        MPI_Finalize();
        exit(2);
	}
}

void check_recv(void* buf, int num_items, MPI_Datatype type, int src, int tag, MPI_Comm comm, char* error_msg, int my_id) {
	MPI_Status status;
	int error = MPI_Recv(buf, num_items, type, src, tag, comm, &status);
    if (error!=MPI_SUCCESS) {
    	fprintf(stderr, "MPI recv error.\n");
		fflush(stderr);
		char string[256];
		int err_len;
	    MPI_Error_string(error, string, &err_len);
        fprintf(stderr, "From process %d to %d: %s\n", src, my_id, error_msg);
        fprintf(stderr, "Error message: %s\n", string);
        MPI_Finalize();
        exit(3);
	}
}
