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

void constructSGTHandlerComm(int sgt_readers, MPI_Comm* sgt_handler_comm) {
	//Create a SGT handler communicator
	int* sgt_handler_ranks = check_malloc(sizeof(int)*sgt_readers);
	int i;
	for (i=0; i<sgt_readers; i++) {
		sgt_handler_ranks[i] = i;
	}
	MPI_Group world_group, sgt_handler_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	int check_rank;
	MPI_Group_rank(world_group, &check_rank);
	MPI_Group_incl(world_group, sgt_readers, sgt_handler_ranks, &sgt_handler_group);
	MPI_Comm_create(MPI_COMM_WORLD, sgt_handler_group, sgt_handler_comm);
	free(sgt_handler_ranks);
}

void check_bcast(void* buf, int num_items, MPI_Datatype type, int root, MPI_Comm comm, char* error_msg, int my_id) {
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
