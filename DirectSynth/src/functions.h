#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//misc.c
void init_sgtfileparams(struct sgtfileparams* sfp);
void construct_sgtmast_datatype(MPI_Datatype* sgtmaster_type);
void construct_sgtindx_datatype(MPI_Datatype* sgtindex_type);
void construct_sgtheader_datatype(MPI_Datatype* sgtheader_type);
void construct_sgtdata_datatype(MPI_Datatype* sgtdata_type, int nt);
void construct_worker_task_datatype(MPI_Datatype* worker_task_type);
void *check_realloc(void *ptr,size_t len);

//mpi_functions.c
void getMPIInfo(int *my_id, int* num_procs);
void constructSGTHandlerComm(int sgt_readers, MPI_Group* sgt_handler_group);
void check_bcast(void* buf, int num_items, MPI_Datatype type, int root, MPI_Comm comm, char* error_msg, int my_id);
void check_send(void* buf, int num_items, MPI_Datatype type, int dest, int tag, MPI_Comm comm, char* error_msg, int my_id);
void check_recv(void* buf, int num_items, MPI_Datatype type, int src, int tag, MPI_Comm comm, char* error_msg, int my_id);

//fp_cache.c
void create_fp_cache(fp_cache_entry** fp_cache);
void find_and_use_fp(fp_cache_entry* fp_cache, char* filename);
int find_fp(fp_cache_entry* fp_cache, char* filename);
void remove_fp_cache(fp_cache_entry** fp_cache);

#endif
