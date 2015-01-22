int master(struct sgtfileparams* sgtfilepar, MPI_Comm* sgt_handler_comm, int num_sgt_readers);
int sgt_handler(struct sgtfileparams* sgtfilepar, int num_comps, MPI_Comm* sgt_handler_comm, int num_sgt_handlers, int my_id);
int task_manager(int num_sgt_handlers, int num_workers, int num_procs, long long MAX_BUFFER_SIZE, int rup_var_spacing, int my_id);
int worker(int num_sgt_handlers, struct sgtfileparams* sgtfilepar, char stat[64], float slat, float slon, int run_id, float det_max_freq, float stoch_max_freq, int run_PSA, int run_rotd, int my_id);

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

//iofunc.c
void *check_malloc(size_t len);
void *check_realloc(void *ptr,size_t len);

//synth.c
void request_sgt(struct sgtheader* sgthead, float* sgtbuf, int request_from_handler_id, long long** sgts_by_handler, int starting_index, int ending_index, int nt, int my_id);

//jbsim3d.c
float** jbsim3d_synth(float*** seis_return, struct seisheader* header, char stat[], float slon, float slat, int ntout, char seis_file[], char rup_geom_file[], struct sgtfileparams* sfp, struct sgtparams* sgtparms, struct sgtmaster sgtmast, struct sgtindex* sgtindx, struct geoprojection geop, long long* indx_master, int num_sgts, long long** sgts_by_handler, int* num_sgts_by_handler, int num_sgt_handlers, int num_rup_vars, struct rupture_variation* rup_vars, int my_id);

//misc_subs.c
double sfrand(int *seed);
void set_ne(float *elon,float *elat,float *slon,float *slat,float *sn,float *se);
void get_master_list_opt(struct sgtparams *sgtp, int np, long long* mindx, int* nm);
