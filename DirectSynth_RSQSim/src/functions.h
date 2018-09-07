int master(struct sgtfileparams* sgtfilepar, MPI_Comm* sgt_handler_comm, int num_sgt_readers, char* stat, int run_id, int run_PSA, int run_rotd, int run_duration);
int sgt_handler(struct sgtfileparams* sgtfilepar, int num_comps, MPI_Comm* sgt_handler_comm, int num_sgt_handlers, MPI_Comm* sgt_readers_comm, int my_id);
int task_manager(int num_sgt_handlers, int num_workers, int num_procs, long long MAX_BUFFER_SIZE, int my_id);
int worker(int argc, char** argv, int num_sgt_handlers, struct sgtfileparams* sgtfilepar, char stat[64], float slat, float slon, int run_id, float det_max_freq, float stoch_max_freq, int run_PSA, int run_rotd, int run_duration, int my_id);

//misc.c
void init_sgtfileparams(struct sgtfileparams* sfp);
void construct_sgtmast_datatype(MPI_Datatype* sgtmaster_type);
void construct_sgtindx_datatype(MPI_Datatype* sgtindex_type);
void construct_sgtheader_datatype(MPI_Datatype* sgtheader_type);
void construct_sgtdata_datatype(MPI_Datatype* sgtdata_type, int nt);
void construct_worker_task_datatype(MPI_Datatype* worker_task_type);
void construct_worker_message_datatype(MPI_Datatype* worker_message_type);
void *check_realloc(void *ptr,size_t len);

//mpi_functions.c
void getMPIInfo(int *my_id, int* num_procs);
void constructSGTHandlerComm(int sgt_readers, MPI_Comm* sgt_handler_comm);
void constructSGTReaderComm(int sgt_readers, MPI_Comm* sgt_readers_comm);
void check_bcast(void* buf, int num_items, MPI_Datatype type, int root, MPI_Comm comm, char* error_msg, int my_id);
void check_send(void* buf, int num_items, MPI_Datatype type, int dest, int tag, MPI_Comm comm, char* error_msg, int my_id);
void check_recv(void* buf, int num_items, MPI_Datatype type, int src, int tag, MPI_Comm comm, char* error_msg, int my_id);

//fp_cache.c
void create_fp_cache();
FILE* find_and_use_fp(char* filename);
int find_fp(char* filename);
void remove_fp_cache();

//iofunc.c
void *check_malloc(size_t len);
void *check_realloc(void *ptr,size_t len);
void check_free(void* ptr);
FILE *fopfile(char *name,char *mode);
void fopfile_ro(char* name, FILE** fp);
int frite(FILE* fp, void *pntr, size_t record_size, size_t num_records);
int freed(FILE* fp, void* pntr, size_t record_size, size_t num_records);
int parse_rup_geom(char* rup_geom_file, struct rup_geom_point** rg_points);

//log.c
void open_log(int my_id);
void write_log(char* string);
void close_log();

//synth.c
int run_synth(task_info* t_info, int* proc_points, int num_sgt_handlers, char stat[64], float slat, float slon, int run_id, float det_max_freq, float stoch_max_freq, int run_PSA, int run_rotd, int run_duration, int my_id);
void request_sgt(struct sgtheader* sgthead, float* sgtbuf, int num_comps, int request_from_handler_id, long long** sgts_by_handler, int starting_index, int ending_index, int nt, int my_id);
void send_data_cluster(char data_filename[256], int src_id, int rup_id, int start_rv, int end_rv, void* data, int data_size_bytes, int my_id);

//jbsim3d.c
#include "srf_structure.h"
float** jbsim3d_synth(float*** seis_return, struct seisheader* header, char stat[], float slon, float slat, int ntout, char seis_file[], struct standrupformat* srf, struct sgtfileparams* sfp, struct sgtparams* sgtparms, struct sgtmaster sgtmast, struct sgtindex* sgtindx, struct geoprojection geop, long long* indx_master, int num_sgts, long long** sgts_by_handler, int* num_sgts_by_handler, int num_sgt_handlers, int num_rup_vars, struct rupture_variation* rup_vars, int my_id);

//misc_subs.c
double sfrand(int *seed);
void set_ne(float *elon,float *elat,float *slon,float *slat,float *sn,float *se);
void get_master_list_opt(struct sgtparams *sgtp, int np, long long* mindx, int* nm);
void get_indx(float *lon,float *lat,float *dep,struct sgtindex *indx,struct geoprojection *gp);
double nt_tol(float fnt,int gnt);

//sgt3d_subs.c
void get_sgtpars(struct sgtfileparams *sgtfpar,struct sgtmaster *sgtmast,struct sgtindex **sgtindx);
void find_sgt(struct sgtparams *sgtpar,struct sgtmaster *sgtmast,struct sgtindex *sgtindx,struct sgtindex *eqindx,struct sgtindex *statindx,float *maxd,float *fwt);

//geoproj_subs.c
void set_geoproj(struct sgtmaster *sgtmast,struct geoprojection *geop);

//rotd.c
int rotd(struct seisheader* header, float* seis_data, struct rotD_record* rotD_records);

//duration.c
#include "duration.h"
int duration(struct seisheader* header, float* full_seis, struct duration_record* out);

//integ_diff.c
void integ_diff(int integ_notdiff, float* seis, int nt, float dt);

extern void spectrad_(struct seisheader* header, int* nx, int* ny, int* npts, float* dt, char* seis_units, char* output_units, char* output_type, char* period, float* filter_high_hz, char* byteswap, char* input_file, char* output_file, float* seis, int* output_option, float* psa_data, int seis_units_len, int output_units_len, int output_type_len, int period_len, int byteswap_len, int input_file_len, int output_file_len);
