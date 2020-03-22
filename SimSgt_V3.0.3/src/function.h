/* Started 07/07/04- attempt to go towards ANSI C */

void mpi_init(int *ac,char ***av,int *np,int *id,char *pname,int *len);
void mpi_final(char *s);

void mpi_sndrcv2(void *,void *,void *,void *,int,int,int,int);
void mpi_global_val(void *,void *,char *,int,MPI_Datatype,MPI_Op);
void mpi_global_vel_minmax(struct vel_minmax *,struct nodeinfo *);
void mpi_exit(int);

void mpi_global_prof1d(struct media_input *,struct nodeinfo *);

void swap_in_place(int,char *);
void swap_in_place_d(int,char *);
void set_tmpdir(char *,int,char *,char *);

void check_procname(char *name,int *len);
void get_ny1ny2(int,int,int,int,int,int *,int *,int *);

void init_field_zero(float *,int);
void init_field_val(float,float *,int);
void floor_c(float *,float *,int);

void *check_malloc(size_t);
void *check_realloc(void *,size_t);

void resample(float *,int,float *,int,int,int,float *,float *,int,float *,float *);
void taper_norm(float *,float *,int,float *);
void norm(float *,float *,int);
double nt_tol(float,int);

void get_dc_par(struct pntsrcs *,float *,float *,int);
void get_pointmt_par(struct pntsrcs *,float *,float *,int);
void get_ff_par(struct pntsrcs *,float *,float *,int,float *);
void get_srf_par(struct pntsrcs *,struct runparams *,int,float *,float *,float *);

void init_outp(struct outputfields *,struct runparams *,char *,struct pntsrcs *,struct restartinfo *,struct dump_output_info *);
void check_rpars(char *,struct runparams *,struct runparams *);
void dump_restart(struct restartinfo *,char *,int,struct runparams *,struct modelstorage *,float *,int,int,int,int);
void get_restart(struct restartinfo *,struct modelstorage *,float *,int,int,int);
void print_rpars(struct runparams *,struct runparams *);
void dump_files(struct outputfields *,int,struct dump_output_info *);
void rm_dump_files(struct outputfields *,int,struct dump_output_info *);

void init_outpP3(struct outputfields *,struct runparamsP3 *,char *,struct pntsrcs *,struct restartinfoP3 *,struct dump_output_info *,struct runparams *rpOLD);
void check_rparsP3(char *,struct runparamsP3 *,struct runparamsP3 *);
void dump_restartP3(struct restartinfoP3 *,char *,int,struct runparamsP3 *,struct modelstorage *,float *,int);
void get_restartP3(struct restartinfoP3 *,float *,int,int,int);
void print_rparsP3(struct runparamsP3 *,struct runparamsP3 *);

float *reed_model(struct modelstorage *,float *,int);

int eff_media(struct media_input *,float *,int,int,int,int,int,struct modelstorage *,int);
void init_media2(struct media_input *,float *,int,int,int,int,int,int,float *,struct qvalues *,struct modelstorage *,int,int);
void init_media(struct media_input *,float *,struct runparams *,struct qvalues *,struct modelstorage *,int);

void reedmedslice(int,int,int,float *,int,int,int,int,struct qvalues *,int,int,float *,float *,int);

size_t reed(int,void *,size_t);
size_t rite(int,void *,size_t);
int opfile_ro(char *);
int opfile(char *);
int cropfile_rw(char *);
int croptrfile(char *);
int fcp_read_write(char *,char *);

void write_mt(struct pntsrcs *,float *,float *,int,int,int,int);

void init_modelc(struct runparams *);

void gcproj(float *,float *,float *,float *,float *,double *,double *,double *,double *,int);
void gen_matrices(double *,double *,float *,float *,float *);

void latlon2km(float *,float *,float *,float *,float *);
void set_g2(float *,float *);
void geocen(float *,double);

void check_nan(float *p,float *m,int nx,int nz,int iy,int it,int j);

void init_dc(struct pntsrcs *,float *);
void init_pointmt(struct pntsrcs *,float *,float *);
void init_mt(struct pntsrcs *,float *,struct modelstorage *,float *,int,int,int,int,int,char *,float *,char *,char *,int);

int readsource(char *,float *,float *,int);

int get_nproc(struct nodeinfo *);
void get_nodetype_neighbor(struct nodeinfo *);
void get_nodetype_neighbor_topol(struct nodeinfo *,char *);
void get_n1n2(int,struct nodeinfo *);
void get_n1n2_indv(int,int,int,int,int,int,int *,int *,int *);

int global2local_srcs(struct pntsrcs *,struct nodeinfo *,int);

void set_model_storage(struct modelstorage *,struct nodeinfo *,int,int,int,int *,int *,char *,char *);
void pad_medslice(float *,struct nodeinfo *,int);
void pad_medslice_ckvpvs(float *,struct nodeinfo *,int,float *);
void pad_medslice_ckvpvsP3_OLD(float *,struct nodeinfo *,int,struct media_input *);
void pad_medslice_ckvpvsP3_smooth1(float *,struct nodeinfo *,int,struct media_input *);
void pad_medslice_ckvpvsP3_smooth2(float *,struct nodeinfo *,int,struct media_input *);

size_t set_model_storageP3(struct modelstorage *,int,int,int);

void init_mediaP3(struct media_input *,float *,struct runparamsP3 *,struct qvalues *,int);
void genmod3dP3(int *nprof,char *,struct runparamsP3 *,struct medstencil *,struct medprof *,int);
void getmedsliceP3(struct medstencil *,struct medprof *,int,float *,struct runparamsP3 *,int,int,float *,float *);
int eff_mediaP3(struct media_input *,float *,struct runparamsP3 *);

void reedmedsliceP3_ms1(int,int,int,float *,int,struct qvalues *,int,float *,float *,int,int,struct nodeinfo *);
void reedmedsliceP3_ms2(struct media_input *,int,float *,int,int,int,struct nodeinfo *);
void reedmedsliceP3_ms3(int,int,int,float *,int,struct qvalues *,struct media_input *,int,int,int,struct nodeinfo *);

void get_1dprofile(struct media_input *,float *,int);

void mvar_coefsP3(float *,struct runparamsP3 *,float *,float *,float *,float *,int);

void set_vminvmax_reedmodP3(float *,struct fdcoefs *,int,int,int,float *);
void set_izord2P3(float *,struct fdcoefs *,int,int,int,float *);
void get_vel_minvmaxP3(float *,struct vel_minmax *,int,int,int,struct runparamsP3 *);
void get_vel_minmaxP3(float *,struct vel_minmax *,int,int,int,struct runparamsP3 *);

void tstepvP3(int,int,float **,float **,float *,struct fdcoefs *,int,int,struct interface *,struct nodeinfo *);
void tstepvbndP3(int,int,float **,float **,float *,struct fdcoefs *,int,int,struct interface *,struct nodeinfo *);
void abs_xzbndP3(int,int,float **,float *,float *,float **,int,int,struct interface *,struct nodeinfo *);
void abs_xzbnd1P3(int,int,float **,float *,float *,float **,int,struct interface *,struct nodeinfo *);
void abs_ybndP3(int,int,float **,float **,float *,struct fdcoefs *,int,int,struct interface *,struct nodeinfo *);

void tsteppP3(int,int,float **,float **,struct fdcoefs *,int,int,int,struct nodeinfo *);
void tsteppbndP3(int,int,float **,float **,struct fdcoefs *,int,int,int,struct nodeinfo *);

void get_dc_parP3(struct pntsrcs *,float *,float *);
void get_pointmt_parP3(struct pntsrcs *,float *,float *);
void get_bforce_parP3(struct pntsrcs *,float *,struct runparamsP3 *);
void get_srf_parP3(struct pntsrcs *,struct runparamsP3 *,int,float *,float *,float *);

void get_dc_parP3_ns(struct pntsrcs *,struct runparamsP3 *,float *,float *);
void get_pointmt_parP3_ns(struct pntsrcs *,struct runparamsP3 *,float *,float *);

void init_dcP3(struct pntsrcs *,struct runparamsP3 *);
void init_pointmtP3(struct pntsrcs *,struct runparamsP3 *);
void init_bforceP3(struct pntsrcs *,struct runparamsP3 *);
void init_srfP3(struct pntsrcs *,float *,char *,char *,char *,struct runparamsP3 *);

void init_dcP3_stf(struct pntsrcs *,struct runparamsP3 *,float *,int);
void init_pointmtP3_stf(struct pntsrcs *,struct runparamsP3 *,float *,int);
void get_MTavg_coefs(struct pntsrcs *,int,struct nodeinfo *);

size_t exchange_pvP3(float *,float *,struct nodeinfo *,int,int);
size_t exchange1_pvP3(float *,float *,struct nodeinfo *,int,int);
size_t exchange1_wptrs_pvP3(float *,float *,struct nodeinfo *,int,int,struct tms *,int);
size_t exchange1_nocpy_pvP3(float *,float *,struct nodeinfo *,int,int);

void get_adjsrc_parP3(struct pntsrcs *,struct runparamsP3 *);
void init_adjsrcP3(struct pntsrcs *,struct runparamsP3 *);

void fsurf_vel(int nx,int nz,float **pvptr,float **medp,struct fdcoefs *fdc,int bndflag,struct media_input *);

void init_report_2015(clock_t,clock_t *,struct tms *,struct tms *,struct tms *,struct tms *,clock_t *,float *,float *,int);
void print_report_2015(clock_t,clock_t *,struct tms *,struct tms *,struct tms *,struct tms *,clock_t *,float *,float *,int);
void dtimes_lite(struct tms,clock_t,struct tms *);

void setcoefs_pvc(float *pvfield,struct nodeinfo *ni,int fs,int order,struct fdcoefs *coefs,float *dt,float *h,float *vsmin);

void tstepvP3_pvc(int nx,int nz,float **pvptr,float **medptr,float *pbnd,struct fdcoefs *fdc,int fs,int ybndflag,struct interface *intfac,struct nodeinfo *ni);

void tstepvbndP3_pvc(int nx,int nz,float **pvptr,float **medptr,float *pbnd,struct fdcoefs *fdc,int fs,int ybndflag,struct interface *intfac,struct nodeinfo *ni);

void tsteppP3_pvc(int nx,int nz,float **pvptr,float **medptr,struct fdcoefs *fdc,int fs,int ybndflag,int eflag,struct nodeinfo *ni);

void tsteppbndP3_pvc(int nx,int nz,float **pvptr,float **medptr,struct fdcoefs *fdc,int fs,int ybndflag,int eflag,struct nodeinfo *ni);

void diff_pvc(int n,float *a,float *c0,float *c1,float *x,float *b,float *c,float *d,float *e,struct fdcoefs *fdc,int ord,int iflag);

void ord4_pvc(int n,float *a,float *x,float *b,float *c,float *d,float *e,float *c0,float *c1);

void diffx_pvc(int n,float *a1,float *a2,float *a3,float *c0,float *c1,float *x1,float *x2,float *b,float *c,float *d,float *e,struct fdcoefs *fdc,int ord,int iflag);

void ord4x_pvc(int n,float *a1,float *a2,float *a3,float *x1,float *x2,float *b,float *c,float *d,float *e,float *c0,float *c1);

void diffx_mv_pvc(int n,float *a1,float *a2,float *a3,float *d0,float *d1,float *c1,float *c2,float *c3,float *x1,float *x2,float *kp,float *ks,float *b,float *c,float *d,float *e,struct fdcoefs *fdc,int ord,int iflag);

void ord4x_mv_pvc(int n,float *a1,float *a2,float *a3,float *c1,float *c2,float *c3,float *x1,float *x2,float *kp,float *ks,float *b,float *c,float *d,float *e,float *d0,float *d1);

void diff_mv_pvc(int n,float *a,float *c0,float *c1,float *cp,float *x,float *ks,float *b,float *c,float *d,float *e,struct fdcoefs *fdc,int ord,int iflag);

void ord4_mv_pvc(int n,float *a,float *cp,float *x,float *ks,float *b,float *c,float *d,float *e,float *c0,float *c1);

void add_srcP3(float *stfp,int nt,float *dt,float *pvf,float *medf,int isrc,int nx,int nz,int it,struct pntsrcs *srcs,int iflag,struct nodeinfo *ni);
