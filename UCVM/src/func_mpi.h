void mpi_init(int *ac,char ***av,int *np,int *id,char *pname,int *len);
void mpi_exit(int);
void check_procname(char *name,int *len);
void mpi_final(char *s);
void mpi_sndrcv2(void *,void *,void *,void *,int,int,int,int);
void mpi_global_val(void *,void *,char *,int,MPI_Datatype,MPI_Op);
