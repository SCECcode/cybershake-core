#ifndef RUPGEN_FUNC_MPI_H
#define RUPGEN_FUNC_MPI_H

#include <mpi.h>

void mpi_exit(int);
void mpi_barrier();
void mpi_init(int *ac,char ***av,int *np,int *id,char *pname,int *len);
void mpi_final(char *s);

int mpi_file_open(MPI_File *fh, char *filename, MPI_Offset offset);
void mpi_file_close(MPI_File *fh);

int mpi_file_write(MPI_File *fh, void *buf, int count, 
		   int num_fields, MPI_Datatype *dt);
int mpi_file_write_at(MPI_File *fh, MPI_Offset offset, 
		      void *buf, int count, 
		      int num_fields, MPI_Datatype *dt);

int mpi_register_rupinfo(MPI_Datatype *MPI_RUP_T, int *num_fields);

#endif
