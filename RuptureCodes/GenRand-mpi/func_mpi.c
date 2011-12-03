#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "func_mpi.h"
#include "defs.h"

void mpi_exit(int val)
{
  MPI_Finalize();
  exit(val);
}


void mpi_barrier()
{
  MPI_Barrier(MPI_COMM_WORLD);
  return;
}


void mpi_init(int *ac,char ***av,int *np,int *id,char *pname,int *len)
{
  MPI_Init(ac,av);
  MPI_Comm_size(MPI_COMM_WORLD,np);
  MPI_Comm_rank(MPI_COMM_WORLD,id);

  MPI_Get_processor_name(pname,len);
}


void mpi_final(char *s)
{
  //fprintf(stderr,"%s\n",s);
  MPI_Finalize();
}


int mpi_file_open(MPI_File *fh, char *filename, MPI_Offset offset)
{
  if (MPI_File_open(MPI_COMM_WORLD, filename,
		    MPI_MODE_CREATE | MPI_MODE_WRONLY, 
		    MPI_INFO_NULL, fh) != MPI_SUCCESS) {
    fprintf(stderr, "Error opening file %s\n", filename);
    return(1);
  }
  if (offset > 0) {
    if (MPI_File_seek(*fh, offset, MPI_SEEK_SET) != MPI_SUCCESS) {
      fprintf(stderr, "Error seeking in file %s\n", filename);
      return(1);
    }
  }
  return(0);
}


void mpi_file_close(MPI_File *fh)
{
  MPI_File_close(fh);
  return;
}


int mpi_file_write(MPI_File *fh, void *buf, int count, 
		   int num_fields, MPI_Datatype *dt)
{
  MPI_Status status;
  int num_wrote;

  if (MPI_File_write(*fh, buf, count, *dt, &status) != MPI_SUCCESS) {
    fprintf(stderr, "Error writing to file\n");
    return(1);
  }

  MPI_Get_count(&status, MPI_INT, &num_wrote);
  if (num_wrote != count * num_fields) {
    fprintf(stderr, "Error writing output, wrote %d of %d\n", 
	    num_wrote, count);
    return(1);
  }

  return(0);
}


int mpi_file_write_at(MPI_File *fh, MPI_Offset offset,
		      void *buf, int count, int num_fields, 
		      MPI_Datatype *dt)
{
  MPI_Status status;
  int num_wrote;

  /* Disable collective IO if directed */
#ifndef RUPGEN_ENABLE_MPI_COLL_IO
  if (MPI_File_write_at(*fh, offset, buf, count, *dt, 
  			&status) != MPI_SUCCESS) {
    fprintf(stderr, "Error writing to file\n");
    return(1);
  }
#else
  if (MPI_File_write_at_all(*fh, offset, buf, count, *dt, 
  			    &status) != MPI_SUCCESS) {
    fprintf(stderr, "Error writing to file\n");
    return(1);
  }
#endif
  
  MPI_Get_count(&status, MPI_INT, &num_wrote);
  if (num_wrote != count * num_fields) {
    fprintf(stderr, "Error writing output, wrote %d of %d\n", 
	    num_wrote, count);
    return(1);
  }

  return(0);
}

int mpi_register_rupinfo(MPI_Datatype *MPI_RUP_T, int *num_fields)
{
  // Register new mesh data type for c*MAX_FILENAME,i,f,i,i,i,i
  *num_fields = 7;
  MPI_Datatype dtype[7] = { MPI_CHAR, MPI_INT, MPI_FLOAT, MPI_INT, MPI_INT,
                                MPI_INT, MPI_INT };
  int blocklen[7] = { MAX_FILENAME, 1, 1, 1, 1, 1, 1 };
  MPI_Aint disp[7] = { 0, MAX_FILENAME, MAX_FILENAME+sizeof(int), 
		       MAX_FILENAME+(sizeof(int))+(sizeof(float)), 
		       MAX_FILENAME+(sizeof(int)*2)+(sizeof(float)), 
		       MAX_FILENAME+(sizeof(int)*3)+(sizeof(float)), 
		       MAX_FILENAME+(sizeof(int)*4)+(sizeof(float)) };
  if (MPI_Type_struct(*num_fields, blocklen, disp, dtype, MPI_RUP_T) != MPI_SUCCESS) {
	return(1);
  }
  if (MPI_Type_commit(MPI_RUP_T) != MPI_SUCCESS) {
  	return(1);
  }
  return(0);
}

