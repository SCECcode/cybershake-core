/*
   memory_slices.c comtains the following functions:
      
      init_model_seek()
      set_model_seek()
      reed_model()
      rite_model()
*/

#include "include.h"

init_model_seek(mds,field_type)
struct modelstorage *mds;
int field_type;
{
if(mds->intmem)
   {
   if(field_type == WAVE_FIELD)
      mds->pv_bufptr = mds->pv_buf;
   else if(field_type == MEDIA_FIELD)
      mds->med_bufptr = mds->med_buf;
   }
else
   {
   if(field_type == WAVE_FIELD)
      {
      lseek(mds->pv_fdr,0,SEEK_SET);
      lseek(mds->pv_fdw,0,SEEK_SET);
      }
   else if(field_type == MEDIA_FIELD)
      {
      lseek(mds->med_fdr,0,SEEK_SET);
      lseek(mds->med_fdw,0,SEEK_SET);
      }
   }
}

set_model_seek(mds,field_type,ip)
struct modelstorage *mds;
int field_type, ip;
{
if(mds->intmem)
   {
   if(field_type == WAVE_FIELD)
      mds->pv_bufptr = mds->pv_buf + ip*mds->ln_pvslice;
   else if(field_type == MEDIA_FIELD)
      mds->med_bufptr = mds->med_buf + ip*mds->ln_medslice;
   }
else
   {
   if(field_type == WAVE_FIELD)
      {
      lseek(mds->pv_fdr,ip*mds->sz_pvslice,SEEK_SET);
      lseek(mds->pv_fdw,ip*mds->sz_pvslice,SEEK_SET);
      }
   else if(field_type == MEDIA_FIELD)
      {
      lseek(mds->med_fdr,ip*mds->sz_medslice,SEEK_SET);
      lseek(mds->med_fdw,ip*mds->sz_medslice,SEEK_SET);
      }
   }
}

float *reed_model(struct modelstorage *mds,float *fptr,int field_type)
{
float *tptr;

if(mds->intmem)
   {
   if(field_type == WAVE_FIELD)
      {
      tptr = mds->pv_bufptr;
      mds->pv_bufptr = mds->pv_bufptr + mds->ln_pvslice;
      }
   else if(field_type == MEDIA_FIELD)
      {
      tptr = mds->med_bufptr;
      mds->med_bufptr = mds->med_bufptr + mds->ln_medslice;
      }
   }
else
   {
   if(field_type == WAVE_FIELD)
      reed(mds->pv_fdr,fptr,mds->sz_pvslice);
   else if(field_type == MEDIA_FIELD)
      reed(mds->med_fdr,fptr,mds->sz_medslice);

   tptr = fptr;
   }
return(tptr);
}
 
rite_model(mds,fptr,field_type)
struct modelstorage *mds;
float *fptr;
int field_type;
{
if(!mds->intmem)
   {
   if(field_type == WAVE_FIELD)
      rite(mds->pv_fdw,fptr,mds->sz_pvslice);
   else if(field_type == MEDIA_FIELD)
      rite(mds->med_fdw,fptr,mds->sz_medslice);
   }
}

/*
   Model space will either be stored internally in core memory or
   in a temporary file on external disk.
*/

void set_model_storage(struct modelstorage *modstor,struct nodeinfo *ni,int nx,int ny,int nz,int *span,int *spanlen,char *name,char *tmpdir)
{
size_t memlen;
char string[1024];

modstor->ln_pvslice = N_WAVE_VARS*nx*nz;
modstor->sz_pvslice = modstor->ln_pvslice*sizeof(float);
modstor->ln_medslice = N_MED_VARS*nx*nz;
modstor->sz_medslice = modstor->ln_medslice*sizeof(float);

if(modstor->intmem)         /* setup internal model storage */
   {
   fprintf(stderr,"**** Internal memory only\n\n");

   memlen = ny*(modstor->sz_pvslice + modstor->sz_medslice);

   if(memlen > modstor->maxmem)
      {
      fprintf(stderr,"**** memlen= %ld, maxmem= %ld\n",memlen,modstor->maxmem);
      fprintf(stderr,"**** Not enough memory available, exiting...\n");
      mpi_exit(-1);
      }

   modstor->pv_buf = (float *) check_malloc (ny*modstor->sz_pvslice);
   modstor->med_buf = (float *) check_malloc (ny*modstor->sz_medslice);
   }
else                       /* open temporary model file */
   {
   fprintf(stderr,"\n");
   fprintf(stderr,"**** Temporary files in dir= %s\n\n",tmpdir);
   fflush(stderr);

   makedir(tmpdir);

   sprintf(string,"%s/%s_pv-%.5d.e3d",tmpdir,name,ni->segmentId);
   modstor->pv_fdr = croptrfile(string);
   modstor->pv_fdw = opfile(string);

   sprintf(string,"%s/%s_med-%.5d.e3d",tmpdir,name,ni->segmentId);
   modstor->med_fdr = croptrfile(string);
   modstor->med_fdw = opfile(string);
   }

/*
   MPI-RWG:
   For distributed computation, spanlen must have 2 additional planes
   to allow for communication

   OLD:  spanlen = 4*span + 2
   NEW:  spanlen = 4*span + 4
*/

if((*span) < 1)
   *span = 1;
*spanlen = 4*(*span) + 4;

while((*spanlen) > ny)
   {
   (*span)--;
   *spanlen = 4*(*span) + 4;
   }

memlen = (*spanlen)*nx*nz*sizeof(float);
while((N_MED_VARS+N_WAVE_VARS)*memlen > modstor->maxmem)
   {
   (*span)--;
   *spanlen = 4*(*span) + 4;
   memlen = (*spanlen)*nx*nz*sizeof(float);

   if((*span) == 0)
      {
      fprintf(stderr,"**** maxmem= %ld\n",modstor->maxmem);
      fprintf(stderr,"**** Not enough memory available, exiting...\n");
      mpi_exit(-1);
      }
   }
}

size_t set_model_storageP3(struct modelstorage *modstor,int nx,int ny,int nz)
{
size_t nbuf, memlen, blen;

modstor->ln_pvslice = N_WAVE_VARS*nx*nz;
modstor->sz_pvslice = modstor->ln_pvslice*sizeof(float);
modstor->ln_medslice = N_MED_VARS*nx*nz;
modstor->sz_medslice = modstor->ln_medslice*sizeof(float);

nbuf = 4*(N_WAVE_VARS-3)*nx*nz;
if(4*(N_WAVE_VARS-3)*nx*ny > nbuf)
   nbuf = 4*(N_WAVE_VARS-3)*nx*ny;
if(4*(N_WAVE_VARS-3)*nz*ny > nbuf)
   nbuf = 4*(N_WAVE_VARS-3)*nz*ny;

fprintf(stderr,"**** Internal memory only\n\n");

memlen = sizeof(float)*((modstor->ln_pvslice + modstor->ln_medslice)*ny + nbuf + 6*(nx+nz));

if(memlen > modstor->maxmem)
   {
   fprintf(stderr,"**** memlen= %ld, maxmem= %ld\n",memlen,modstor->maxmem);
   fprintf(stderr,"**** Not enough memory available, exiting...\n");
   mpi_exit(-1);
   }

blen = ny*modstor->ln_pvslice + nbuf + 6*(nx+nz);
modstor->pv_buf = (float *) check_malloc (blen*sizeof(float));

blen = ny*modstor->ln_pvslice;
modstor->med_buf = (float *) check_malloc (blen*sizeof(float));

return(memlen);
}
