#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

#include <libmemcached/memcached.h>

//void timeshift_sgt(float *seis,int ntout,float *gf,struct sgtheader *gfh,int ntsum,float *t0,float *bt0,int nsgt);
void timeshift_sgt(float* seis, int ntout, float* gf, int ntsum, float* t0, float* bt0, struct sgtparams* sgtpar);

void sgt_subset(struct sgtfileparams *sgtfpar,struct sgtfileparams *sgtextract,struct sgtmaster *sgtmast,struct sgtindex *sgtindx,int nm,long long *mindx,char *dir)
{
struct sgtmaster exmast;
struct sgtindex *exindx;
struct sgtheader sgthead;
float *sgt;
int im, ip, nreed;
off_t blen, off;
char ofile[512];

exmast.geoproj = sgtmast->geoproj;
exmast.modellon = sgtmast->modellon;
exmast.modellat = sgtmast->modellat;
exmast.modelrot = sgtmast->modelrot;
exmast.xshift = sgtmast->xshift;
exmast.yshift = sgtmast->yshift;
exmast.globnp = nm;
exmast.localnp = nm;
exmast.nt = sgtmast->nt;

exindx = (struct sgtindex *) check_malloc ((exmast.globnp)*sizeof(struct sgtindex));

for(im=0;im<nm;im++)
   {
   ip = 0;
   while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
      ip++;

   exindx[im].indx = sgtindx[ip].indx;
   exindx[im].xsgt = sgtindx[ip].xsgt;
   exindx[im].ysgt = sgtindx[ip].ysgt;
   exindx[im].zsgt = sgtindx[ip].zsgt;
   exindx[im].h = sgtindx[ip].h;

   if(mindx[im] != exindx[im].indx)
      {
      fprintf(stderr,"houston, we have a problem...\n");
      exit(-1);
      }
   }

sgt = (float *) check_malloc (6*(sgtmast->nt)*sizeof(float));
blen = (off_t)(sizeof(struct sgtheader)) + (off_t)(6*(sgtmast->nt)*sizeof(float));

if(strcmp(dir,".") != 0) {
   makedir(dir);
}

//try creating big buffer to accumulate extracted sgt data to reduce # of writes
char* extracted_sgt = check_malloc(nm*blen);
memset(extracted_sgt, 0, blen*nm);
int extracted_entry_size = blen;

if(sgtfpar->xfile[0] != '\0')
   {
   if(sgtextract->xfile[0] == '/')
      sprintf(ofile,"%s",sgtextract->xfile);
   else
      sprintf(ofile,"%s/%s",dir,sgtextract->xfile);

   fcroptrfile(ofile, &(sgtextract->xfp));
   //sgtextract->xfp = croptrfile(ofile);

   frite(sgtextract->xfp,&exmast,sizeof(struct sgtmaster), 1);
   //rite(sgtextract->xfdr,&exmast,sizeof(struct sgtmaster));
   frite(sgtextract->xfp,exindx,sizeof(struct sgtindex), exmast.globnp);
   //rite(sgtextract->xfdr,exindx,(exmast.globnp)*sizeof(struct sgtindex));

   if (sgtfpar->xfile_header[0] != '\0') {
	//Change because the headers aren't in the SGT file
	blen = (off_t)(6*(sgtmast->nt)*sizeof(float));
	//set head off past sgtmaster and sgtindices
        sgtfpar->x_head_off = sizeof(struct sgtmaster) + sgtmast->globnp*sizeof(struct sgtindex);
        //also, open SGT file
        fopfile_ro(sgtfpar->xfile, &(sgtfpar->xfp));
        //sgtfpar->xfdr = opfile_ro(sgtfpar->xfile);
	sgtfpar->xcur_off = 0;
   }

   for(im=0;im<nm;im++)
      {
      ip = 0;
      while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
         ip++;

      if (sgtfpar->xfile_header[0] != '\0') {
        off = (off_t)(ip)*blen - sgtfpar->xcur_off;
        fseek(sgtfpar->x_head_fp, sgtfpar->x_head_off + ip*sizeof(struct sgtheader), SEEK_SET);
        //lseek(sgtfpar->x_head_fdr, sgtfpar->x_head_off + ip*sizeof(struct sgtheader), SEEK_SET);
      } else {
        off = sgtfpar->head_off + (off_t)(ip)*blen - sgtfpar->xcur_off;
      }
      fseek(sgtfpar->xfp,off,SEEK_CUR);
      //lseek(sgtfpar->xfdr,off,SEEK_CUR);
      sgtfpar->xcur_off = sgtfpar->head_off + (off_t)(ip)*blen;


      if (sgtfpar->xfile_header[0] != '\0') {
	//read sgtheaders from the separate header file
        nreed = freed(sgtfpar->x_head_fp,&sgthead,sizeof(struct sgtheader), 1);
	//nreed = reed(sgtfpar->x_head_fdr,&sgthead,sizeof(struct sgtheader));
      } else {
	nreed = freed(sgtfpar->xfp,&sgthead,sizeof(struct sgtheader), 1);
        //nreed = reed(sgtfpar->xfdr,&sgthead,sizeof(struct sgtheader));
        sgtfpar->xcur_off = sgtfpar->xcur_off + (off_t)(nreed);
      }
      nreed = freed(sgtfpar->xfp,sgt,sizeof(float),6*sgtmast->nt);
      //nreed = reed(sgtfpar->xfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      sgtfpar->xcur_off = sgtfpar->xcur_off + (off_t)(nreed);

      //frite(sgtextract->xfp,&sgthead,sizeof(struct sgtheader), 1);
      //rite(sgtextract->xfdr,&sgthead,sizeof(struct sgtheader));
      //frite(sgtextract->xfp,sgt,sizeof(float),6*sgtmast->nt);
      //rite(sgtextract->xfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      memcpy(extracted_sgt+im*extracted_entry_size, &sgthead, sizeof(struct sgtheader));
      memcpy(extracted_sgt+im*extracted_entry_size+sizeof(struct sgtheader), sgt, 6*sgtmast->nt*sizeof(float));
      }

   frite(sgtextract->xfp,extracted_sgt,1,nm*extracted_entry_size);
   fflush(sgtextract->xfp);

   if (sgtfpar->xfile_header[0] != '\0') {
        fclose(sgtfpar->x_head_fp);
	//close(sgtfpar->x_head_fdr);
   }
   fclose(sgtfpar->xfp);
   //close(sgtfpar->xfdr);
   fclose(sgtextract->xfp);
   //close(sgtextract->xfdr);
   }

//set blen back
blen = (off_t)(sizeof(struct sgtheader)) + (off_t)(6*(sgtmast->nt)*sizeof(float));
sgtfpar->head_off = (off_t)(sizeof(struct sgtmaster)) + (off_t)((sgtmast->globnp)*sizeof(struct sgtindex));
memset(extracted_sgt, 0, nm*blen);

if(sgtfpar->yfile[0] != '\0')
   {
   if(sgtextract->yfile[0] == '/')
      sprintf(ofile,"%s",sgtextract->yfile);
   else
      sprintf(ofile,"%s/%s",dir,sgtextract->yfile);

   fcroptrfile(ofile, &(sgtextract->yfp));
   //sgtextract->yfdr = croptrfile(ofile);

   frite(sgtextract->yfp,&exmast,sizeof(struct sgtmaster), 1);
   //rite(sgtextract->yfdr,&exmast,sizeof(struct sgtmaster));
   frite(sgtextract->yfp,exindx,sizeof(struct sgtindex), exmast.globnp);
   //rite(sgtextract->yfdr,exindx,(exmast.globnp)*sizeof(struct sgtindex));

   if (sgtfpar->yfile_header[0] != '\0') {
	blen = (off_t)(6*(sgtmast->nt)*sizeof(float));
	sgtfpar->head_off = 0;
	//Set offset in file header past sgtmaster and sgtindexes
        sgtfpar->y_head_off = sizeof(struct sgtmaster) + sgtmast->globnp*sizeof(struct sgtindex);
        //also, open SGT file
        fopfile_ro(sgtfpar->yfile,&(sgtfpar->yfp));
        //sgtfpar->yfdr = opfile_ro(sgtfpar->yfile);
        sgtfpar->ycur_off = 0;
   }

   for(im=0;im<nm;im++)
      {
      ip = 0;
      while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
         ip++;

      if(sgtfpar->yfile_header[0] != '\0') {
        off = (off_t)(ip)*blen - sgtfpar->ycur_off;
        fseek(sgtfpar->y_head_fp, sgtfpar->y_head_off + ip*sizeof(struct sgtheader), SEEK_SET);
        //lseek(sgtfpar->y_head_fdr, sgtfpar->y_head_off + ip*sizeof(struct sgtheader), SEEK_SET);
      } else {
        off = sgtfpar->head_off + (off_t)(ip)*blen - sgtfpar->ycur_off;
      }
      fseek(sgtfpar->yfp,off,SEEK_CUR);
      //lseek(sgtfpar->yfdr,off,SEEK_CUR);
      sgtfpar->ycur_off = sgtfpar->head_off + (off_t)(ip)*blen;

      if (sgtfpar->yfile_header[0] != '\0') {
	nreed = freed(sgtfpar->y_head_fp,&sgthead,sizeof(struct sgtheader), 1);
        //nreed = reed(sgtfpar->y_head_fdr,&sgthead,sizeof(struct sgtheader));
      } else {
	nreed = freed(sgtfpar->yfp,&sgthead,sizeof(struct sgtheader), 1);
        //nreed = reed(sgtfpar->yfdr,&sgthead,sizeof(struct sgtheader));
        sgtfpar->ycur_off = sgtfpar->ycur_off + (off_t)(nreed);
      }

      nreed = freed(sgtfpar->yfp,sgt,sizeof(float),6*sgtmast->nt);
      //nreed = reed(sgtfpar->yfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      sgtfpar->ycur_off = sgtfpar->ycur_off + (off_t)(nreed);

      //frite(sgtextract->yfp,&sgthead,sizeof(struct sgtheader),1);
      //rite(sgtextract->yfdr,&sgthead,sizeof(struct sgtheader));
      //frite(sgtextract->yfp,sgt,sizeof(float),6*sgtmast->nt);
      //rite(sgtextract->yfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      memcpy(extracted_sgt+im*extracted_entry_size, &sgthead, sizeof(struct sgtheader));
      memcpy(extracted_sgt+im*extracted_entry_size+sizeof(struct sgtheader), sgt, 6*sgtmast->nt*sizeof(float));
      }
   
   frite(sgtextract->yfp, extracted_sgt, 1, nm*extracted_entry_size);
   fflush(sgtextract->yfp);

   if (sgtfpar->yfile_header[0] != '\0') {
        fclose(sgtfpar->y_head_fp);
	//close(sgtfpar->y_head_fdr);
   }

   fclose(sgtfpar->yfp);
   //close(sgtfpar->yfdr);
   fclose(sgtextract->yfp);
   //close(sgtextract->yfdr);
   }

//set blen back
blen = (off_t)(sizeof(struct sgtheader)) + (off_t)(6*(sgtmast->nt)*sizeof(float));
sgtfpar->head_off = (off_t)(sizeof(struct sgtmaster)) + (off_t)((sgtmast->globnp)*sizeof(struct sgtindex));
memset(extracted_sgt, 0, nm*blen);

if(sgtfpar->zfile[0] != '\0')
   {
   if(sgtextract->zfile[0] == '/')
      sprintf(ofile,"%s",sgtextract->zfile);
   else
      sprintf(ofile,"%s/%s",dir,sgtextract->zfile);

   fcroptrfile(ofile, &(sgtextract->zfp));
   //sgtextract->zfdr = croptrfile(ofile);

   frite(sgtextract->zfp,&exmast,sizeof(struct sgtmaster), 1);
   //rite(sgtextract->zfdr,&exmast,sizeof(struct sgtmaster));
   frite(sgtextract->zfp,exindx,sizeof(struct sgtindex), exmast.globnp);
   //rite(sgtextract->zfdr,exindx,(exmast.globnp)*sizeof(struct sgtindex));

   if (sgtfpar->zfile_header[0] != '\0') {
        blen = (off_t)(6*(sgtmast->nt)*sizeof(float));
        sgtfpar->head_off = 0;
        //Set offset in file header past sgtmaster and sgtindexes
        sgtfpar->z_head_off = sizeof(struct sgtmaster) + sgtmast->globnp*sizeof(struct sgtindex);
        //also, open SGT file
	fopfile_ro(sgtfpar->zfile, &(sgtfpar->zfp));
        //sgtfpar->zfdr = opfile_ro(sgtfpar->zfile);
        sgtfpar->zcur_off = 0;
   }

   for(im=0;im<nm;im++)
      {
      ip = 0;
      while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
         ip++;

      if (sgtfpar->zfile_header[0] != '\0') {
	off = (off_t)(ip)*blen - sgtfpar->zcur_off;
	fseek(sgtfpar->z_head_fp,sgtfpar->z_head_off + ip*sizeof(struct sgtheader), SEEK_SET);
	//lseek(sgtfpar->z_head_fdr, sgtfpar->z_head_off + ip*sizeof(struct sgtheader), SEEK_SET);
      } else {
        off = sgtfpar->head_off + (off_t)(ip)*blen - sgtfpar->zcur_off;
      }
      fseek(sgtfpar->zfp,off,SEEK_CUR);
      //lseek(sgtfpar->zfdr,off,SEEK_CUR);
      sgtfpar->zcur_off = sgtfpar->head_off + (off_t)(ip)*blen;

      if (sgtfpar->zfile_header[0] != '\0') {
	nreed = freed(sgtfpar->z_head_fp,&sgthead,sizeof(struct sgtheader),1);
	//nreed = reed(sgtfpar->z_head_fdr,&sgthead,sizeof(struct sgtheader));
      } else {
	nreed = freed(sgtfpar->zfp,&sgthead,sizeof(struct sgtheader), 1);
      	//nreed = reed(sgtfpar->zfdr,&sgthead,sizeof(struct sgtheader));
        sgtfpar->zcur_off = sgtfpar->zcur_off + (off_t)(nreed);
      }

      nreed = freed(sgtfpar->zfp,sgt,sizeof(float),6*(sgtmast->nt));
      //nreed = reed(sgtfpar->zfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      sgtfpar->zcur_off = sgtfpar->zcur_off + (off_t)(nreed);

      //frite(sgtextract->zfp,&sgthead,sizeof(struct sgtheader), 1);
      //rite(sgtextract->zfdr,&sgthead,sizeof(struct sgtheader));
      //frite(sgtextract->zfp,sgt,sizeof(float),6*sgtmast->nt);
      //rite(sgtextract->zfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      memcpy(extracted_sgt+im*extracted_entry_size, &sgthead, sizeof(struct sgtheader));
      memcpy(extracted_sgt+im*extracted_entry_size+sizeof(struct sgtheader), sgt, 6*sgtmast->nt*sizeof(float));
      }

   frite(sgtextract->zfp,extracted_sgt,1,nm*extracted_entry_size);
   fflush(sgtextract->zfp);

   if (sgtfpar->zfile_header[0] != '\0') {
	fclose(sgtfpar->z_head_fp);
	//close(sgtfpar->z_head_fdr);
   }

   fclose(sgtfpar->zfp);
   //close(sgtfpar->zfdr);
   fclose(sgtextract->zfp);
   //close(sgtextract->zfdr);
   }

   free(extracted_sgt);
}

void get_sgtpars(struct sgtfileparams *sgtfpar,struct sgtmaster *sgtmast,struct sgtindex **sgtindx)
{
struct sgtmaster* tmast;
struct sgtindex *tindx;
struct sgtmaster *sgtmast_tmp;
struct sgtmaster *tmast_tmp;
struct sgtindex *tindx_tmp;

int ip;

int eflag = 0;
int rflag = 0;

printf("In getpars.\n");

memcached_return ret;
memcached_st* mst;
int i;

mst = memcached_create(NULL);
memcached_server_st *server;
char server_name[] = "localhost";
server = memcached_server_list_append(NULL, server_name, 11211, &ret);

ret = memcached_server_push(mst, server);
memcached_server_free(server);

size_t sgtmaster_size = (size_t)sizeof(struct sgtmaster);
size_t incoming_size;
size_t sgtindex_size;
//This is what memcached thinks 1 MB is
int one_mb = 1000*1000;
uint32_t flags = 0;

if(sgtfpar->xfile[0] != '\0')
   {
   char* filename_with_header;
   FILE* fp_with_header;

   if(sgtfpar->xfile_header[0] != '\0') {
	printf("Using separate header file %s for x.\n", sgtfpar->xfile_header);
	filename_with_header = sgtfpar->xfile_header;
	//It would be nice if we could avoid this, but since we're not caching the individual point headers, we'll have to open it anyway
	fopfile_ro(sgtfpar->xfile_header, &(sgtfpar->x_head_fp));
        //sgtfpar->x_head_fdr = opfile_ro(sgtfpar->xfile_header);
	fp_with_header = sgtfpar->x_head_fp;
	//fp_with_header = sgtfpar->x_head_fdr;
   } else {
	printf("Header is part of SGTs.\n");
	filename_with_header = sgtfpar->xfile;
	fopfile_ro(sgtfpar->xfile, &(sgtfpar->xfp));
	//sgtfpar->xfdr = opfile_ro(sgtfpar->xfile);
	fp_with_header = sgtfpar->xfp;
	//fp_with_header = sgtfpar->xfdr;
   }
   
   char* mast_key = check_malloc(sizeof(char) * (strlen(filename_with_header) + strlen("_sgtmaster") + 1));
   sprintf(mast_key, "%s%s", filename_with_header, "_sgtmaster");
   
   sgtmast_tmp = (struct sgtmaster *) memcached_get(mst, mast_key, strlen(mast_key), &incoming_size, &flags, &ret);
   
   if (sgtmast_tmp==NULL) {
   	printf("Key %s hasn't been cached, adding.\n", mast_key);
	freed(fp_with_header,sgtmast,sizeof(struct sgtmaster), 1);
	//reed(fp_with_header,sgtmast,sizeof(struct sgtmaster));
        ret = memcached_set(mst, mast_key, strlen(mast_key), (char*) sgtmast, sgtmaster_size, 0, flags);
	if (ret!=MEMCACHED_SUCCESS) {
		//Continue running even if caching fails
        	fprintf(stderr, "Caching of %s failed.\n", mast_key);
		fprintf(stderr,  "%s\n", memcached_strerror(mst, ret));
        } else {
                printf("Key %s added.\n", mast_key);
        }
   } else {
        printf("Retrieved key %s from cache.\n", mast_key);
	//Have to put it in the sgtmast ptr since we will be passing it around in jbsim3d.c
        memcpy(sgtmast, sgtmast_tmp, sizeof(struct sgtmaster));
   }
   free(mast_key);

   sgtindex_size = (size_t) (sgtmast->globnp * sizeof(struct sgtindex));
   *sgtindx = check_malloc(sgtindex_size);

   int num_pieces = sgtindex_size/one_mb;
   if (num_pieces * one_mb < sgtindex_size) {
   	num_pieces++;
   }

   char* buf;
   char* index_key = check_malloc(sizeof(char) * (strlen(filename_with_header) + strlen("_sgtheader") + 4));

   for (i=0; i<num_pieces; i++) {
   	int receive = one_mb;
        if (sgtindex_size-i*one_mb<one_mb) {
   	     receive = sgtindex_size - i*one_mb;
        }
        sprintf(index_key, "%s_sgtheader.%d", filename_with_header, i);
        buf = memcached_get(mst, index_key, strlen(index_key), &incoming_size, &flags, &ret);
        if (buf==NULL) {
	        //need to add
                printf("%s is missing, adding to cache.\n", index_key);
                buf = check_malloc(receive);
		fseek(fp_with_header, sgtmaster_size + i*one_mb, SEEK_SET);
                //lseek(fp_with_header, sgtmaster_size + i*one_mb, SEEK_SET);
		freed(fp_with_header, buf, 1, receive);
                //reed(fp_with_header, buf, receive);
                ret = memcached_set(mst, index_key, strlen(index_key), buf, receive, 0, flags);
                if (ret!=MEMCACHED_SUCCESS) {
	                fprintf(stderr, "Caching %s failed.\n", index_key);
                        fprintf(stderr,  "%s\n", memcached_strerror(mst, ret));
                } else {
			fprintf(stderr, "Caching %s successful.\n", index_key);
		}
                //Cast to char* so our pointer arithmetic works correctly
                memcpy((char*)(*sgtindx)+i*one_mb, buf, receive);
                free(buf);
	} else {
                printf("Retrieved key %s from cache.\n", index_key);
                memcpy((char*)(*sgtindx)+i*one_mb, buf, receive);
        }
   }
   fseek(fp_with_header, sgtmaster_size + sgtindex_size, SEEK_SET);
   //lseek(fp_with_header, sgtmaster_size + sgtindex_size, SEEK_SET);
   free(index_key);

   if(sgtfpar->xfile_header[0] != '\0') {
	sgtfpar->head_off = 0;
   } else {
	sgtfpar->head_off = (off_t)(sgtmaster_size + sgtindex_size);
   }
   sgtfpar->xcur_off = sgtfpar->head_off;
 
   rflag = 1;
}

if(sgtfpar->yfile[0] != '\0')
   {
   char* filename_with_header;
   FILE* fp_with_header;

   if(sgtfpar->yfile_header[0] != '\0') {
        printf("Using separate header file for y.\n");
        filename_with_header = sgtfpar->yfile_header;
        //It would be nice if we could avoid this, but since we're not caching the individual point headers, we'll have to open it anyway
	fopfile_ro(sgtfpar->yfile_header, &(sgtfpar->y_head_fp));
        //sgtfpar->y_head_fdr = opfile_ro(sgtfpar->yfile_header);
	fp_with_header = sgtfpar->y_head_fp;
        //fp_with_header = sgtfpar->y_head_fdr;
   } else {
        printf("Header is part of SGTs.\n");
        filename_with_header = sgtfpar->yfile;
	fopfile_ro(sgtfpar->yfile, &(sgtfpar->yfp));
        //sgtfpar->yfdr = opfile_ro(sgtfpar->yfile);
	fp_with_header = sgtfpar->yfp;
        //fp_with_header = sgtfpar->yfdr;
   }

   char* mast_key = check_malloc(sizeof(char) * (strlen(filename_with_header) + strlen("_sgtmaster") + 1));
   sprintf(mast_key, "%s%s", filename_with_header, "_sgtmaster");

   tmast = check_malloc(sgtmaster_size);
   tmast_tmp = (struct sgtmaster *) memcached_get(mst, mast_key, strlen(mast_key), &incoming_size, &flags, &ret);

   if (tmast_tmp==NULL) {
   	printf("Key %s hasn't been cached, adding.\n", mast_key);
	freed(fp_with_header,tmast,sgtmaster_size, 1);
        //reed(fp_with_header,tmast,sgtmaster_size);
        ret = memcached_set(mst, mast_key, strlen(mast_key), (char*) tmast, sgtmaster_size, 0, flags);
        if (ret!=MEMCACHED_SUCCESS) {
	        fprintf(stderr, "Caching of %s failed.\n", mast_key);
		fprintf(stderr,  "%s\n", memcached_strerror(mst, ret));
        } else {
                printf("Key %s added.\n", mast_key);
        }
   } else {
        printf("Retrieved key %s from cache.\n", mast_key);
	//Try getting rid of this copy
        memcpy(tmast, tmast_tmp, sgtmaster_size);
   }
   free(mast_key);

   if(rflag == 0) {
	//copy to sgtmast
	memcpy(sgtmast, tmast, sgtmaster_size);
   } else {
      eflag = 0;
      if(tmast->geoproj != sgtmast->geoproj)
         eflag = -1;
      if(tmast->modellon > sgtmast->modellon + 0.001 || tmast->modellon < sgtmast->modellon - 0.001)
         eflag = -1;
      if(tmast->modellat > sgtmast->modellat + 0.001 || tmast->modellat < sgtmast->modellat - 0.001)
         eflag = -1;
      if(tmast->modelrot > sgtmast->modelrot + 0.001 || tmast->modelrot < sgtmast->modelrot - 0.001)
         eflag = -1;
      if(tmast->xshift > sgtmast->xshift + 0.001 || tmast->xshift < sgtmast->xshift - 0.001)
         eflag = -1;
      if(tmast->yshift > sgtmast->yshift + 0.001 || tmast->yshift < sgtmast->yshift - 0.001)
         eflag = -1;
      if(tmast->globnp != sgtmast->globnp)
         eflag = -1;
      if(tmast->localnp != sgtmast->localnp)
         eflag = -1;
      if(tmast->nt != sgtmast->nt)
         eflag = -1;

      if(eflag != 0)
         {
         fprintf(stderr,"sgtmaster inconsistency in yfile= %s, exiting ...\n",sgtfpar->yfile);
         exit(-1);
         }
   }

   sgtindex_size = (size_t) (sgtmast->globnp * sizeof(struct sgtindex));
   tindx = check_malloc(sgtindex_size);

   int num_pieces = sgtindex_size/one_mb;
   if (num_pieces * one_mb < sgtindex_size) {
        num_pieces++;
   }

   char* buf;
   char* index_key = check_malloc(sizeof(char) * (strlen(filename_with_header) + strlen("_sgtheader") + 4));

   for (i=0; i<num_pieces; i++) {
        int receive = one_mb;
        if (sgtindex_size-i*one_mb<one_mb) {
             receive = sgtindex_size - i*one_mb;
        }
        sprintf(index_key, "%s_sgtheader.%d", filename_with_header, i);
        buf = memcached_get(mst, index_key, strlen(index_key), &incoming_size, &flags, &ret);
        if (buf==NULL) {
                //need to add
                printf("%s is missing, adding to cache.\n", index_key);
                buf = check_malloc(receive);
		fseek(fp_with_header, sgtmaster_size + i*one_mb, SEEK_SET);
                //lseek(fp_with_header, sgtmaster_size + i*one_mb, SEEK_SET);
		freed(fp_with_header, buf, 1, receive);
                //reed(fp_with_header, buf, receive);
                ret = memcached_set(mst, index_key, strlen(index_key), buf, receive, 0, flags);
                if (ret!=MEMCACHED_SUCCESS) {
                        fprintf(stderr, "Caching %s failed.\n", index_key);
                        fprintf(stderr,  "%s\n", memcached_strerror(mst, ret));
                }
                //Cast to char* so our pointer arithmetic works correctly
                memcpy((char*)tindx+i*one_mb, buf, receive);
                free(buf);
        } else {
                printf("Retrieved key %s from cache.\n", index_key);
                memcpy((char*)tindx+i*one_mb, buf, receive);
        }
   }
   fseek(fp_with_header, sgtmaster_size + sgtindex_size, SEEK_SET);
   //lseek(fp_with_header, sgtmaster_size + sgtindex_size, SEEK_SET);
   free(index_key);

   if (rflag==0) {
	memcpy((char*)(*sgtindx), tindx, sgtindex_size);
	rflag = 1;
   } else {
      //Do index check
      for(ip=0;ip<sgtmast->globnp;ip++)
         {
         if(tindx[ip].indx != (*sgtindx)[ip].indx) {
            eflag = -1;
         }
         }

      if(eflag != 0)
         {
         fprintf(stderr,"sgtindex inconsistency in yfile= %s, exiting ...\n",sgtfpar->yfile);
         exit(-1);
         }
   }
   free(tindx);
   free(tmast);

   if(sgtfpar->yfile_header[0] != '\0') {
        sgtfpar->head_off = 0;
   } else {
        sgtfpar->head_off = (off_t)(sgtmaster_size + sgtindex_size);
   }
   sgtfpar->ycur_off = sgtfpar->head_off;
}


if(sgtfpar->zfile[0] != '\0')
   {
   char* filename_with_header;
   FILE* fp_with_header;
   
   if(sgtfpar->zfile_header[0] != '\0') {
        printf("Using separate header file for z.\n");
        filename_with_header = sgtfpar->zfile_header;
        //It would be nice if we could avoid this, but since we're not caching the individual point headers, we'll have to open it anyway
	fopfile_ro(sgtfpar->zfile_header, &(sgtfpar->z_head_fp));
        //sgtfpar->z_head_fdr = opfile_ro(sgtfpar->zfile_header);
	fp_with_header = sgtfpar->z_head_fp;
        //fp_with_header = sgtfpar->z_head_fdr;
   } else {
        printf("Header is part of SGTs.\n");
        filename_with_header = sgtfpar->zfile;
	fopfile_ro(sgtfpar->zfile, &(sgtfpar->zfp));
        //sgtfpar->zfdr = opfile_ro(sgtfpar->zfile);
	fp_with_header = sgtfpar->zfp;
        //fp_with_header = sgtfpar->zfdr;
   }

   char* mast_key = check_malloc(sizeof(char) * (strlen(filename_with_header) + strlen("_sgtmaster") + 1));
   sprintf(mast_key, "%s%s", filename_with_header, "_sgtmaster");

   tmast = check_malloc(sgtmaster_size);
   tmast_tmp = (struct sgtmaster *) memcached_get(mst, mast_key, strlen(mast_key), &incoming_size, &flags, &ret);

   if (tmast_tmp==NULL) {
        printf("Key %s hasn't been cached, adding.\n", mast_key);
	freed(fp_with_header,tmast,sgtmaster_size, 1);
        //reed(fp_with_header,tmast,sgtmaster_size);
        ret = memcached_set(mst, mast_key, strlen(mast_key), (char*) tmast, sgtmaster_size, 0, flags);
        if (ret!=MEMCACHED_SUCCESS) {
                fprintf(stderr, "Caching of %s failed.\n", mast_key);
		fprintf(stderr,  "%s\n", memcached_strerror(mst, ret));
        } else {
                printf("Key %s added.\n", mast_key);
        }
   } else {
        printf("Retrieved key %s from cache.\n", mast_key);
        //Try getting rid of this copy
        memcpy(tmast, tmast_tmp, sgtmaster_size);
   }
   free(mast_key);

   if(rflag == 0) {
        //copy to sgtmast
        memcpy(sgtmast, tmast, sgtmaster_size);
   } else {
      eflag = 0;
      if(tmast->geoproj != sgtmast->geoproj)
         eflag = -1;
      if(tmast->modellon > sgtmast->modellon + 0.001 || tmast->modellon < sgtmast->modellon - 0.001)
         eflag = -1;
      if(tmast->modellat > sgtmast->modellat + 0.001 || tmast->modellat < sgtmast->modellat - 0.001)
         eflag = -1;
      if(tmast->modelrot > sgtmast->modelrot + 0.001 || tmast->modelrot < sgtmast->modelrot - 0.001)
         eflag = -1;
      if(tmast->xshift > sgtmast->xshift + 0.001 || tmast->xshift < sgtmast->xshift - 0.001)
         eflag = -1;
      if(tmast->yshift > sgtmast->yshift + 0.001 || tmast->yshift < sgtmast->yshift - 0.001)
         eflag = -1;
      if(tmast->globnp != sgtmast->globnp)
         eflag = -1;
      if(tmast->localnp != sgtmast->localnp)
         eflag = -1;
      if(tmast->nt != sgtmast->nt)
         eflag = -1;

      if(eflag != 0)
         {
         fprintf(stderr,"sgtmaster inconsistency in zfile= %s, exiting ...\n",sgtfpar->zfile);
         exit(-1);
         }
   }

   sgtindex_size = (size_t) (sgtmast->globnp * sizeof(struct sgtindex));
   tindx = check_malloc(sgtindex_size);
   
   int num_pieces = sgtindex_size/one_mb;
   if (num_pieces * one_mb < sgtindex_size) {
        num_pieces++;
   }
   
   char* buf;
   char* index_key = check_malloc(sizeof(char) * (strlen(filename_with_header) + strlen("_sgtheader") + 4));

   for (i=0; i<num_pieces; i++) {
        int receive = one_mb;
        if (sgtindex_size-i*one_mb<one_mb) {
             receive = sgtindex_size - i*one_mb;
        }
        sprintf(index_key, "%s_sgtheader.%d", filename_with_header, i);
        buf = memcached_get(mst, index_key, strlen(index_key), &incoming_size, &flags, &ret);
        if (buf==NULL) {
                //need to add
                printf("%s is missing, adding to cache.\n", index_key);
                buf = check_malloc(receive);
		fseek(fp_with_header, sgtmaster_size + i*one_mb, SEEK_SET);
                //lseek(fp_with_header, sgtmaster_size + i*one_mb, SEEK_SET);
		freed(fp_with_header, buf, receive, 1);
                //reed(fp_with_header, buf, receive);
                ret = memcached_set(mst, index_key, strlen(index_key), buf, receive, 0, flags);
                if (ret!=MEMCACHED_SUCCESS) {
                        fprintf(stderr, "Caching %s failed, aborting.\n", index_key);
                        fprintf(stderr,  "%s\n", memcached_strerror(mst, ret));
                }
                //Cast to char* so our pointer arithmetic works correctly
                memcpy((char*)tindx+i*one_mb, buf, receive);
                free(buf);
        } else {
                printf("Retrieved key %s from cache.\n", index_key);
                memcpy((char*)tindx+i*one_mb, buf, receive);
        }
   }
   fseek(fp_with_header, sgtmaster_size + sgtindex_size, SEEK_SET);
   //lseek(fp_with_header, sgtmaster_size + sgtindex_size, SEEK_SET);
   free(index_key);

   if (rflag==0) {
        memcpy((char*)*sgtindx, tindx, sgtindex_size);
        rflag = 1;
   } else {
      //Do index check
      for(ip=0;ip<sgtmast->globnp;ip++)
         {
         if(tindx[ip].indx != (*sgtindx)[ip].indx) {
            eflag = -1;
         }
         } 

      if(eflag != 0)
         {
         fprintf(stderr,"sgtindex inconsistency in zfile= %s, exiting ...\n",sgtfpar->zfile);
         exit(-1);
         }
   }
   free(tindx);
   free(tmast);
   
   if(sgtfpar->zfile_header[0] != '\0') {
        sgtfpar->head_off = 0;
   } else { 
        sgtfpar->head_off = (off_t)(sgtmaster_size + sgtindex_size);
   }
   sgtfpar->zcur_off = sgtfpar->head_off;
}

memcached_free(mst);
printf("sgtindx[0].h = %f\n", (*sgtindx)[0].h);
fflush(stdout);
fflush(stderr);
}

void find_sgt(struct sgtparams *sgtpar,struct sgtmaster *sgtmast,struct sgtindex *sgtindx,struct sgtindex *eqindx,struct sgtindex *statindx,float *maxd,float *fwt)
{
float xx, yy, zz, rng, zdp, rexact, zexact, del, delta[4];
float sum, mind, xwt;
int ip, i;
int p0, p1, p2, p3;
int zflag;

/* first see if there is an exact match */

ip = 0;
//while(eqindx->indx > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
//   ip++;
//Do a binary search instead
int start_index = 0;
int end_index = sgtmast->globnp-1;
int mid_index = (start_index + end_index)/2;
int iterations = 0;
while (end_index>=start_index) {
	mid_index = (start_index + end_index)/2;
	/*if (iterations>30) {
		fprintf(stderr, "Error in finding SGTs, aborting.\n");
		if (debug) {
			char buf[512];
			sprintf(buf, "Error looking for %ld. start_index = %d, mid_index = %d, end_index = %d.  start_index.indx = %ld, mid_index.indx = %ld, end_index.indx = %ld.", eqindx->indx, start_index, mid_index, end_index, sgtindx[start_index].indx, sgtindx[mid_index].indx, sgtindx[end_index].indx);
			write_log(buf);
			close_log();
		}
		MPI_Finalize();
		exit(10);
	}*/
	if (eqindx->indx > sgtindx[mid_index].indx) {
		start_index = mid_index+1;
	} else if (eqindx->indx < sgtindx[mid_index].indx) {
		end_index = mid_index-1;
	} else {
		break;
	}
	iterations++;
}
ip = mid_index;

if(eqindx->indx == sgtindx[ip].indx)   /* great, exact match, this will be easy */
   {
   sgtpar->nsgt = 1;
   sgtpar->indx[0] = sgtindx[ip].indx;
   sgtpar->wt[0] = 1.0;

   
//   fprintf(stderr,"SGT EXACT: eqindx= %Ld ip= %d\n",eqindx->indx,ip);
   
   }
else  /* more difficult, find up to 4 SGT that bracket point in range and depth */
   {
   if (debug) {
	char buf[512];
	sprintf(buf, "No exact match found for %ld.\n", eqindx->indx);
	write_log(buf);
   }
   sgtpar->nsgt = 0;
   p0 = -1; p1 = -1; p2 = -1; p3 = -1;

   delta[0] = 1.0e+15;
   delta[1] = 1.0e+15;
   delta[2] = 1.0e+15;
   delta[3] = 1.0e+15;

   xx = (eqindx->xsgt - statindx->xsgt);
   yy = (eqindx->ysgt - statindx->ysgt);
   rexact = xx*xx + yy*yy;
   zexact = (eqindx->zsgt - statindx->zsgt);

   for(ip=0;ip<sgtmast->globnp;ip++)
      {
      xx = (sgtindx[ip].xsgt - statindx->xsgt);
      yy = (sgtindx[ip].ysgt - statindx->ysgt);
      rng = xx*xx + yy*yy;
      zdp = (sgtindx[ip].zsgt - statindx->zsgt);

      xx = (sgtindx[ip].xsgt - eqindx->xsgt);
      yy = (sgtindx[ip].ysgt - eqindx->ysgt);
      zz = (sgtindx[ip].zsgt - eqindx->zsgt);
      del = xx*xx + yy*yy + zz*zz;

      if(rng <= rexact && zdp <= zexact)
         {
	 if(p0 < 0)   /* first time here */
	    {
	    p0 = sgtpar->nsgt;
	    sgtpar->nsgt = sgtpar->nsgt + 1;
	    }

	 if(del < delta[p0])
	    {
	    delta[p0] = del;
	    sgtpar->indx[p0] = sgtindx[ip].indx;
	    sgtpar->wt[p0] = sqrt(del);
	    }
	 }
      else if(rng > rexact && zdp <= zexact)
         {
	 if(p1 < 0)   /* first time here */
	    {
	    p1 = sgtpar->nsgt;
	    sgtpar->nsgt = sgtpar->nsgt + 1;
	    }

	 if(del < delta[p1])
	    {
	    delta[p1] = del;
	    sgtpar->indx[p1] = sgtindx[ip].indx;
	    sgtpar->wt[p1] = sqrt(del);
	    }
	 }
      else if(rng <= rexact && zdp > zexact)
         {
	 if(p2 < 0)   /* first time here */
	    {
	    p2 = sgtpar->nsgt;
	    sgtpar->nsgt = sgtpar->nsgt + 1;
	    }

	 if(del < delta[p2])
	    {
	    delta[p2] = del;
	    sgtpar->indx[p2] = sgtindx[ip].indx;
	    sgtpar->wt[p2] = sqrt(del);
	    }
	 }
      else /* should be (rng > rexact && zdp > zexact)  */
         {
	 if(p3 < 0)   /* first time here */
	    {
	    p3 = sgtpar->nsgt;
	    sgtpar->nsgt = sgtpar->nsgt + 1;
	    }

	 if(del < delta[p3])
	    {
	    delta[p3] = del;
	    sgtpar->indx[p3] = sgtindx[ip].indx;
	    sgtpar->wt[p3] = sqrt(del);
	    }
	 }
      }

   zflag = -1;
   for(i=0;i<sgtpar->nsgt;i++)
      {
      if(sgtpar->wt[i] == 0.0)
         zflag = i;
      }

   if(zflag >= 0)
      {
      for(i=0;i<sgtpar->nsgt;i++)
         sgtpar->wt[i] = 0.0;

      sgtpar->wt[zflag] = 1.0;
      }
   else
      {
      sum = 0.0;
      for(i=0;i<sgtpar->nsgt;i++)
         sum = sum + 1.0/sgtpar->wt[i];

      sum = 1.0/sum;
      for(i=0;i<sgtpar->nsgt;i++)
         sgtpar->wt[i] = sum/sgtpar->wt[i];
      }

   mind = 1.0e+15;
   for(i=0;i<sgtpar->nsgt;i++)
      {
      if(delta[i] < mind)
	 {
         mind = delta[i];
	 xwt = sgtpar->wt[i];
	 }
      }
   if(mind > *maxd)
      {
      *maxd = mind;
      *fwt = xwt;
      }

   /*
   for(i=0;i<sgtpar->nsgt;i++)
      {
      if(delta[i] > *maxd)
         *maxd = delta[i];
      }
      */

   if(sgtpar->nsgt < 4)
      fprintf(stderr,"*** tried to find 4 SGT, but only found %d: eq.zsgt= %d\n",sgtpar->nsgt,eqindx->zsgt);

      
   //   {
/*
for(i=0;i<sgtpar->nsgt;i++)
   {
   fprintf(stderr,"%d) sgti=%Ld eqi=%Ld delta=%13.5e maxd=%13.5e\n",i,sgtpar->indx[i],eqindx->indx,delta[i],*maxd);
   }*/
/*
   exit(-1);
   }
*/

   }
}

void read_part_sgt(struct sgtfileparams*sgtfpar, struct sgtmaster* sgtmast, struct sgtindex *sgtindx, struct sgtheader *sgthead, float* sgtbuf, int starting_pt, int ending_pt) {
	int ip;
	int num_pts_to_read = ending_pt - starting_pt;

	printf("Reading from point %d to point %d.\n", starting_pt, ending_pt);

	float xmom = 0.0;
	float ymom = 0.0;

	//3 components
	memset(sgtbuf, 0, (long long)18*num_pts_to_read*sgtmast->nt*sizeof(float));
	memset(sgthead, 0, num_pts_to_read*sizeof(struct sgtheader));
	int file_entry_size = 6*sgtmast->nt*sizeof(float) + sizeof(struct sgtheader);
	char* rawdata = check_malloc((long long)num_pts_to_read*file_entry_size);
	memset(rawdata, 0, (long long)num_pts_to_read*file_entry_size);

	if (sgtfpar->xfile[0] != '\0') {
		if (sgtfpar->xcur_off==0) {
			//This is our first read
			long loc = fseek(sgtfpar->xfp,(sgtfpar->head_off),SEEK_SET);
		}
		//otherwise, just start where we left off
		freed(sgtfpar->xfp,rawdata,1,((size_t)num_pts_to_read)*file_entry_size);
		for (ip=0; ip<num_pts_to_read; ip++) {
			memcpy(sgthead+ip, rawdata+(long long)ip*file_entry_size, sizeof(struct sgtheader));
			memcpy(sgtbuf+(long long)ip*18*sgtmast->nt, rawdata+(long long)ip*file_entry_size+sizeof(struct sgtheader), 6*sgtmast->nt*sizeof(float));
		}
		sgtfpar->xcur_off = ftell(sgtfpar->xfp);
		xmom = sgthead[0].xmom;
	}
	memset(rawdata, 0, (long long)num_pts_to_read*file_entry_size);
	//y
        if (sgtfpar->yfile[0] != '\0') {
                if (sgtfpar->ycur_off==0) {
			long loc = fseek(sgtfpar->yfp,(sgtfpar->head_off),SEEK_SET);
		}
                freed(sgtfpar->yfp,rawdata,1,((size_t)num_pts_to_read)*file_entry_size);
                for (ip=0; ip<num_pts_to_read; ip++) {
                        memcpy(sgthead+ip, rawdata+(long long)ip*file_entry_size, sizeof(struct sgtheader));
                        memcpy(sgtbuf+(long long)(ip*18+6)*sgtmast->nt, rawdata+(long long)ip*file_entry_size+sizeof(struct sgtheader), 6*sgtmast->nt*sizeof(float));
			sgthead[ip].xmom = xmom;
                }
                sgtfpar->ycur_off = ftell(sgtfpar->yfp);
                ymom = sgthead[0].ymom;
	}
        memset(rawdata, 0, (long long)num_pts_to_read*file_entry_size);
	//z
	if (sgtfpar->zfile[0] != '\0') {
                if (sgtfpar->zcur_off==0) {
                        long loc = fseek(sgtfpar->zfp,(sgtfpar->head_off),SEEK_SET);
                }
                freed(sgtfpar->zfp,rawdata,1,((size_t)num_pts_to_read)*file_entry_size);
                for (ip=0; ip<num_pts_to_read; ip++) {
                        memcpy(sgthead+ip, rawdata+(long long)ip*file_entry_size, sizeof(struct sgtheader));
                        memcpy(sgtbuf+(long long)(ip*18+12)*sgtmast->nt, rawdata+(long long)ip*file_entry_size+sizeof(struct sgtheader), 6*sgtmast->nt*sizeof(float));
        	        sgthead[ip].xmom = xmom;
	                sgthead[ip].ymom = ymom;
                }
                sgtfpar->zcur_off = ftell(sgtfpar->zfp);
        }
	free(rawdata);
}

void read_sgt(struct sgtfileparams *sgtfpar,struct sgtmaster *sgtmast,struct sgtindex *sgtindx,struct sgtheader *sgthead,float *sgtbuf)
{
long long ip;
//float *sgtptr;

float xmom = 0.0;
float ymom = 0.0;

//for(ip=0;ip<18*(sgtmast->globnp)*(sgtmast->nt);ip++)
//   sgtbuf[ip] = 0.0;

memset(sgtbuf, 0, (long long)sgtmast->globnp*sgtmast->nt*sizeof(float));
int full_point_size = 18*sgtmast->nt*sizeof(float) + sizeof(struct sgtheader);
long long file_entry_size = 6*sgtmast->nt*sizeof(float) + sizeof(struct sgtheader);
char* rawdata = check_malloc(sgtmast->globnp*file_entry_size);
memset(rawdata, 0, sgtmast->globnp*file_entry_size);

if(sgtfpar->xfile[0] != '\0')
   {
   long loc = fseek(sgtfpar->xfp,(sgtfpar->head_off),SEEK_SET);
   //long loc = lseek(sgtfpar->xfdr,(sgtfpar->head_off),SEEK_SET);

   freed(sgtfpar->xfp,rawdata,1,((size_t)sgtmast->globnp)*file_entry_size);

   for(ip=0;ip<(sgtmast->globnp);ip++)
      {
      memcpy(sgthead+ip, rawdata+ip*file_entry_size, sizeof(struct sgtheader));
      memcpy(sgtbuf+ip*18*sgtmast->nt, rawdata+ip*file_entry_size+sizeof(struct sgtheader), 6*sgtmast->nt*sizeof(float));
      //sgtptr = sgtbuf + ip*18*(sgtmast->nt);
      //freed(sgtfpar->xfp,&sgthead[ip],sizeof(struct sgtheader), 1);
      //reed(sgtfpar->xfdr,&sgthead[ip],sizeof(struct sgtheader));
      //freed(sgtfpar->xfp,sgtptr,sizeof(float),6*sgtmast->nt);
      //reed(sgtfpar->xfdr,sgtptr,6*(sgtmast->nt)*sizeof(float));
      }
   fclose(sgtfpar->xfp);
   //close(sgtfpar->xfdr);

   xmom = sgthead[0].xmom;
   }

memset(rawdata, 0, sgtmast->globnp*file_entry_size);

if(sgtfpar->yfile[0] != '\0')
   {
   fseek(sgtfpar->yfp,(sgtfpar->head_off),SEEK_SET);
   //lseek(sgtfpar->yfdr,(sgtfpar->head_off),SEEK_SET);

   freed(sgtfpar->yfp,rawdata,1,((size_t)sgtmast->globnp)*file_entry_size);	

   for(ip=0;ip<(sgtmast->globnp);ip++)
      {
      memcpy(sgthead+ip, rawdata+ip*file_entry_size, sizeof(struct sgtheader));
      memcpy(sgtbuf+(ip*18+6)*sgtmast->nt, rawdata+ip*file_entry_size+sizeof(struct sgtheader), 6*sgtmast->nt*sizeof(float));
      //sgtptr = sgtbuf + (ip*18 + 6)*(sgtmast->nt);
      //freed(sgtfpar->yfp,&sgthead[ip],sizeof(struct sgtheader),1);
      //reed(sgtfpar->yfdr,&sgthead[ip],sizeof(struct sgtheader));
      //freed(sgtfpar->yfp,sgtptr,sizeof(float),6*sgtmast->nt);
      //reed(sgtfpar->yfdr,sgtptr,6*(sgtmast->nt)*sizeof(float));

      sgthead[ip].xmom = xmom;
      }
   fclose(sgtfpar->yfp);
   //close(sgtfpar->yfdr);

   ymom = sgthead[0].ymom;
   }

memset(rawdata, 0, sgtmast->globnp*file_entry_size);

if(sgtfpar->zfile[0] != '\0')
   {
   fseek(sgtfpar->zfp,sgtfpar->head_off,SEEK_SET);
   //lseek(sgtfpar->zfdr,(sgtfpar->head_off),SEEK_SET);

   freed(sgtfpar->zfp,rawdata,1,((size_t)sgtmast->globnp)*file_entry_size);

   for(ip=0;ip<(sgtmast->globnp);ip++)
      {
      memcpy(sgthead+ip, rawdata+ip*file_entry_size, sizeof(struct sgtheader));
      memcpy(sgtbuf+(ip*18+12)*sgtmast->nt, rawdata+ip*file_entry_size+sizeof(struct sgtheader), 6*sgtmast->nt*sizeof(float));
      //sgtptr = sgtbuf + (ip*18 + 12)*(sgtmast->nt);
      //freed(sgtfpar->zfp,&sgthead[ip],sizeof(struct sgtheader), 1);
      //reed(sgtfpar->zfdr,&sgthead[ip],sizeof(struct sgtheader));
      //freed(sgtfpar->zfp,sgtptr,sizeof(float),6*sgtmast->nt);
      //reed(sgtfpar->zfdr,sgtptr,6*(sgtmast->nt)*sizeof(float));

      sgthead[ip].xmom = xmom;
      sgthead[ip].ymom = ymom;
      }
   fclose(sgtfpar->zfp);
   //close(sgtfpar->zfdr);
   }
   free(rawdata);
}

void sum_sgt(float *seis,int ntout,float *gfmech,struct sgtparams *sgtpar,int ntsum,float *rupt,float *tstart,struct mechparam mp)
{
int ig, ip, it, im;
float pbar, maxgft, backt0, t0[4], gft[4];
float *sptr, *gfptr;
struct sgtheader* sgtheadptr;

pbar = 0.0;
for(ig=0;ig<sgtpar->nsgt;ig++)
   {
   //ip = sgtpar->master_ip[ig];
   sgtheadptr = sgtpar->sgt_head_ptrs[ig];
   //pbar = pbar + sqrt(sgthead[ip].rho/sgthead[ip].mu);  /* slowness */
   pbar = pbar + sqrt(sgtheadptr->rho/sgtheadptr->mu);
   }
pbar = pbar/(float)(sgtpar->nsgt);

maxgft = -1.0e+15;
for(ig=0;ig<sgtpar->nsgt;ig++)
   {
   //ip = sgtpar->master_ip[ig];
   sgtheadptr = sgtpar->sgt_head_ptrs[ig];
   //gft[ig] = pbar*sgthead[ip].cdist;
   gft[ig] = pbar*sgtheadptr->cdist;
   if(gft[ig] > maxgft)
      maxgft = gft[ig];
   }

backt0 = 0.0;
for(ig=0;ig<sgtpar->nsgt;ig++)
   {
   //ip = sgtpar->master_ip[ig];
   sgtheadptr = sgtpar->sgt_head_ptrs[ig];
   t0[ig] = maxgft - gft[ig];
   backt0 = backt0 - t0[ig]*sgtpar->wt[ig];
   //t0[ig] = t0[ig] + sgthead[ip].tst;
   t0[ig] = t0[ig] + sgtheadptr->tst;
   }
backt0 = backt0 + *rupt - *tstart;

for(im=0;im<mp.nmech;im++)
   {
   sptr = seis + 3*im*ntout;
   gfptr = gfmech + 12*im*ntsum;
   //timeshift_sgt(sptr,ntout,gfptr,sgthead,ntsum,t0,&backt0,sgtpar->nsgt);
   timeshift_sgt(sptr,ntout,gfptr,ntsum,t0,&backt0,sgtpar);
   }
}

//void timeshift_sgt(float *seis,int ntout,float *gf,struct sgtheader *gfh,int ntsum,float *t0,float *bt0,int nsgt)
void timeshift_sgt(float* seis, int ntout, float* gf, int ntsum, float* t0, float* bt0, struct sgtparams* sgtpar)
{
struct complex *gc0, *gc1, *gc2;
float *gf0, *gf1, *gf2, *sv, *sn, *se, *gfv, *gfn, *gfe;
float cosA, sinA, arg, fac, norm, tmpre, scale, tsh;
int i, ig, tapst, it, nts3, nts6, nts9;
int itshift, it0, nf2;

int nsgt = sgtpar->nsgt;
struct sgtheader* sgtheadptr;

int taplen = 10;
float zap = 0.0;
float half = 0.5;
float one = 1.0;
float two = 2.0;
float pi = 3.141592654;

sv = seis;
sn = seis + ntout;
se = seis + 2*ntout;

for(ig=0;ig<nsgt;ig++)
   {
   gf0 = gf + 3*ig*ntsum;
   gf1 = gf + 3*ig*ntsum + ntsum;
   gf2 = gf + 3*ig*ntsum + 2*ntsum;

   /* taper */
   sgtheadptr = sgtpar->sgt_head_ptrs[ig];
   //tapst = gfh[ig].nt - taplen;
   tapst = sgtheadptr->nt - taplen;
   arg = pi/(float)(taplen);
   //for(it=tapst+1;it<gfh[ig].nt;it++)
   for(it=tapst+1;it<sgtheadptr->nt;it++)
      {
      fac = half*(one + cos(arg*(it-tapst)));
      //fprintf(stderr, "gf0[%d] = %e, gf1[%d] = %e, gf2[%d] = %e\n", it, gf0[it], it, gf1[it], it, gf2[it]);
      gf0[it] = fac*gf0[it];
      gf1[it] = fac*gf1[it];
      gf2[it] = fac*gf2[it];
      }

/* apply time shift */

   tsh = t0[ig] + *bt0;
   if(tsh >= 0.0) {
      //itshift = (int)(tsh/gfh[ig].dt + 0.5);
      itshift = (int)(tsh/sgtheadptr->dt + 0.5);
   } else {
      //itshift = (int)(tsh/gfh[ig].dt - 0.5);
      itshift = (int)(tsh/sgtheadptr->dt - 0.5);
   }

   //it0 = gfh[ig].nt + itshift;
   it0 = sgtheadptr->nt + itshift;

   if(it0 > ntsum)
      it0 = ntsum;

   if(it0 > ntout)
      it0 = ntout;

   if(itshift < 0)
      {
      for(i=0;i<it0;i++)
         {
	 sv[i] = sv[i] + gf0[i-itshift];
	 sn[i] = sn[i] + gf1[i-itshift];
	 se[i] = se[i] + gf2[i-itshift];
	 }
      }
   else
      {
      for(i=it0-1;i>=itshift;i--)
         {
	 //fprintf(stderr, "gf0[%d] = %e, gf1[%d] = %e, gf2[%d] = %e\n", i-itshift, gf0[i-itshift], i-itshift, gf1[i-itshift], i-itshift, gf2[i-itshift]);
	 sv[i] = sv[i] + gf0[i-itshift];
	 sn[i] = sn[i] + gf1[i-itshift];
	 se[i] = se[i] + gf2[i-itshift];
	 }
      }
   }
}

void mech_sgt(float* gfmech, struct sgtparams* sgtpar, int nts, struct mechparam mp, float* scl)
//void mech_sgt(float *gfmech,float *sgtbuf,struct sgtheader *sgthead,struct sgtparams *sgtpar,int nts,struct mechparam mp,float *scl)
{
struct sgtheader *sgtheadptr;
float *sgtbufptr;
float *zdd, *rdd, *zds, *rds, *tds, *zss, *rss, *tss;
float *axx, *ayy, *azz, *axy, *axz, *ayz;
float *bxx, *byy, *bzz, *bxy, *bxz, *byz;
float *cxx, *cyy, *czz, *cxy, *cxz, *cyz;
float *gfn, *gfe, *gfv, *gfmptr;
float f1, f2, f3, f4, f5;
float cxS, sxS, cxD, sxD, cx2D, sx2D, cxL, sxL;
float cxT, sxT, cx2T, sx2T;
float arg, cosA, sinA, rad, tan, scale;
float xamp, yamp, zamp, sx, sy, sz;
float mxx, myy, mzz, mxy, mxz, myz;
float sum, rake;
float u1, u2, u3, vx, vy, vz, l2m;
float us, ud, ux, uy, uz;
int it, ig, im;

float half = 0.5;
float two = 2.0;
float rperd = 0.017453293;

arg = (mp.dip)*rperd;
cxD = cos(arg);
sxD = sin(arg);

cx2D = cxD*cxD - sxD*sxD;
sx2D = two*sxD*cxD;

for(im=0;im<mp.nmech;im++)
   {
   u1 = u2 = u3 = 0;
   if(mp.flag[im] == U1FLAG)
      u1 = 1;
   else if(mp.flag[im] == U2FLAG)
      u2 = 1;
   else if(mp.flag[im] == U3FLAG)
      u3 = 1;

   gfmptr = gfmech + im*12*nts;
   zapit(gfmptr,12*nts);

   arg = (mp.rak)*rperd;
   cxL = cos(arg);
   sxL = sin(arg);

   sum = 0.0;
   for(ig=0;ig<(sgtpar->nsgt);ig++)
      {
      //sgtheadptr = sgthead + sgtpar->master_ip[ig];
      sgtheadptr = sgtpar->sgt_head_ptrs[ig];
      //fprintf(stderr, "mech_sgt: head_ptr is at %ld, xmom is %f.\n", sgtheadptr, sgtheadptr->xmom);
      //sgtbufptr = sgtbuf + 18*sgtheadptr->nt*sgtpar->master_ip[ig];
      sgtbufptr = sgtpar->sgtbuf_ptrs[ig];
      //fprintf(stderr, "mech_sgt: buf_ptr is at %ld.\n", sgtbufptr);

      arg = (mp.stk - sgtheadptr->xazim)*rperd;
      cxT = cos(arg);
      sxT = sin(arg);

      vx = -sxD*sxT;
      vy =  sxD*cxT;
      vz = -cxD;

      us = u1*cxL - u2*sxL;
      ud = u1*sxL + u2*cxL;

      ux = -(u3*sxD - ud*cxD)*sxT + us*cxT;
      uy =  (u3*sxD - ud*cxD)*cxT + us*sxT;
      uz = -(u3*cxD + ud*sxD);

      l2m = sgtheadptr->lam + two*sgtheadptr->mu;

      mxx = l2m*vx*ux + (sgtheadptr->lam)*vy*uy + (sgtheadptr->lam)*vz*uz;
      myy = (sgtheadptr->lam)*vx*ux + l2m*vy*uy + (sgtheadptr->lam)*vz*uz;
      mzz = (sgtheadptr->lam)*vx*ux + (sgtheadptr->lam)*vy*uy + l2m*vz*uz;
      mxy = (sgtheadptr->mu)*(vx*uy + vy*ux);
      mxz = (sgtheadptr->mu)*(vx*uz + vz*ux);
      myz = (sgtheadptr->mu)*(vy*uz + vz*uy);

      arg = sgtheadptr->xazim*rperd;
      cosA = cos(arg);
      sinA = sin(arg);

      gfv = gfmptr + 3*ig*nts;
      gfn = gfmptr + 3*ig*nts + nts;
      gfe = gfmptr + 3*ig*nts + 2*nts;

      axx = sgtbufptr;
      ayy = sgtbufptr + (sgtheadptr->nt);
      azz = sgtbufptr + 2*(sgtheadptr->nt);
      axy = sgtbufptr + 3*(sgtheadptr->nt);
      axz = sgtbufptr + 4*(sgtheadptr->nt);
      ayz = sgtbufptr + 5*(sgtheadptr->nt);
      bxx = sgtbufptr + 6*(sgtheadptr->nt);
      byy = sgtbufptr + 7*(sgtheadptr->nt);
      bzz = sgtbufptr + 8*(sgtheadptr->nt);
      bxy = sgtbufptr + 9*(sgtheadptr->nt);
      bxz = sgtbufptr + 10*(sgtheadptr->nt);
      byz = sgtbufptr + 11*(sgtheadptr->nt);
      cxx = sgtbufptr + 12*(sgtheadptr->nt);
      cyy = sgtbufptr + 13*(sgtheadptr->nt);
      czz = sgtbufptr + 14*(sgtheadptr->nt);
      cxy = sgtbufptr + 15*(sgtheadptr->nt);
      cxz = sgtbufptr + 16*(sgtheadptr->nt);
      cyz = sgtbufptr + 17*(sgtheadptr->nt);

      sum = sum + (*scl)*(sgtheadptr->mu)*(sgtpar->wt[ig]);


 	//fprintf(stderr,"area= %13.5e mu= %13.5e\n",(*scl),(sgtheadptr->mu));


      xamp = 0.0;
      if(sgtheadptr->xmom > 0.0) {
         xamp = (*scl)*(sgtpar->wt[ig])/(sgtheadptr->xmom);
	 //fprintf(stderr, "scl=%f, wg=%f, xmom=%f, xamp=%f\n", *scl, sgtpar->wt[ig], sgtheadptr->xmom, xamp);
      }
	 
      yamp = 0.0;
      if(sgtheadptr->ymom > 0.0)
         yamp = (*scl)*(sgtpar->wt[ig])/(sgtheadptr->ymom);
	 
      zamp = 0.0;
      if(sgtheadptr->zmom > 0.0)
         zamp = (*scl)*(sgtpar->wt[ig])/(sgtheadptr->zmom);

      for(it=0;it<sgtheadptr->nt;it++)
         {

         sx = xamp*(axx[it]*mxx + ayy[it]*myy + azz[it]*mzz
               + axy[it]*mxy + axz[it]*mxz + ayz[it]*myz);

         sy = yamp*(bxx[it]*mxx + byy[it]*myy + bzz[it]*mzz
               + bxy[it]*mxy + bxz[it]*mxz + byz[it]*myz);

         sz = zamp*(cxx[it]*mxx + cyy[it]*myy + czz[it]*mzz
               + cxy[it]*mxy + cxz[it]*mxz + cyz[it]*myz);

	 /*if (it<10) {
		 fprintf(stderr, "it=%d, sx=%f, sy=%f, sz=%f\n", it, sx, sy, sz);
		 fprintf(stderr, "axx=%f, mxx=%f, ayy=%f, myy=%f, azz=%f, mzz=%f, axy=%f, mxy=%f, axz=%f, mxz=%f, ayx=%f, myz=%f\n", axx[it], mxx, ayy[it], myy, azz[it], mzz, axy[it], mxy, axz[it], mxz, ayz[it], myz);
		 fprintf(stderr, "xamp=%f\n", xamp);
	 }*/
	
         gfe[it] = sx*sinA + sy*cosA;
         gfn[it] = sx*cosA - sy*sinA;
         gfv[it] = -sz;
         }
      }
   }

*scl = sum;  /* scl now contains the moment released for this point source */
//fprintf(stderr,"sum = %e\n", sum);
}
