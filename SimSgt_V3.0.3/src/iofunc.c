/*
   iofunc.c contains the following functions:
      
      wrout()
      fopfile()
      opfile()
      croptrfile()
      reed()
      rite()
*/

#include "include.h"

/*
      wrout() writes selected output for the velocities vx, vy and vz.
*/

wrout(outp,pvf,medf,iy,glob_iy,nseis,nx,nz,it)
struct outputfields *outp;
float *pvf, *medf;
int iy, glob_iy, nseis, nx, nz, it;
{
struct tsoutputparams *tsoutptr;
struct seisoutputparams *soutptr;
struct seisheader *shptr;
struct sgtoutputparams *sgtoutptr;
struct sgtmaster *sgtmastptr;
struct sgtheader *sgtheadptr;
float *pin, *pout;
int i, j, np, ip;
off_t off, bblen, nrite;
int its, kp;
int size_float = sizeof(float);
float den, num;
float one = 1.0;
float onep5 = 1.5;
float two = 2.0;
float three = 3.0;
float quart = 0.25;
float normf = 1.0e-05;

np = nx*nz;

if(nseis)
   {
   soutptr = &(outp->spnts);

   if(soutptr->iflag)  /* all_in_one output */
      {
      if(soutptr->flushbuf == 1)
         {
         bblen = 3*(soutptr->np)*(soutptr->nbufrite)*size_float;
         off = sizeof(int) + (soutptr->np)*sizeof(struct seisheader) + bblen
                        - soutptr->cur_off;
         lseek(soutptr->fdw,off,SEEK_CUR);

         bblen = 3*(soutptr->np)*(soutptr->buflen)*size_float;
         nrite = rite(soutptr->fdw,soutptr->s,bblen);

         soutptr->cur_off = soutptr->cur_off + off + nrite;
         soutptr->nbufrite = soutptr->nbufrite + soutptr->buflen;
         soutptr->buflen = 0;
         }
      else
         {
         if(it%(soutptr->tinc) == 0)
            {
            shptr = soutptr->shead;
            for(j=0;j<(soutptr->np);j++)
               {
               if(iy == shptr[j].iy)
                  {
                  if(j == 0) /* increment buffer counter */
                     soutptr->buflen = soutptr->buflen + 1;

                  ip = shptr[j].iz*nx + shptr[j].ix;

                  its = it/(soutptr->tinc) - (soutptr->nbufrite);
                  kp = 3*(soutptr->np)*its + 3*j;

                  for(i=0;i<3;i++)
                     soutptr->s[i+kp] = pvf[i*np+ip];
                  }
               }
            }
         }
      }
   }   

tsoutptr = &(outp->xyslice);
if(tsoutptr->iflag)
   {
   if(tsoutptr->flushbuf == 1 && tsoutptr->buflen)
      {
      bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*size_float;
      off = bblen + tsoutptr->head_off - tsoutptr->cur_off;
      lseek(tsoutptr->fdw,off,SEEK_CUR);

      bblen = 3*(tsoutptr->np)*(tsoutptr->buflen)*size_float;
      nrite = rite(tsoutptr->fdw,tsoutptr->vbuf,bblen);

      tsoutptr->cur_off = tsoutptr->cur_off + off + nrite;
      tsoutptr->nbufrite = tsoutptr->nbufrite + tsoutptr->buflen;
      tsoutptr->buflen = 0;
      }
   else
      {
      if(glob_iy >= tsoutptr->iyleft && glob_iy <= tsoutptr->iyright && glob_iy%(tsoutptr->idy) == 0 && it%(tsoutptr->idt) == 0 && it >= 0)
         {
	 its = it/(tsoutptr->idt) - (tsoutptr->nbufrite);
	 tsoutptr->buflen = 1 + its;

         for(i=0;i<3;i++)
            {
            pin = pvf + i*np + tsoutptr->idz*nx;
	    pout = tsoutptr->vbuf + (3*its + i)*(tsoutptr->np) + tsoutptr->nx*((glob_iy - tsoutptr->iyleft)/tsoutptr->idy);

            copy(pin,pout,tsoutptr->nx,tsoutptr->idx);
            }
         }
      }
   }

tsoutptr = &(outp->xzslice);
if(tsoutptr->iflag)
   {
   if(tsoutptr->flushbuf == 1 && tsoutptr->buflen)
      {
      bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*size_float;
      off = bblen + tsoutptr->head_off - tsoutptr->cur_off;
      lseek(tsoutptr->fdw,off,SEEK_CUR);

      bblen = 3*(tsoutptr->np)*(tsoutptr->buflen)*size_float;
      nrite = rite(tsoutptr->fdw,tsoutptr->vbuf,bblen);

      tsoutptr->cur_off = tsoutptr->cur_off + off + nrite;
      tsoutptr->nbufrite = tsoutptr->nbufrite + tsoutptr->buflen;
      tsoutptr->buflen = 0;
      }
   else
      {
      if(glob_iy == tsoutptr->idy && it%(tsoutptr->idt) == 0 && it >= 0)
         {
         its = it/(tsoutptr->idt) - (tsoutptr->nbufrite);
         tsoutptr->buflen = 1 + its;

         for(i=0;i<3;i++)
            {
            for(j=0;j<tsoutptr->nz;j++)
               {
               pin = pvf + i*np + (tsoutptr->iz0 + j*tsoutptr->idz)*nx;
               pout = tsoutptr->vbuf + (3*its + i)*(tsoutptr->np) + j*tsoutptr->nx;
               copy(pin,pout,tsoutptr->nx,tsoutptr->idx);
               }
            }
         }
      }
   }

tsoutptr = &(outp->yzslice);
if(tsoutptr->iflag)
   {
   if(tsoutptr->flushbuf == 1 && tsoutptr->buflen)
      {
      bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*size_float;
      off = bblen + tsoutptr->head_off - tsoutptr->cur_off;
      lseek(tsoutptr->fdw,off,SEEK_CUR);

      bblen = 3*(tsoutptr->np)*(tsoutptr->buflen)*size_float;
      nrite = rite(tsoutptr->fdw,tsoutptr->vbuf,bblen);

      tsoutptr->cur_off = tsoutptr->cur_off + off + nrite;
      tsoutptr->nbufrite = tsoutptr->nbufrite + tsoutptr->buflen;
      tsoutptr->buflen = 0;
      }
   else
      {
      if(glob_iy >= tsoutptr->iyleft && glob_iy <= tsoutptr->iyright && glob_iy%(tsoutptr->idy) == 0 && it%(tsoutptr->idt) == 0 && it >= 0)
         {
         its = it/(tsoutptr->idt) - (tsoutptr->nbufrite);
         tsoutptr->buflen = 1 + its;

         for(i=0;i<3;i++)
            {
            pin = pvf + i*np + tsoutptr->idx;
            pout = tsoutptr->vbuf + (3*its + i)*(tsoutptr->np) + tsoutptr->nz*((glob_iy - tsoutptr->iyleft)/tsoutptr->idy);

            copy(pin,pout,tsoutptr->nz,tsoutptr->idz*nx);
            }
         }
      }
   }

sgtoutptr = &(outp->sgtpnts);
if(sgtoutptr->iflag)
   {
   sgtmastptr = &(sgtoutptr->sgtmast);
   if(sgtoutptr->flushbuf == 1)
      {
      bblen = 6*(sgtmastptr->localnp)*(sgtoutptr->nbufrite)*size_float;
      off = sizeof(struct sgtmaster) + (sgtmastptr->globnp)*sizeof(struct sgtindex)
             + (sgtmastptr->localnp)*sizeof(struct sgtheader) + bblen - sgtoutptr->cur_off;
      lseek(sgtoutptr->fdw,off,SEEK_CUR);

      bblen = 6*(sgtmastptr->localnp)*(sgtoutptr->buflen)*size_float;
      nrite = rite(sgtoutptr->fdw,sgtoutptr->sgt,bblen);

      sgtoutptr->cur_off = sgtoutptr->cur_off + off + nrite;
      sgtoutptr->nbufrite = sgtoutptr->nbufrite + sgtoutptr->buflen;
      sgtoutptr->buflen = 0;
      }
   else
      {
      if(it%(sgtoutptr->tinc) == 0)
         {
         sgtheadptr = sgtoutptr->sgthead;
         for(j=0;j<(sgtmastptr->localnp);j++)
            {
            if(glob_iy == sgtheadptr[j].ysgt)  /* remember, ysgt is global index */
               {
	       if(j == 0) /* increment buffer counter */
                  sgtoutptr->buflen = sgtoutptr->buflen + 1;

	       ip = sgtheadptr[j].zsgt*nx + sgtheadptr[j].xsgt;

	       num = two*(medf[11*np+ip] + medf[12*np+ip])/medf[11*np+ip];
	       den = normf*quart*medf[11*np+ip]/(medf[12*np+ip]*(onep5*medf[11*np+ip] + medf[12*np+ip]));

	       if(it == 0) /* set lambda,rigidity,density and rewrite header */
	          {
	          sgtheadptr[j].lam = (1.0e+10)*medf[11*np+ip];
	          sgtheadptr[j].mu = (1.0e+10)*medf[12*np+ip];
	          sgtheadptr[j].rho = three/(medf[5*np+ip] + medf[6*np+ip] + medf[7*np+ip]);

                  off = sizeof(struct sgtmaster) + (sgtmastptr->globnp)*sizeof(struct sgtindex)
		         + j*sizeof(struct sgtheader) - sgtoutptr->cur_off;
                  lseek(sgtoutptr->fdw,off,SEEK_CUR);

	          nrite = rite(sgtoutptr->fdw,&sgtheadptr[j],sizeof(struct sgtheader));
                  sgtoutptr->cur_off = sgtoutptr->cur_off + off + nrite;
	          }

	       its = it/(sgtoutptr->tinc) - (sgtoutptr->nbufrite);
               kp = 6*(sgtmastptr->localnp)*its + 6*j;

/*
            Note that the strains are NOT normalized by the factor
            of (h).  This factor would balance the (1/h) factor
	    of the body force used to specify the moment.  Since these
	    factors cancel, neither is applied.  See main.c
*/

	       /* mxx */
	       sgtoutptr->sgt[kp]   = den*(num*pvf[3*np+ip]
				        - pvf[4*np+ip] - pvf[5*np+ip]);
	       /* myy */
	       sgtoutptr->sgt[1+kp] = den*(num*pvf[4*np+ip]
				        - pvf[3*np+ip] - pvf[5*np+ip]);
	       /* mzz */
	       sgtoutptr->sgt[2+kp] = den*(num*pvf[5*np+ip]
				        - pvf[3*np+ip] - pvf[4*np+ip]);

	       /* mxy */
	       sgtoutptr->sgt[3+kp] = normf*pvf[6*np+ip]/medf[12*np+ip];
	       /* mxz */
	       sgtoutptr->sgt[4+kp] = normf*pvf[7*np+ip]/medf[12*np+ip];
	       /* myz */
	       sgtoutptr->sgt[5+kp] = normf*pvf[8*np+ip]/medf[12*np+ip];
               }
            }
         }
      }
   }
}

FILE *fopfile(name,mode)
char *name, *mode;
{
FILE *fp, *fopen();
int j;
char path[512];

if((fp = fopen(name,mode)) == NULL && errno == ENOENT)
   {
   j = strlen(name);
   while(name[j] != '/' && j > 1)
      j--;

   strncpy(path,name,j);
   path[j] = '\0';
   makedir(path);
   }

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s ERRNO= %d\n", name, mode,errno);
   perror(NULL);
   exit(-1);
   }
return(fp);
}

FILE *fopdirfile(name,mode)
char *name, *mode;
{
FILE *fp, *fopen();
int j;
char path[512];

if((fp = fopen(name,mode)) == NULL && errno == ENOENT)
   {
   j = strlen(name);
   while(name[j] != '/' && j > 1)
      j--;

   strncpy(path,name,j);
   path[j] = '\0';
   makedir(path);
   }

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

opfile_ro (name)
char *name;
{
int fd;
if ((fd = open (name, RDONLY_FLAGS, 0444)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

opfile (name)
char *name;
{
int fd;
if ((fd = open (name, RDWR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int cropfile_rw(char *name)
{
int fd;
if ((fd = open(name,CROP_FLAGS,0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return(fd);
}

croptrfile (name)
char *name;
{
int fd;
if ((fd = open (name, CROPTR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

size_t reed(int fd, void *pntr, size_t length)
{
size_t temp;
if ((temp = read(fd, pntr, length)) < length)
   {
   fprintf (stderr, "READ ERROR\n");
   fprintf (stderr, "%Ld attempted  %Ld read\n", length, temp);
   exit(-1);
   }
return(temp);
}

size_t rite(int fd, void *pntr, size_t length)
{
size_t temp;
if ((temp = write(fd, pntr, length)) < length)
   {
   fprintf (stderr, "WRITE ERROR\n");
   fprintf (stderr, "%Ld attempted  %Ld written\n", length, temp);
   exit(-1);
   }
return(temp);
}

void init_outp(struct outputfields *outp,struct runparams *rpars,char *name,struct pntsrcs *srcs,struct restartinfo *rst,struct dump_output_info *di)
{
FILE *fpr;
struct tsoutputparams *tsoutptr;
struct tsheader_proc tshead_p;
struct tsheader tshead;
struct seisoutputparams *soutptr;
struct seisheader *sheadptr;
struct sgtoutputparams *sgtoutptr;
struct sgtmaster *sgtmastptr, *rst_sgtmast;
struct sgtindex *sgtindxptr;
struct sgtheader *sgtheadptr;
char string[2048], reed_file[2048];

int k, it, iy, yend, i, iseis, iyleft, iyright, npad;
int xsgt, ysgt, zsgt;
float lonsgt, latsgt, depsgt;
float xx, yy;
int npread;
off_t off, bblen;
float *initval;

makedir(outp->seisdir);

if(outp->nseis)
   {
   fprintf(stderr,"**** Individual station locations:\n");
   fflush(stderr);

   iyleft = rpars->ny1 + 2;
   if(rpars->ny1 == 0) iyleft = 0;

   iyright = rpars->ny2 - 3;
   if(rpars->ny2 == rpars->globny) iyright = rpars->ny2 - 1;

   fpr = fopfile(outp->seiscords,"r");
   fscanf(fpr,"%d",&outp->nseis);

   soutptr = &(outp->spnts);

   soutptr->iflag = 1;  /* backward compatible with all_in_one=1 */

   if(soutptr->iflag == 1)
      {
      soutptr->buflen = 0;
      soutptr->nbufrite = 0;
      soutptr->flushbuf = 0;

      soutptr->ntout = 1 + (((int)((rpars->nt-1)/rpars->span) + 1)*rpars->span - 1)/(soutptr->tinc);

      soutptr->shead = (struct seisheader *) check_malloc (outp->nseis*sizeof(struct seisheader));

      sheadptr = soutptr->shead;
      iseis = 0;
      for(i=0;i<outp->nseis;i++)
         {
         fscanf(fpr,"%d %d %d %s",&(sheadptr[iseis].ix),
                                  &(sheadptr[iseis].iy),
                                  &(sheadptr[iseis].iz),
                                  string);

	 if(sheadptr[iseis].iy >= iyleft && sheadptr[iseis].iy <= iyright)
	    {
            strncpy(sheadptr[iseis].name,string,NBYTE_STATNAME-1);
            sheadptr[iseis].name[NBYTE_STATNAME-1] = '\0';

            sheadptr[iseis].indx = i;
            sheadptr[iseis].nt = soutptr->ntout;
            sheadptr[iseis].dt = rpars->dt*soutptr->tinc;
            sheadptr[iseis].h = rpars->h;
            sheadptr[iseis].modelrot = rpars->modelrot;

	    xx = sheadptr[iseis].ix*rpars->h;
	    yy = sheadptr[iseis].iy*rpars->h;

	    gcproj(&xx,&yy,&(sheadptr[iseis].slon),&(sheadptr[iseis].slat),&rpars->erad,&rpars->g0,&rpars->b0,rpars->amat,rpars->ainv,0);

            iseis++;
	    }
         }
      fclose(fpr);

      outp->nseis = iseis;

      if(outp->nseis)
         {
         soutptr->np = outp->nseis;

         soutptr->shead = (struct seisheader *) check_realloc (soutptr->shead,(soutptr->np)*sizeof( struct seisheader));
         soutptr->s = (float *) check_malloc (3*(soutptr->np)*rpars->span*sizeof(float));

         sprintf(soutptr->local_file,"%s/%s_seis-%.5d.e3d",outp->seisdir,name,rpars->segmentId);
         sprintf(soutptr->main_file,"%s/%s_seis-%.5d.e3d",di->main_dir,di->name,rpars->segmentId);

	 if(rst->read_flag == 0)
            soutptr->fdw = croptrfile(soutptr->local_file);
	 else
            {
            fcp_read_write(soutptr->main_file,soutptr->local_file);
            soutptr->fdw = opfile(soutptr->local_file);
            }

         fprintf(stderr,"     ALL seismograms written into single output file\n\n");

         fprintf(stderr,"     local file= %s\n",soutptr->local_file);
         fprintf(stderr,"     main file= %s\n",soutptr->main_file);
         fprintf(stderr,"     nseis= %d\n",soutptr->np);
         fprintf(stderr,"     ntout= %d\n",soutptr->ntout);
         fprintf(stderr,"     dtout= %f\n",rpars->dt*soutptr->tinc);
         fprintf(stderr,"     seismogram output buffer is %.2f Mbytes\n",(3.0*(soutptr->np)*rpars->span*sizeof(float))/1.0e+06);
         fprintf(stderr,"     Output file is %.2f Mbytes\n",(soutptr->np)*(sizeof(struct seisheader) + 3.0*(soutptr->ntout)*sizeof(float))/1.0e+06);
         fflush(stderr);

	 if(rst->read_flag == 0)
            rite(soutptr->fdw,&(soutptr->np),sizeof(int));
	 else
	    {
            reed(soutptr->fdw,&npread,sizeof(int));
	    if(npread != soutptr->np)
	       {
	       fprintf(stderr,"***** npread (%d) != soutptr->np (%d), file:\n",npread,soutptr->np);
	       fprintf(stderr,"      %s,\t\texiting ...\n",soutptr->local_file);
	       exit(-1);
	       }
	    }

         sheadptr = soutptr->shead;
         for(i=0;i<(soutptr->np);i++)
	    {
	    /* re-write headers just in case nt has changed */
            rite(soutptr->fdw,&sheadptr[i],sizeof(struct seisheader));

	    /* now adjust header iy into local system */
	    sheadptr[i].iy = sheadptr[i].iy - rpars->ny1;
	    }

	 if(rst->read_flag == 1)
	    {
            bblen = sizeof(int) + (off_t)((soutptr->np)*(sizeof(struct seisheader) + 3*(soutptr->ntout)*sizeof(float)));
            lseek(soutptr->fdw,0,SEEK_SET);
            off = lseek(soutptr->fdw,0,SEEK_END);

	    npad = (bblen - off)/(off_t)(3*sizeof(float)*(soutptr->np));
	    if(npad > 0)
	       {
               initval = (float *) check_malloc (3*npad*sizeof(float));
               for(i=0;i<3*npad;i++)
                  initval[i] = 0;

               for(i=0;i<(soutptr->np);i++)
                  rite(soutptr->fdw,initval,3*npad*sizeof(float));

	       free(initval);
	       }

	    for(it=0;it<=rst->it;it=it+soutptr->tinc)
	       soutptr->nbufrite = soutptr->nbufrite + 1;

	    bblen = 3*(soutptr->np)*(soutptr->nbufrite)*sizeof(float);
            soutptr->cur_off = sizeof(int) + (soutptr->np)*sizeof(struct seisheader) + bblen;
	    }
	 else
	    {
            initval = (float *) check_malloc (3*(soutptr->ntout)*sizeof(float));
            for(i=0;i<3*(soutptr->ntout);i++)
               initval[i] = 0;

            for(i=0;i<(soutptr->np);i++)
               rite(soutptr->fdw,initval,3*(soutptr->ntout)*sizeof(float));

	    free(initval);

            soutptr->cur_off = 0;
	    }

         lseek(soutptr->fdw,soutptr->cur_off,SEEK_SET);
	 }
      else
         {
         soutptr->np = outp->nseis;
         soutptr->iflag = 0;
         free(soutptr->shead);

         fprintf(stderr,"     NO seismograms written from this node\n");
         fflush(stderr);
	 }
      }
   fprintf(stderr,"\n");
   fflush(stderr);

   if(di->enable_flag && soutptr->iflag)
      strcpy(di->local_outputdir,outp->seisdir);
   }

if(outp->ts_xy || outp->ts_xz || outp->ts_yz)
   {
   fprintf(stderr,"**** Time slices:\n");
   fflush(stderr);
   }

tsoutptr = &(outp->xyslice);
if(outp->ts_xy)
   {
   fprintf(stderr,"     xy-plane time slice:\n");
   fflush(stderr);

   tsoutptr->iflag = 1;

   iyleft = rpars->ny1 + 2;
   if(rpars->ny1 == 0) iyleft = 0;

   iyright = rpars->ny2 - 3;
   if(rpars->ny2 == rpars->globny) iyright = rpars->ny2 - 1;

   tshead_p.iyleft = -1;
   tshead_p.localny = 0;
   for(iy=0;iy<rpars->globny;iy++)
      {
      if(iy >= iyleft && iy <= iyright && iy%(outp->dyts) == 0)
         {
	 if(tshead_p.iyleft < 0)
	    tshead_p.iyleft = iy;

	 tshead_p.iyright = iy;

         tshead_p.localny = tshead_p.localny + 1;
	 }
      }

   tshead_p.ix0 = 0;
   tshead_p.iy0 = 0;
   tshead_p.iz0 = outp->iz_ts;
   tshead_p.it0 = 0;
   tshead_p.nx = (int)(((rpars->nx-1)/(outp->dxts))+1);
   tshead_p.ny = (int)(((rpars->globny-1)/(outp->dyts))+1);
   tshead_p.nz = 1;
   tshead_p.nt = (int)(((rpars->nt-1)/(outp->dtts))+1);
   tshead_p.dx = rpars->h*outp->dxts;
   tshead_p.dy = rpars->h*outp->dyts;
   tshead_p.dz = rpars->h;
   tshead_p.dt = rpars->dt*outp->dtts;
   tshead_p.modelrot = rpars->modelrot;
   tshead_p.modellat = rpars->modellat;
   tshead_p.modellon = rpars->modellon;

   tsoutptr->ix0 = 0;
   tsoutptr->iy0 = 0;
   tsoutptr->iz0 = outp->iz_ts;
   tsoutptr->it0 = 0;
   tsoutptr->idx = outp->dxts;
   tsoutptr->idy = outp->dyts;
   tsoutptr->idz = outp->iz_ts;
   tsoutptr->idt = outp->dtts;
   tsoutptr->nx = (int)(((rpars->nx-1)/(tsoutptr->idx))+1);
   tsoutptr->ny = tshead_p.localny;
   tsoutptr->np = (tsoutptr->nx)*(tsoutptr->ny);
   tsoutptr->ntout = tshead_p.nt;

   tsoutptr->iyleft = iyleft;
   tsoutptr->iyright = iyright;

   tsoutptr->buflen = 0;
   tsoutptr->nbufrite = 0;
   tsoutptr->flushbuf = 0;
   tsoutptr->head_off = sizeof(struct tsheader_proc);

   tsoutptr->vbuf = (float *) check_malloc (3*(tsoutptr->np)*(rpars->span)*sizeof(float));

   sprintf(tsoutptr->local_file,"%s/%s_xyts-%.5d.e3d",outp->seisdir,name,rpars->segmentId);
   sprintf(tsoutptr->main_file,"%s/%s_xyts-%.5d.e3d",di->main_dir,di->name,rpars->segmentId);

   if(rst->read_flag == 0)
      {
      tsoutptr->fdw = croptrfile(tsoutptr->local_file);

      rite(tsoutptr->fdw,&tshead_p,sizeof(struct tsheader_proc));

      initval = (float *) check_malloc (3*tshead_p.nx*sizeof(float));
      for(iy=0;iy<3*tshead_p.nx;iy++)
         initval[iy] = 0;

      for(it=0;it<tshead_p.nt*tshead_p.localny;it++)
         rite(tsoutptr->fdw,initval,3*tshead_p.nx*sizeof(float));

      free(initval);

      tsoutptr->cur_off = sizeof(struct tsheader_proc);
      }
   else
      {
      fcp_read_write(tsoutptr->main_file,tsoutptr->local_file);

      tsoutptr->fdw = opfile(tsoutptr->local_file);

      /* re-write header just in case nt has changed */
      rite(tsoutptr->fdw,&tshead_p,sizeof(struct tsheader_proc));

      if(tshead_p.localny > 0)
         {
         bblen = (off_t)(3*tsoutptr->ntout*tsoutptr->np*sizeof(float)) + (off_t)(sizeof(struct tsheader_proc));
         lseek(tsoutptr->fdw,0,SEEK_SET);
         off = lseek(tsoutptr->fdw,0,SEEK_END);

         npad = (bblen - off)/(off_t)(3*sizeof(float)*tsoutptr->np);
         if(npad > 0)
            {
            initval = (float *) check_malloc (3*tshead_p.nx*sizeof(float));
            for(i=0;i<3*tshead_p.nx;i++)
               initval[i] = 0;

            for(it=0;it<npad*tshead_p.localny;it++)
               rite(tsoutptr->fdw,initval,3*tshead_p.nx*sizeof(float));

            free(initval);
            }

         for(it=0;it<=rst->it;it=it+tsoutptr->idt)
            tsoutptr->nbufrite = tsoutptr->nbufrite + 1;

         bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*sizeof(float);
         tsoutptr->cur_off = sizeof(struct tsheader_proc) + bblen;
	 }
      }

   lseek(tsoutptr->fdw,tsoutptr->cur_off,SEEK_SET);

   fprintf(stderr,"   local file= %s\n",tsoutptr->local_file);
   fprintf(stderr,"    main file= %s\n",tsoutptr->main_file);
   fprintf(stderr,"         nxts= %d\n",tshead_p.nx);
   fprintf(stderr,"         nyts= %d\n",tshead_p.localny);
   fprintf(stderr,"         ntts= %d\n",tshead_p.nt);
   fprintf(stderr,"         dtts= %d\n",outp->dtts);
   fprintf(stderr,"     xy-plane time slice output buffer is %.2f Mbytes\n",(3.0*(tsoutptr->np)*rpars->span*sizeof(float))/1.0e+06);
   fprintf(stderr,"     Output file is %.2f Mbytes\n\n",(3.0*tshead_p.nt*tshead_p.localny*tshead_p.nx*sizeof(float) + sizeof(struct tsheader_proc))/1.0e+06);
   fflush(stderr);

   if(di->enable_flag)
      strcpy(di->local_outputdir,outp->seisdir);
   }
else
   tsoutptr->iflag = 0;

tsoutptr = &(outp->xzslice);
if(outp->ts_xz)
   {
   iyleft = rpars->ny1 + 2;
   if(rpars->ny1 == 0) iyleft = 0;

   iyright = rpars->ny2 - 3;
   if(rpars->ny2 == rpars->globny) iyright = rpars->ny2 - 1;

   tsoutptr->iflag = 0;
   if(outp->iy_ts >= iyleft && outp->iy_ts <= iyright)
      {
      fprintf(stderr,"     xz-plane time slice:\n");
      fflush(stderr);

      tsoutptr->iflag = 1;

      tshead.ix0 = 0;
      tshead.iy0 = outp->iy_ts;
      tshead.iz0 = 0;
      tshead.it0 = 0;
      tshead.nx = (int)(((rpars->nx-1)/(outp->dxts))+1);
      tshead.ny = 1;
      tshead.nz = (int)(((rpars->nz-1)/(outp->dzts))+1);
      tshead.nt = (int)(((rpars->nt-1)/(outp->dtts))+1);
      tshead.dx = rpars->h*outp->dxts;
      tshead.dy = rpars->h;
      tshead.dz = rpars->h*outp->dxts;
      tshead.dt = rpars->dt*outp->dtts;
      tshead.modelrot = rpars->modelrot;
      tshead.modellat = rpars->modellat;
      tshead.modellon = rpars->modellon;

      tsoutptr->ix0 = 0;
      tsoutptr->iy0 = outp->iy_ts;
      tsoutptr->iz0 = 0;
      if(rpars->freesurf == 1)
         tsoutptr->iz0 = 1;

      tsoutptr->it0 = 0;
      tsoutptr->idx = outp->dxts;
      tsoutptr->idy = outp->iy_ts;
      tsoutptr->idz = outp->dzts;
      tsoutptr->idt = outp->dtts;
      tsoutptr->nx = (int)(((rpars->nx-1)/(tsoutptr->idx))+1);
      tsoutptr->ny = 1;
      tsoutptr->nz = (int)(((rpars->nz-1)/(tsoutptr->idz))+1);
      tsoutptr->np = (tsoutptr->nx)*(tsoutptr->nz);
      tsoutptr->ntout = tshead.nt;

      tsoutptr->iyleft = iyleft;
      tsoutptr->iyright = iyright;

      tsoutptr->buflen = 0;
      tsoutptr->nbufrite = 0;
      tsoutptr->flushbuf = 0;
      tsoutptr->head_off = sizeof(struct tsheader);

      tsoutptr->vbuf = (float *) check_malloc (3*(tsoutptr->np)*(rpars->span)*sizeof(float));

      sprintf(tsoutptr->local_file,"%s/%s_xzts.e3d",outp->seisdir,name);
      sprintf(tsoutptr->main_file,"%s/%s_xzts.e3d",di->main_dir,di->name);

      if(rst->read_flag == 0)
         {
         tsoutptr->fdw = croptrfile(tsoutptr->local_file);

         rite(tsoutptr->fdw,&tshead,sizeof(struct tsheader));

         initval = (float *) check_malloc (3*tshead.nx*sizeof(float));
         for(iy=0;iy<3*tshead.nx;iy++)
            initval[iy] = 0;

         for(it=0;it<tshead.nt*tshead.nz;it++)
            rite(tsoutptr->fdw,initval,3*tshead.nx*sizeof(float));

         free(initval);

         tsoutptr->cur_off = sizeof(struct tsheader);
         }
      else
         {
         fcp_read_write(tsoutptr->main_file,tsoutptr->local_file);

         tsoutptr->fdw = opfile(tsoutptr->local_file);

      /* re-write header just in case nt has changed */
         rite(tsoutptr->fdw,&tshead,sizeof(struct tsheader));

         bblen = (off_t)(3*tsoutptr->ntout*tsoutptr->np*sizeof(float)) + (off_t)(sizeof(struct tsheader));
         lseek(tsoutptr->fdw,0,SEEK_SET);
         off = lseek(tsoutptr->fdw,0,SEEK_END);

         npad = (bblen - off)/(off_t)(3*sizeof(float)*tsoutptr->np);
         if(npad > 0)
            {
            initval = (float *) check_malloc (3*tshead.nx*sizeof(float));
            for(i=0;i<3*tshead.nx;i++)
               initval[i] = 0;

            for(it=0;it<npad*tshead.nz;it++)
               rite(tsoutptr->fdw,initval,3*tshead.nx*sizeof(float));

            free(initval);
            }

         for(it=0;it<=rst->it;it=it+tsoutptr->idt)
            tsoutptr->nbufrite = tsoutptr->nbufrite + 1;

         bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*sizeof(float);
         tsoutptr->cur_off = sizeof(struct tsheader) + bblen;
         }

      lseek(tsoutptr->fdw,tsoutptr->cur_off,SEEK_SET);

      fprintf(stderr,"   local file= %s\n",tsoutptr->local_file);
      fprintf(stderr,"    main file= %s\n",tsoutptr->main_file);
      fprintf(stderr,"         nxts= %d\n",tshead.nx);
      fprintf(stderr,"         nzts= %d\n",tshead.nz);
      fprintf(stderr,"         ntts= %d\n",tshead.nt);
      fprintf(stderr,"         dtts= %d\n",outp->dtts);
      fprintf(stderr,"     xz-plane time slice output buffer is %.2f Mbytes\n",(3.0*(tsoutptr->np)*rpars->span*sizeof(float))/1.0e+06);
      fprintf(stderr,"     Output file is %.2f Mbytes\n\n",(3.0*tshead.nt*tshead.nx*sizeof(float) + sizeof(struct tsheader))/1.0e+06);
      fflush(stderr);

      if(di->enable_flag)
         strcpy(di->local_outputdir,outp->seisdir);
      }
   }
else
   tsoutptr->iflag = 0;

tsoutptr = &(outp->yzslice);
if(outp->ts_yz)
   {
   fprintf(stderr,"     yz-plane time slice:\n");
   fflush(stderr);

   tsoutptr->iflag = 1;

   iyleft = rpars->ny1 + 2;
   if(rpars->ny1 == 0) iyleft = 0;

   iyright = rpars->ny2 - 3;
   if(rpars->ny2 == rpars->globny) iyright = rpars->ny2 - 1;

   tshead_p.iyleft = -1;
   tshead_p.localny = 0;
   for(iy=0;iy<rpars->globny;iy++)
      {
      if(iy >= iyleft && iy <= iyright && iy%(outp->dyts) == 0)
         {
         if(tshead_p.iyleft < 0)
            tshead_p.iyleft = iy;

         tshead_p.iyright = iy;

         tshead_p.localny = tshead_p.localny + 1;
         }
      }

   tshead_p.ix0 = outp->ix_ts;
   tshead_p.iy0 = 0;
   tshead_p.iz0 = 0;
   tshead_p.it0 = 0;
   tshead_p.nx = 1;
   tshead_p.ny = (int)(((rpars->globny-1)/(outp->dyts))+1);
   tshead_p.nz = (int)(((rpars->nz-1)/(outp->dzts))+1);
   tshead_p.nx = tshead_p.nz;  /* to make merge_ts work 2011-11-09 */
   tshead_p.nt = (int)(((rpars->nt-1)/(outp->dtts))+1);
   tshead_p.dx = rpars->h;
   tshead_p.dy = rpars->h*outp->dyts;
   tshead_p.dz = rpars->h*outp->dxts;
   tshead_p.dt = rpars->dt*outp->dtts;
   tshead_p.modelrot = rpars->modelrot;
   tshead_p.modellat = rpars->modellat;
   tshead_p.modellon = rpars->modellon;

   tsoutptr->ix0 = outp->ix_ts;
   tsoutptr->iy0 = 0;
   tsoutptr->iz0 = 0;
   tsoutptr->it0 = 0;
   tsoutptr->idx = outp->ix_ts;
   tsoutptr->idy = outp->dyts;
   tsoutptr->idz = outp->dzts;
   tsoutptr->idt = outp->dtts;
   tsoutptr->nz = (int)(((rpars->nz-1)/(tsoutptr->idz))+1);
   tsoutptr->ny = tshead_p.localny;
   tsoutptr->np = (tsoutptr->nz)*(tsoutptr->ny);
   tsoutptr->ntout = tshead_p.nt;

   tsoutptr->iyleft = iyleft;
   tsoutptr->iyright = iyright;

   tsoutptr->buflen = 0;
   tsoutptr->nbufrite = 0;
   tsoutptr->flushbuf = 0;
   tsoutptr->head_off = sizeof(struct tsheader_proc);

   tsoutptr->vbuf = (float *) check_malloc (3*(tsoutptr->np)*(rpars->span)*sizeof(float));

   sprintf(tsoutptr->local_file,"%s/%s_yzts-%.4d.e3d",outp->seisdir,name,rpars->segmentId);
   sprintf(tsoutptr->main_file,"%s/%s_yzts-%.4d.e3d",di->main_dir,di->name,rpars->segmentId);

   if(rst->read_flag == 0)
      {
      tsoutptr->fdw = croptrfile(tsoutptr->local_file);

      rite(tsoutptr->fdw,&tshead_p,sizeof(struct tsheader_proc));

      initval = (float *) check_malloc (3*tshead_p.nz*sizeof(float));
      for(iy=0;iy<3*tshead_p.nz;iy++)
         initval[iy] = 0;

      for(it=0;it<tshead_p.nt*tshead_p.localny;it++)
         rite(tsoutptr->fdw,initval,3*tshead_p.nz*sizeof(float));

      free(initval);

      tsoutptr->cur_off = sizeof(struct tsheader_proc);
      }
   else
      {
      fcp_read_write(tsoutptr->main_file,tsoutptr->local_file);

      tsoutptr->fdw = opfile(tsoutptr->local_file);

      /* re-write header just in case nt has changed */
      rite(tsoutptr->fdw,&tshead_p,sizeof(struct tsheader_proc));

      if(tshead_p.localny > 0)
         {
         bblen = (off_t)(3*tsoutptr->ntout*tsoutptr->np*sizeof(float)) + (off_t)(sizeof(struct tsheader_proc));
         lseek(tsoutptr->fdw,0,SEEK_SET);
         off = lseek(tsoutptr->fdw,0,SEEK_END);

         npad = (bblen - off)/(off_t)(3*sizeof(float)*tsoutptr->np);
         if(npad > 0)
            {
            initval = (float *) check_malloc (3*tshead_p.nz*sizeof(float));
            for(i=0;i<3*tshead_p.nz;i++)
               initval[i] = 0;

            for(it=0;it<npad*tshead_p.localny;it++)
               rite(tsoutptr->fdw,initval,3*tshead_p.nz*sizeof(float));

            free(initval);
            }

         for(it=0;it<=rst->it;it=it+tsoutptr->idt)
            tsoutptr->nbufrite = tsoutptr->nbufrite + 1;

         bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*sizeof(float);
         tsoutptr->cur_off = sizeof(struct tsheader_proc) + bblen;
         }
      }

   lseek(tsoutptr->fdw,tsoutptr->cur_off,SEEK_SET);

   fprintf(stderr,"   local file= %s\n",tsoutptr->local_file);
   fprintf(stderr,"    main file= %s\n",tsoutptr->main_file);
   fprintf(stderr,"         nzts= %d\n",tshead_p.nz);
   fprintf(stderr,"         nyts= %d\n",tshead_p.localny);
   fprintf(stderr,"         ntts= %d\n",tshead_p.nt);
   fprintf(stderr,"         dtts= %d\n",outp->dtts);
   fprintf(stderr,"     yz-plane time slice output buffer is %.2f Mbytes\n",(3.0*(tsoutptr->np)*rpars->span*sizeof(float))/1.0e+06);
   fprintf(stderr,"     Output file is %.2f Mbytes\n\n",(3.0*tshead_p.nt*tshead_p.localny*tshead_p.nz*sizeof(float) + sizeof(struct tsheader_proc))/1.0e+06);
   fflush(stderr);

   if(di->enable_flag)
      strcpy(di->local_outputdir,outp->seisdir);
   }
else
   tsoutptr->iflag = 0;

sgtoutptr = &(outp->sgtpnts);
if(outp->sgtout)
   {
   fprintf(stderr,"**** Strain Green's tensor output:\n\n");
   fprintf(stderr,"     xmom=%13.5e\n",rpars->xmom);
   fprintf(stderr,"     ymom=%13.5e\n",rpars->ymom);
   fprintf(stderr,"     zmom=%13.5e\n\n",rpars->zmom);
   fflush(stderr);

   sgtmastptr = &(sgtoutptr->sgtmast);

   iyleft = rpars->ny1 + 2;
   if(rpars->ny1 == 0) iyleft = 0;

   iyright = rpars->ny2 - 3;
   if(rpars->ny2 == rpars->globny) iyright = rpars->ny2 - 1;

   fpr = fopfile(outp->sgtcords,"r");
   fprintf(stderr,"     Output coordinate locations file: %s\n",outp->sgtcords);

   sgtoutptr->iflag = 1;
   sgtoutptr->buflen = 0;
   sgtoutptr->nbufrite = 0;
   sgtoutptr->flushbuf = 0;

   sgtmastptr->geoproj = rpars->geoproj;
   sgtmastptr->modellon = rpars->modellon;
   sgtmastptr->modellat = rpars->modellat;
   sgtmastptr->modelrot = rpars->modelrot;
   sgtmastptr->xshift = rpars->xshift;
   sgtmastptr->yshift = rpars->yshift;
   sgtmastptr->nt = 1 + (((int)((rpars->nt-1)/rpars->span) + 1)*rpars->span - 1)/(sgtoutptr->tinc);

   fgets(string,1024,fpr);
   while(strncmp(string,"#",1) == 0)  /* skip comment lines */
      {
      fprintf(stderr,"     %s",string);
      fgets(string,1024,fpr);
      }
   fprintf(stderr,"\n");
   fflush(stderr);

   sscanf(string,"%d",&(sgtmastptr->globnp));

   sgtoutptr->sgtindx = (struct sgtindex *) check_malloc ((sgtmastptr->globnp)*sizeof(struct sgtindex));
   sgtoutptr->sgthead = (struct sgtheader *) check_malloc ((sgtmastptr->globnp)*sizeof(struct sgtheader));

   sgtindxptr = sgtoutptr->sgtindx;
   sgtheadptr = sgtoutptr->sgthead;
   iseis = 0;
   for(i=0;i<(sgtmastptr->globnp);i++)
      {
      fgets(string,1024,fpr);
      sscanf(string,"%d %d %d %Ld %f %f %f",&sgtindxptr[i].xsgt,&sgtindxptr[i].ysgt,&sgtindxptr[i].zsgt,&sgtindxptr[i].indx,&lonsgt,&latsgt,&depsgt);

      sgtindxptr[i].h = rpars->h;

      if(sgtindxptr[i].ysgt >= iyleft && sgtindxptr[i].ysgt <= iyright)
         {
         sgtheadptr[iseis].indx = sgtindxptr[i].indx;
         sgtheadptr[iseis].geoproj = rpars->geoproj;
         sgtheadptr[iseis].modellon = rpars->modellon;
         sgtheadptr[iseis].modellat = rpars->modellat;
         sgtheadptr[iseis].modelrot = rpars->modelrot;
         sgtheadptr[iseis].xshift = rpars->xshift;
         sgtheadptr[iseis].yshift = rpars->yshift;
         sgtheadptr[iseis].xazim = 90.0 + rpars->modelrot;

         sgtheadptr[iseis].nt = sgtmastptr->nt;
         sgtheadptr[iseis].dt = rpars->dt*sgtoutptr->tinc;
         sgtheadptr[iseis].h = rpars->h;

         sgtheadptr[iseis].xsrc = srcs->ix[0];
         sgtheadptr[iseis].ysrc = srcs->iy[0] + rpars->ny1; /* store as global index */
         sgtheadptr[iseis].zsrc = srcs->iz[0];

         sgtheadptr[iseis].xsgt = sgtindxptr[i].xsgt;
         sgtheadptr[iseis].ysgt = sgtindxptr[i].ysgt;
         sgtheadptr[iseis].zsgt = sgtindxptr[i].zsgt;

         sgtheadptr[iseis].xmom = rpars->xmom;
         sgtheadptr[iseis].ymom = rpars->ymom;
         sgtheadptr[iseis].zmom = rpars->zmom;

         sgtheadptr[iseis].tst = -rpars->tdelay;

	 /*
	 k = 0;
	 if(sgtheadptr[iseis].xsrc >= 0 && sgtheadptr[iseis].xsrc < rpars->nx)
	    k = sgtheadptr[iseis].xsrc + (sgtheadptr[iseis].ysrc - iyleft)*rpars->nx;

         sgtheadptr[iseis].src_lat = rpars->mc[k].lat;
         sgtheadptr[iseis].src_lon = rpars->mc[k].lon;
	 */

	 xx = sgtheadptr[iseis].xsrc*rpars->h;
	 yy = sgtheadptr[iseis].ysrc*rpars->h;

	 gcproj(&xx,&yy,&(sgtheadptr[iseis].src_lon),&(sgtheadptr[iseis].src_lat),&rpars->erad,&rpars->g0,&rpars->b0,rpars->amat,rpars->ainv,0);
         sgtheadptr[iseis].src_dep = (srcs->iz[0] - 1)*rpars->h;

	 /*
	 k = 0;
	 if(sgtheadptr[iseis].xsgt >= 0 && sgtheadptr[iseis].xsgt < rpars->nx)
	    k = sgtheadptr[iseis].xsgt + (sgtheadptr[iseis].ysgt - iyleft)*rpars->nx;

         sgtheadptr[iseis].sgt_lat = rpars->mc[k].lat;
         sgtheadptr[iseis].sgt_lon = rpars->mc[k].lon;
	 */

	 xx = sgtheadptr[iseis].xsgt*rpars->h;
	 yy = sgtheadptr[iseis].ysgt*rpars->h;

	 gcproj(&xx,&yy,&(sgtheadptr[iseis].sgt_lon),&(sgtheadptr[iseis].sgt_lat),&rpars->erad,&rpars->g0,&rpars->b0,rpars->amat,rpars->ainv,0);
         sgtheadptr[iseis].sgt_dep = depsgt;

	 xsgt = sgtheadptr[iseis].xsgt - sgtheadptr[iseis].xsrc;
	 ysgt = sgtheadptr[iseis].ysgt - sgtheadptr[iseis].ysrc;
	 zsgt = sgtheadptr[iseis].zsgt - sgtheadptr[iseis].zsrc;

	 sgtheadptr[iseis].cdist = (rpars->h)*sqrt(xsgt*xsgt + ysgt*ysgt + zsgt*zsgt);

         iseis++;
         }
      }
   fclose(fpr);

   sgtmastptr->localnp = iseis;

   if(sgtmastptr->localnp)
      {
      sgtoutptr->sgthead = (struct sgtheader *) check_realloc (sgtoutptr->sgthead,(sgtmastptr->localnp)*sizeof(struct sgtheader));

      sgtoutptr->sgt = (float *) check_malloc (6*(sgtmastptr->localnp)*rpars->span*sizeof(float));

      sprintf(sgtoutptr->local_file,"%s/%s_sgt-%.5d.e3d",outp->seisdir,name,rpars->segmentId);
      sprintf(sgtoutptr->main_file,"%s/%s_sgt-%.5d.e3d",di->main_dir,di->name,rpars->segmentId);

      if(rst->read_flag == 0)
         sgtoutptr->fdw = croptrfile(sgtoutptr->local_file);
      else
         {
         fcp_read_write(sgtoutptr->main_file,sgtoutptr->local_file);

         sgtoutptr->fdw = opfile(sgtoutptr->local_file);
         }

      fprintf(stderr,"local file= %s\n",sgtoutptr->local_file);
      fprintf(stderr," main file= %s\n",sgtoutptr->main_file);
      fprintf(stderr,"     npnts= %d\n",sgtmastptr->localnp);
      fprintf(stderr,"     ntout= %d\n",sgtmastptr->nt);
      fprintf(stderr,"     dtout= %f\n",rpars->dt*sgtoutptr->tinc);
      fprintf(stderr,"     Strain Green's Tensor buffer is %.2f Mbytes\n",(6.0*(sgtmastptr->localnp)*rpars->span*sizeof(float))/1.0e+06);
      fprintf(stderr,"     Output file is %.2f Mbytes\n\n",(sizeof(struct sgtmaster) + (sgtmastptr->globnp)*sizeof(struct sgtindex) + (sgtmastptr->localnp)*(sizeof(struct sgtheader) + 6.0*(sgtmastptr->nt)*sizeof(float)))/1.0e+06);
      fflush(stderr);

   /* fill up file with header info and dummy time histories  */

      if(rst->read_flag == 0)
         {
         rite(sgtoutptr->fdw,sgtmastptr,sizeof(struct sgtmaster));
         rite(sgtoutptr->fdw,sgtindxptr,(sgtmastptr->globnp)*sizeof(struct sgtindex));

         sgtheadptr = sgtoutptr->sgthead;
         for(i=0;i<(sgtmastptr->localnp);i++)
            rite(sgtoutptr->fdw,&sgtheadptr[i],sizeof(struct sgtheader));

         initval = (float *) check_malloc (6*(sgtmastptr->nt)*sizeof(float));
         for(i=0;i<6*(sgtmastptr->nt);i++)
            initval[i] = 0;

         for(i=0;i<(sgtmastptr->localnp);i++)
            rite(sgtoutptr->fdw,initval,6*(sgtmastptr->nt)*sizeof(float));

         free(initval);

         lseek(sgtoutptr->fdw,0,SEEK_SET);
         sgtoutptr->cur_off = 0;
	 }
      else
	 {
	 rst_sgtmast = (struct sgtmaster *)check_malloc(sizeof(struct sgtmaster));
         reed(sgtoutptr->fdw,rst_sgtmast,sizeof(struct sgtmaster));

	 if(rst_sgtmast->globnp != sgtmastptr->globnp)
	    {
	    fprintf(stderr,"***** rst_sgtmast->globnp (%d) != sgtmastptr->globnp (%d), file:\n",rst_sgtmast->globnp,sgtmastptr->globnp);
	    fprintf(stderr,"      %s,\t\texiting ...\n",sgtoutptr->local_file);
	    fflush(stderr);
	    exit(-1);
	    }

	 if(rst_sgtmast->localnp != sgtmastptr->localnp)
	    {
	    fprintf(stderr,"***** rst_sgtmast->localnp (%d) != sgtmastptr->localnp (%d), file:\n",rst_sgtmast->localnp,sgtmastptr->localnp);
	    fprintf(stderr,"      %s,\t\texiting ...\n",sgtoutptr->local_file);
	    fflush(stderr);
	    exit(-1);
	    }

	 free(rst_sgtmast);

/* if we made it here, everything is OK, re-write headers just in case nt has changed */

         lseek(sgtoutptr->fdw,0,SEEK_SET);

         rite(sgtoutptr->fdw,sgtmastptr,sizeof(struct sgtmaster));
         rite(sgtoutptr->fdw,sgtindxptr,(sgtmastptr->globnp)*sizeof(struct sgtindex));

         sgtheadptr = sgtoutptr->sgthead;
         for(i=0;i<(sgtmastptr->localnp);i++)
            rite(sgtoutptr->fdw,&sgtheadptr[i],sizeof(struct sgtheader));

/* check if we need to pad for new nt */

         bblen = sizeof(struct sgtmaster) + (sgtmastptr->globnp)*sizeof(struct sgtindex)
                + (sgtmastptr->localnp)*sizeof(struct sgtheader)
	        + 6*(sgtmastptr->localnp)*(sgtmastptr->nt)*sizeof(float);

         lseek(sgtoutptr->fdw,0,SEEK_SET);
         off = lseek(sgtoutptr->fdw,0,SEEK_END);

	 npad = (bblen - off)/(off_t)(3*sizeof(float)*(sgtmastptr->localnp));
	 if(npad > 0)
	    {
            initval = (float *) check_malloc (3*npad*sizeof(float));
            for(i=0;i<3*npad;i++)
               initval[i] = 0;

            for(i=0;i<(sgtmastptr->localnp);i++)
               rite(sgtoutptr->fdw,initval,3*npad*sizeof(float));

	    free(initval);
	    }

	 for(it=0;it<=rst->it;it=it+sgtoutptr->tinc)
	    sgtoutptr->nbufrite = sgtoutptr->nbufrite + 1;

         sgtoutptr->cur_off = sizeof(struct sgtmaster)
	        + (sgtmastptr->globnp)*sizeof(struct sgtindex)
                + (sgtmastptr->localnp)*sizeof(struct sgtheader)
	        + 6*(sgtmastptr->localnp)*(sgtoutptr->nbufrite)*sizeof(float);

         lseek(sgtoutptr->fdw,sgtoutptr->cur_off,SEEK_SET);
	 }
      }
   else
      {
      sgtmastptr->localnp = 0;
      sgtoutptr->iflag = 0;
      free(sgtoutptr->sgthead);
      }

   if(di->enable_flag)
      strcpy(di->local_outputdir,outp->seisdir);
   }
else
   sgtoutptr->iflag = 0;
}

void check_rpars(char *file,struct runparams *rp0,struct runparams *rp1)
{
int fail_flag;

fail_flag = 0;
if(rp0->segmentId != rp1->segmentId)
   fail_flag = 1;
if(rp0->nodeType != rp1->nodeType)
   fail_flag = 2;
if(rp0->nx != rp1->nx)
   fail_flag = 3;
if(rp0->globny != rp1->globny)
   fail_flag = 4;
if(rp0->nz != rp1->nz)
   fail_flag = 5;

/*   now we allow for different NT
if(rp0->nt != rp1->nt)
   fail_flag = 6;
*/

if(rp0->ny1 != rp1->ny1)
   fail_flag = 7;
if(rp0->ny2 != rp1->ny2)
   fail_flag = 8;
if(rp0->span != rp1->span)
   fail_flag = 9;
if(rp0->freesurf != rp1->freesurf)
   fail_flag = 10;
if(rp0->h != rp1->h)
   fail_flag = 11;
if(rp0->dt != rp1->dt)
   fail_flag = 12;
if(rp0->geoproj != rp1->geoproj)
   fail_flag = 13;
if(rp0->modelrot != rp1->modelrot)
   fail_flag = 14;
if(rp0->modellon != rp1->modellon)
   fail_flag = 15;
if(rp0->modellat != rp1->modellat)
   fail_flag = 16;
if(rp0->xmom != rp1->xmom)
   fail_flag = 17;
if(rp0->ymom != rp1->ymom)
   fail_flag = 18;
if(rp0->zmom != rp1->zmom)
   fail_flag = 19;
if(rp0->tdelay != rp1->tdelay)
   fail_flag = 20;

if(fail_flag != 0)
   {
   fprintf(stderr,"Problem with input parameters in restart file:\n");
   fprintf(stderr,"   %s\n",file);
   fprintf(stderr,"     fail_flag= %d\n\n",fail_flag);
   print_rpars(rp0,rp1);
   fprintf(stderr,"\n    exiting ...\n");
   fflush(stderr);
   mpi_exit(-1);
   }
}

void get_restart(struct restartinfo *rst,struct modelstorage *mds,float *pvf,int nx,int ny,int nz)
{
float *pvf_ptr;
int iy;

if(mds->intmem)
   {
   for(iy=0;iy<ny;iy++)
      {
      pvf_ptr  = pvf + iy*N_WAVE_VARS*nx*nz;
      reed(rst->fdr,pvf_ptr,N_WAVE_VARS*nx*nz*sizeof(float));
      }
   }
else
   {
   for(iy=0;iy<ny;iy++)
      {
      reed(rst->fdr,pvf,N_WAVE_VARS*nx*nz*sizeof(float));
      rite_model(mds,pvf,WAVE_FIELD);
      }
   init_model_seek(mds,WAVE_FIELD);
   }
close(rst->fdr);
}

void dump_restart(struct restartinfo *rst,char *name,int sid,struct runparams *rpars,struct modelstorage *mds,float *pvf,int it,int nx,int ny,int nz)
{
float *pvf_ptr;
int iy, rmflag;
char string[2048];
off_t file_len, reed_len, buf_len, nreed, default_buffer_length;
int fdr, fdw;
char *cbuf;

default_buffer_length = 10000000;   /* 10 Mbytes */

rst->it = it + rpars->span - 1;

if(strcmp(rst->dumpfile,"NOT_STARTED_YET") == 0)
   {
   sprintf(rst->dumpfile,"%s/%s_rst-%.5d.e3d",rst->dir,name,sid);
   rst->fdw = croptrfile(rst->dumpfile);

   sprintf(rst->bkupfile,"%s.backup",rst->dumpfile);
   rmflag = 0;
   }
else
   { /* RWG 06/28/07 NEW WAY: use reed() and rite() pairs to copy files */
   fcp_read_write(rst->dumpfile,rst->bkupfile);

   rmflag = 1;
   }

fprintf(stderr,"**** Dumping wavefield to restart file:\n");
fprintf(stderr,"     %s\n",rst->dumpfile);
fflush(stderr);

lseek(rst->fdw,0,SEEK_SET);
rite(rst->fdw,&rst->it,sizeof(int));
rite(rst->fdw,rpars,sizeof(struct runparams));

if(mds->intmem)
   {
   for(iy=0;iy<ny;iy++)
      {
      pvf_ptr  = pvf + iy*N_WAVE_VARS*nx*nz;
      rite(rst->fdw,pvf_ptr,N_WAVE_VARS*nx*nz*sizeof(float));
      }
   }
else
   {
   init_model_seek(mds,WAVE_FIELD);
   for(iy=0;iy<ny;iy++)
      {
      pvf  = reed_model(mds,pvf,WAVE_FIELD);
      rite(rst->fdw,pvf,N_WAVE_VARS*nx*nz*sizeof(float));
      }
   }

/* OLD WAY: use system call to explicitly remove file

if(rmflag != 0)
   {
   sprintf(string,"\\rm %s",rst->bkupfile);
   iy = system(string);
   }

sprintf(string,"ls -lt %s",rst->dir);
iy = system(string);
*/

/* RWG 06/28/07 NEW WAY: use unlink() to remove file */

if(rmflag != 0)
   unlink(rst->bkupfile);
}

void print_rpars(struct runparams *rp0,struct runparams *rp1)
{
fprintf(stderr,"rp0->segmentId= %d\trp1->segmentId= %d\n",rp0->segmentId,rp1->segmentId);
fprintf(stderr,"rp0->nodeType= %d\trp1->nodeType= %d\n",rp0->nodeType,rp1->nodeType);
fprintf(stderr,"rp0->nx= %d\trp1->nx= %d\n",rp0->nx,rp1->nx);
fprintf(stderr,"rp0->globny= %d\trp1->globny= %d\n",rp0->globny,rp1->globny);
fprintf(stderr,"rp0->nz= %d\trp1->nz= %d\n",rp0->nz,rp1->nz);
fprintf(stderr,"rp0->nt= %d\trp1->nt= %d\n",rp0->nt,rp1->nt);
fprintf(stderr,"rp0->ny1= %d\trp1->ny1= %d\n",rp0->ny1,rp1->ny1);
fprintf(stderr,"rp0->ny2= %d\trp1->ny2= %d\n",rp0->ny2,rp1->ny2);
fprintf(stderr,"rp0->span= %d\trp1->span= %d\n",rp0->span,rp1->span);
fprintf(stderr,"rp0->freesurf= %d\trp1->freesurf= %d\n",rp0->freesurf,rp1->freesurf);
fprintf(stderr,"rp0->h= %13.5e\trp1->h= %13.5e\n",rp0->h,rp1->h);
fprintf(stderr,"rp0->dt= %13.5e\trp1->dt= %13.5e\n",rp0->dt,rp1->dt);
fprintf(stderr,"rp0->modelrot= %13.5e\trp1->modelrot= %13.5e\n",rp0->modelrot,rp1->modelrot);
fprintf(stderr,"rp0->modellon= %13.5e\trp1->modellon= %13.5e\n",rp0->modellon,rp1->modellon);
fprintf(stderr,"rp0->modellat= %13.5e\trp1->modellat= %13.5e\n",rp0->modellat,rp1->modellat);
fprintf(stderr,"rp0->xmom= %13.5e\trp1->xmom= %13.5e\n",rp0->xmom,rp1->xmom);
fprintf(stderr,"rp0->ymom= %13.5e\trp1->ymom= %13.5e\n",rp0->ymom,rp1->ymom);
fprintf(stderr,"rp0->zmom= %13.5e\trp1->zmom= %13.5e\n",rp0->zmom,rp1->zmom);
fprintf(stderr,"rp0->tdelay= %13.5e\trp1->tdelay= %13.5e\n",rp0->tdelay,rp1->tdelay);
}

void dump_files(struct outputfields *outp,int sid,struct dump_output_info *di)
{
off_t file_len, reed_len, buf_len, nreed, default_buffer_length;
int fdr, fdw;
char *cbuf;
char string[2048], reed_file[2048], rite_file[2048];
struct tsoutputparams *tsoutptr;
struct seisoutputparams *soutptr;
struct sgtoutputparams *sgtoutptr;

default_buffer_length = 50000000;   /* 10 Mbytes */

fprintf(stderr,"**** Dumping all output files to:\n");
fprintf(stderr,"\t%s\n",di->main_dir);
fprintf(stderr,"...\t");
fflush(stderr);

/* RWG 06/28/07 NEW WAY: use reed() and rite() pairs to copy files */

soutptr = &(outp->spnts);
if(outp->nseis && soutptr->iflag)
   fcp_read_write(soutptr->local_file,soutptr->main_file);

tsoutptr = &(outp->xyslice);
if(outp->ts_xy && tsoutptr->iflag)
   fcp_read_write(tsoutptr->local_file,tsoutptr->main_file);

tsoutptr = &(outp->xzslice);
if(outp->ts_xz && tsoutptr->iflag)
   fcp_read_write(tsoutptr->local_file,tsoutptr->main_file);

tsoutptr = &(outp->yzslice);
if(outp->ts_yz && tsoutptr->iflag)
   fcp_read_write(tsoutptr->local_file,tsoutptr->main_file);

sgtoutptr = &(outp->sgtpnts);
if(outp->sgtout && sgtoutptr->iflag)
   fcp_read_write(sgtoutptr->local_file,sgtoutptr->main_file);

fprintf(stderr,"DONE\n\n");
fflush(stderr);
}

void rm_dump_files(struct outputfields *outp,int sid,struct dump_output_info *di)
{
char string[2048];
struct tsoutputparams *tsoutptr;
struct seisoutputparams *soutptr;
struct sgtoutputparams *sgtoutptr;

fprintf(stderr,"**** Removing local output files ... ");
fflush(stderr);

/* RWG 06/28/07 NEW WAY: use unlink() to remove files after closing file descriptor */

soutptr = &(outp->spnts);
if(outp->nseis && soutptr->iflag)
   {
   close(soutptr->fdw);
   unlink(soutptr->local_file);
   }

tsoutptr = &(outp->xyslice);
if(outp->ts_xy && tsoutptr->iflag)
   {
   close(tsoutptr->fdw);
   unlink(tsoutptr->local_file);
   }

tsoutptr = &(outp->xzslice);
if(outp->ts_xz && tsoutptr->iflag)
   {
   close(tsoutptr->fdw);
   unlink(tsoutptr->local_file);
   }

tsoutptr = &(outp->yzslice);
if(outp->ts_yz && tsoutptr->iflag)
   {
   close(tsoutptr->fdw);
   unlink(tsoutptr->local_file);
   }

sgtoutptr = &(outp->sgtpnts);
if(outp->sgtout && sgtoutptr->iflag)
   {
   close(sgtoutptr->fdw);
   unlink(sgtoutptr->local_file);
   }

fprintf(stderr,"DONE\n\n");
fflush(stderr);
}

int fcp_read_write(char *reed_file,char *rite_file)
{
off_t file_len, reed_len, buf_len, nreed, default_buffer_length;
int fdr, fdw;
char *cbuf;
struct stat sbuf;

fdr = stat(reed_file,&sbuf); /* stat reed_file to see if it actually exists */

if(fdr == -1 && errno == ENOENT) /* file doesn't exist */
   return(-1);

default_buffer_length = 25000000;   /* 25 Mbytes */

cbuf = (char *)check_malloc(default_buffer_length);

fdr = opfile_ro(reed_file);
fdw = cropfile_rw(rite_file);

file_len = lseek(fdr,(off_t)(0),SEEK_END);  /* brute force find file length */
lseek(fdr,(off_t)(0),SEEK_SET);

buf_len = default_buffer_length;
if(buf_len > file_len)
   buf_len = file_len;

reed_len = 0;
while(reed_len < file_len)
   {
   nreed = reed(fdr,cbuf,buf_len);
   reed_len = reed_len + nreed;

   rite(fdw,cbuf,buf_len);

   if((reed_len+buf_len) > file_len)
      buf_len = file_len - reed_len;
   }

close(fdr);
close(fdw);
free(cbuf);

return(0);
}

/*
eventually need to remove rpOLD from init_outpP3

options not working/checked yet:

        point seismos: programmed 20111116, still needs checking
        xy time slice: not yet
        xz time slice: not yet
        yz time slice: not yet
        sgt output: not yet
*/

void init_outpP3(struct outputfields *outp,struct runparamsP3 *rpars,char *name,struct pntsrcs *srcs,struct restartinfoP3 *rst,struct dump_output_info *di,struct runparams *rpOLD)
{
FILE *fpr;
struct tsoutputparams *tsoutptr;
struct tsheader_proc tshead_p;
struct tsheader tshead;
struct seisoutputparams *soutptr;
struct seisheader *sheadptr;
struct sgtoutputparams *sgtoutptr;
struct sgtmaster *sgtmastptr, *rst_sgtmast;
struct sgtindex *sgtindxptr;
struct sgtheader *sgtheadptr;
char string[2048], reed_file[2048];

int k, it, iy, yend, i, iseis, iyleft, iyright, npad;
int xsgt, ysgt, zsgt;
float lonsgt, latsgt, depsgt;
float xx, yy;
int npread;
off_t off, bblen;
float *initval;

struct nodeinfo *ni;
int nseis_comps = 3;

/* change 3 components to 9 components for nseis output */
nseis_comps = 9;

ni = &(rpars->ni);

makedir(outp->seisdir);

if(outp->nseis)
   {
   fprintf(stderr,"**** Individual station locations:\n");
   fflush(stderr);

   fpr = fopfile(outp->seiscords,"r");
   fscanf(fpr,"%d",&outp->nseis);

   soutptr = &(outp->spnts);

   soutptr->iflag = 1;  /* backward compatible with all_in_one=1 */

   if(soutptr->iflag == 1)
      {
      soutptr->buflen = 0;
      soutptr->nbufrite = 0;
      soutptr->flushbuf = 0;

      soutptr->ntout = 1 + (rpars->nt-1)/(soutptr->tinc);

      soutptr->shead = (struct seisheader *) check_malloc (outp->nseis*sizeof(struct seisheader));

      sheadptr = soutptr->shead;
      iseis = 0;
      for(i=0;i<outp->nseis;i++)
         {
         fscanf(fpr,"%d %d %d %s",&(sheadptr[iseis].ix),
                                  &(sheadptr[iseis].iy),
                                  &(sheadptr[iseis].iz),
                                  string);

	 if(sheadptr[iseis].ix >= ni->ixminus && sheadptr[iseis].ix <= ni->ixplus &&
	       sheadptr[iseis].iy >= ni->iyminus && sheadptr[iseis].iy <= ni->iyplus &&
	          sheadptr[iseis].iz >= ni->izminus && sheadptr[iseis].iz <= ni->izplus)
	    {
            strncpy(sheadptr[iseis].name,string,NBYTE_STATNAME-1);
            sheadptr[iseis].name[NBYTE_STATNAME-1] = '\0';

            sheadptr[iseis].indx = i;
            sheadptr[iseis].nt = soutptr->ntout;
            sheadptr[iseis].dt = rpars->dt*soutptr->tinc;
            sheadptr[iseis].h = rpars->h;
            sheadptr[iseis].modelrot = rpars->modelrot;

	    xx = sheadptr[iseis].ix*rpars->h;
	    yy = sheadptr[iseis].iy*rpars->h;

	    gcproj(&xx,&yy,&(sheadptr[iseis].slon),&(sheadptr[iseis].slat),&rpars->erad,&rpars->g0,&rpars->b0,rpars->amat,rpars->ainv,0);

            iseis++;
	    }
         }
      fclose(fpr);

      outp->nseis = iseis;

      if(outp->nseis)
         {
         soutptr->np = outp->nseis;

         soutptr->shead = (struct seisheader *) check_realloc (soutptr->shead,(soutptr->np)*sizeof( struct seisheader));
         soutptr->s = (float *) check_malloc (nseis_comps*(soutptr->np)*sizeof(float));

         sprintf(soutptr->local_file,"%s/%s_seis-%.5d.e3d",outp->seisdir,name,ni->segmentId);
         sprintf(soutptr->main_file,"%s/%s_seis-%.5d.e3d",di->main_dir,di->name,ni->segmentId);

	 if(rst->read_flag == 0)
            soutptr->fdw = croptrfile(soutptr->local_file);
	 else
            {
            fcp_read_write(soutptr->main_file,soutptr->local_file);
            soutptr->fdw = opfile(soutptr->local_file);
            }

         fprintf(stderr,"     ALL seismograms written into single output file\n\n");

         fprintf(stderr,"     local file= %s\n",soutptr->local_file);
         fprintf(stderr,"     main file= %s\n",soutptr->main_file);
         fprintf(stderr,"     nseis= %d\n",soutptr->np);
         fprintf(stderr,"     ntout= %d\n",soutptr->ntout);
         fprintf(stderr,"     dtout= %f\n",rpars->dt*soutptr->tinc);
         fprintf(stderr,"     seismogram output buffer is %.2f Mbytes\n",(1.0*nseis_comps*(soutptr->np)*sizeof(float))/1.0e+06);
         fprintf(stderr,"     Output file is %.2f Mbytes\n\n",(soutptr->np)*(sizeof(struct seisheader) + 1.0*nseis_comps*(soutptr->ntout)*sizeof(float))/1.0e+06);
         fflush(stderr);

         fprintf(stderr,"     stations in this node:");
         for(i=0;i<(soutptr->np);i++)
            fprintf(stderr,"\t%s",sheadptr[i].name);
         fprintf(stderr,"\n");
         fflush(stderr);

	 if(rst->read_flag == 0)
            rite(soutptr->fdw,&(soutptr->np),sizeof(int));
	 else
	    {
            reed(soutptr->fdw,&npread,sizeof(int));
	    if(npread != soutptr->np)
	       {
	       fprintf(stderr,"***** npread (%d) != soutptr->np (%d), file:\n",npread,soutptr->np);
	       fprintf(stderr,"      %s,\t\texiting ...\n",soutptr->local_file);
	       exit(-1);
	       }
	    }

         sheadptr = soutptr->shead;
         for(i=0;i<(soutptr->np);i++)
	    {
	    /* re-write headers just in case nt has changed */
            rite(soutptr->fdw,&sheadptr[i],sizeof(struct seisheader));

	    /* now adjust header indices into local system */
	    sheadptr[i].ix = sheadptr[i].ix - ni->nx1;
	    sheadptr[i].iy = sheadptr[i].iy - ni->ny1;
	    sheadptr[i].iz = sheadptr[i].iz - ni->nz1;
	    }

	 if(rst->read_flag == 1)
	    {
            bblen = sizeof(int) + (off_t)((soutptr->np)*(sizeof(struct seisheader) + nseis_comps*(soutptr->ntout)*sizeof(float)));
            lseek(soutptr->fdw,0,SEEK_SET);
            off = lseek(soutptr->fdw,0,SEEK_END);

	    npad = (bblen - off)/(off_t)(nseis_comps*sizeof(float)*(soutptr->np));
	    if(npad > 0)
	       {
               initval = (float *) check_malloc (nseis_comps*npad*sizeof(float));
               for(i=0;i<nseis_comps*npad;i++)
                  initval[i] = 0;

               for(i=0;i<(soutptr->np);i++)
                  rite(soutptr->fdw,initval,nseis_comps*npad*sizeof(float));

	       free(initval);
	       }

	    for(it=0;it<=rst->it;it=it+soutptr->tinc)
	       soutptr->nbufrite = soutptr->nbufrite + 1;

	    bblen = nseis_comps*(soutptr->np)*(soutptr->nbufrite)*sizeof(float);
            soutptr->cur_off = sizeof(int) + (soutptr->np)*sizeof(struct seisheader) + bblen;
	    }
	 else
	    {
            initval = (float *) check_malloc (nseis_comps*(soutptr->ntout)*sizeof(float));
            for(i=0;i<nseis_comps*(soutptr->ntout);i++)
               initval[i] = 0;

            for(i=0;i<(soutptr->np);i++)
               rite(soutptr->fdw,initval,nseis_comps*(soutptr->ntout)*sizeof(float));

	    free(initval);

            soutptr->cur_off = 0;
	    }

         lseek(soutptr->fdw,soutptr->cur_off,SEEK_SET);
	 }
      else
         {
         soutptr->np = outp->nseis;
         soutptr->iflag = 0;
         free(soutptr->shead);

         fprintf(stderr,"     NO seismograms written from this node\n");
         fflush(stderr);
	 }
      }
   fprintf(stderr,"\n");
   fflush(stderr);

   if(di->enable_flag && soutptr->iflag)
      strcpy(di->local_outputdir,outp->seisdir);
   }

/* rpOLD placeholder */

if(outp->ts_xy || outp->ts_xz || outp->ts_yz)
   {
   fprintf(stderr,"**** Time slices:\n");
   fflush(stderr);
   }

tsoutptr = &(outp->xyslice);
if(outp->ts_xy)
   {
   fprintf(stderr,"     xy-plane time slice:\n");
   fflush(stderr);

   tsoutptr->iflag = 1;

   iyleft = rpOLD->ny1 + 2;
   if(rpOLD->ny1 == 0) iyleft = 0;

   iyright = rpOLD->ny2 - 3;
   if(rpOLD->ny2 == rpOLD->globny) iyright = rpOLD->ny2 - 1;

   tshead_p.iyleft = -1;
   tshead_p.localny = 0;
   for(iy=0;iy<rpOLD->globny;iy++)
      {
      if(iy >= iyleft && iy <= iyright && iy%(outp->dyts) == 0)
         {
	 if(tshead_p.iyleft < 0)
	    tshead_p.iyleft = iy;

	 tshead_p.iyright = iy;

         tshead_p.localny = tshead_p.localny + 1;
	 }
      }

   tshead_p.ix0 = 0;
   tshead_p.iy0 = 0;
   tshead_p.iz0 = outp->iz_ts;
   tshead_p.it0 = 0;
   tshead_p.nx = (int)(((rpOLD->nx-1)/(outp->dxts))+1);
   tshead_p.ny = (int)(((rpOLD->globny-1)/(outp->dyts))+1);
   tshead_p.nz = 1;
   tshead_p.nt = (int)(((rpOLD->nt-1)/(outp->dtts))+1);
   tshead_p.dx = rpOLD->h*outp->dxts;
   tshead_p.dy = rpOLD->h*outp->dyts;
   tshead_p.dz = rpOLD->h;
   tshead_p.dt = rpOLD->dt*outp->dtts;
   tshead_p.modelrot = rpOLD->modelrot;
   tshead_p.modellat = rpOLD->modellat;
   tshead_p.modellon = rpOLD->modellon;

   tsoutptr->ix0 = 0;
   tsoutptr->iy0 = 0;
   tsoutptr->iz0 = outp->iz_ts;
   tsoutptr->it0 = 0;
   tsoutptr->idx = outp->dxts;
   tsoutptr->idy = outp->dyts;
   tsoutptr->idz = outp->iz_ts;
   tsoutptr->idt = outp->dtts;
   tsoutptr->nx = (int)(((rpOLD->nx-1)/(tsoutptr->idx))+1);
   tsoutptr->ny = tshead_p.localny;
   tsoutptr->np = (tsoutptr->nx)*(tsoutptr->ny);
   tsoutptr->ntout = tshead_p.nt;

   tsoutptr->iyleft = iyleft;
   tsoutptr->iyright = iyright;

   tsoutptr->buflen = 0;
   tsoutptr->nbufrite = 0;
   tsoutptr->flushbuf = 0;
   tsoutptr->head_off = sizeof(struct tsheader_proc);

   tsoutptr->vbuf = (float *) check_malloc (3*(tsoutptr->np)*(rpOLD->span)*sizeof(float));

   sprintf(tsoutptr->local_file,"%s/%s_xyts-%.5d.e3d",outp->seisdir,name,rpOLD->segmentId);
   sprintf(tsoutptr->main_file,"%s/%s_xyts-%.5d.e3d",di->main_dir,di->name,rpOLD->segmentId);

   if(rst->read_flag == 0)
      {
      tsoutptr->fdw = croptrfile(tsoutptr->local_file);

      rite(tsoutptr->fdw,&tshead_p,sizeof(struct tsheader_proc));

      initval = (float *) check_malloc (3*tshead_p.nx*sizeof(float));
      for(iy=0;iy<3*tshead_p.nx;iy++)
         initval[iy] = 0;

      for(it=0;it<tshead_p.nt*tshead_p.localny;it++)
         rite(tsoutptr->fdw,initval,3*tshead_p.nx*sizeof(float));

      free(initval);

      tsoutptr->cur_off = sizeof(struct tsheader_proc);
      }
   else
      {
      fcp_read_write(tsoutptr->main_file,tsoutptr->local_file);

      tsoutptr->fdw = opfile(tsoutptr->local_file);

      /* re-write header just in case nt has changed */
      rite(tsoutptr->fdw,&tshead_p,sizeof(struct tsheader_proc));

      if(tshead_p.localny > 0)
         {
         bblen = (off_t)(3*tsoutptr->ntout*tsoutptr->np*sizeof(float)) + (off_t)(sizeof(struct tsheader_proc));
         lseek(tsoutptr->fdw,0,SEEK_SET);
         off = lseek(tsoutptr->fdw,0,SEEK_END);

         npad = (bblen - off)/(off_t)(3*sizeof(float)*tsoutptr->np);
         if(npad > 0)
            {
            initval = (float *) check_malloc (3*tshead_p.nx*sizeof(float));
            for(i=0;i<3*tshead_p.nx;i++)
               initval[i] = 0;

            for(it=0;it<npad*tshead_p.localny;it++)
               rite(tsoutptr->fdw,initval,3*tshead_p.nx*sizeof(float));

            free(initval);
            }

         for(it=0;it<=rst->it;it=it+tsoutptr->idt)
            tsoutptr->nbufrite = tsoutptr->nbufrite + 1;

         bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*sizeof(float);
         tsoutptr->cur_off = sizeof(struct tsheader_proc) + bblen;
	 }
      }

   lseek(tsoutptr->fdw,tsoutptr->cur_off,SEEK_SET);

   fprintf(stderr,"   local file= %s\n",tsoutptr->local_file);
   fprintf(stderr,"    main file= %s\n",tsoutptr->main_file);
   fprintf(stderr,"         nxts= %d\n",tshead_p.nx);
   fprintf(stderr,"         nyts= %d\n",tshead_p.localny);
   fprintf(stderr,"         ntts= %d\n",tshead_p.nt);
   fprintf(stderr,"         dtts= %d\n",outp->dtts);
   fprintf(stderr,"     xy-plane time slice output buffer is %.2f Mbytes\n",(3.0*(tsoutptr->np)*rpOLD->span*sizeof(float))/1.0e+06);
   fprintf(stderr,"     Output file is %.2f Mbytes\n\n",(3.0*tshead_p.nt*tshead_p.localny*tshead_p.nx*sizeof(float) + sizeof(struct tsheader_proc))/1.0e+06);
   fflush(stderr);

   if(di->enable_flag)
      strcpy(di->local_outputdir,outp->seisdir);
   }
else
   tsoutptr->iflag = 0;

tsoutptr = &(outp->xzslice);
if(outp->ts_xz)
   {
   iyleft = rpOLD->ny1 + 2;
   if(rpOLD->ny1 == 0) iyleft = 0;

   iyright = rpOLD->ny2 - 3;
   if(rpOLD->ny2 == rpOLD->globny) iyright = rpOLD->ny2 - 1;

   tsoutptr->iflag = 0;
   if(outp->iy_ts >= iyleft && outp->iy_ts <= iyright)
      {
      fprintf(stderr,"     xz-plane time slice:\n");
      fflush(stderr);

      tsoutptr->iflag = 1;

      tshead.ix0 = 0;
      tshead.iy0 = outp->iy_ts;
      tshead.iz0 = 0;
      tshead.it0 = 0;
      tshead.nx = (int)(((rpOLD->nx-1)/(outp->dxts))+1);
      tshead.ny = 1;
      tshead.nz = (int)(((rpOLD->nz-1)/(outp->dzts))+1);
      tshead.nt = (int)(((rpOLD->nt-1)/(outp->dtts))+1);
      tshead.dx = rpOLD->h*outp->dxts;
      tshead.dy = rpOLD->h;
      tshead.dz = rpOLD->h*outp->dxts;
      tshead.dt = rpOLD->dt*outp->dtts;
      tshead.modelrot = rpOLD->modelrot;
      tshead.modellat = rpOLD->modellat;
      tshead.modellon = rpOLD->modellon;

      tsoutptr->ix0 = 0;
      tsoutptr->iy0 = outp->iy_ts;
      tsoutptr->iz0 = 0;
      if(rpOLD->freesurf == 1)
         tsoutptr->iz0 = 1;

      tsoutptr->it0 = 0;
      tsoutptr->idx = outp->dxts;
      tsoutptr->idy = outp->iy_ts;
      tsoutptr->idz = outp->dzts;
      tsoutptr->idt = outp->dtts;
      tsoutptr->nx = (int)(((rpOLD->nx-1)/(tsoutptr->idx))+1);
      tsoutptr->ny = 1;
      tsoutptr->nz = (int)(((rpOLD->nz-1)/(tsoutptr->idz))+1);
      tsoutptr->np = (tsoutptr->nx)*(tsoutptr->nz);
      tsoutptr->ntout = tshead.nt;

      tsoutptr->iyleft = iyleft;
      tsoutptr->iyright = iyright;

      tsoutptr->buflen = 0;
      tsoutptr->nbufrite = 0;
      tsoutptr->flushbuf = 0;
      tsoutptr->head_off = sizeof(struct tsheader);

      tsoutptr->vbuf = (float *) check_malloc (3*(tsoutptr->np)*(rpOLD->span)*sizeof(float));

      sprintf(tsoutptr->local_file,"%s/%s_xzts.e3d",outp->seisdir,name);
      sprintf(tsoutptr->main_file,"%s/%s_xzts.e3d",di->main_dir,di->name);

      if(rst->read_flag == 0)
         {
         tsoutptr->fdw = croptrfile(tsoutptr->local_file);

         rite(tsoutptr->fdw,&tshead,sizeof(struct tsheader));

         initval = (float *) check_malloc (3*tshead.nx*sizeof(float));
         for(iy=0;iy<3*tshead.nx;iy++)
            initval[iy] = 0;

         for(it=0;it<tshead.nt*tshead.nz;it++)
            rite(tsoutptr->fdw,initval,3*tshead.nx*sizeof(float));

         free(initval);

         tsoutptr->cur_off = sizeof(struct tsheader);
         }
      else
         {
         fcp_read_write(tsoutptr->main_file,tsoutptr->local_file);

         tsoutptr->fdw = opfile(tsoutptr->local_file);

      /* re-write header just in case nt has changed */
         rite(tsoutptr->fdw,&tshead,sizeof(struct tsheader));

         bblen = (off_t)(3*tsoutptr->ntout*tsoutptr->np*sizeof(float)) + (off_t)(sizeof(struct tsheader));
         lseek(tsoutptr->fdw,0,SEEK_SET);
         off = lseek(tsoutptr->fdw,0,SEEK_END);

         npad = (bblen - off)/(off_t)(3*sizeof(float)*tsoutptr->np);
         if(npad > 0)
            {
            initval = (float *) check_malloc (3*tshead.nx*sizeof(float));
            for(i=0;i<3*tshead.nx;i++)
               initval[i] = 0;

            for(it=0;it<npad*tshead.nz;it++)
               rite(tsoutptr->fdw,initval,3*tshead.nx*sizeof(float));

            free(initval);
            }

         for(it=0;it<=rst->it;it=it+tsoutptr->idt)
            tsoutptr->nbufrite = tsoutptr->nbufrite + 1;

         bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*sizeof(float);
         tsoutptr->cur_off = sizeof(struct tsheader) + bblen;
         }

      lseek(tsoutptr->fdw,tsoutptr->cur_off,SEEK_SET);

      fprintf(stderr,"   local file= %s\n",tsoutptr->local_file);
      fprintf(stderr,"    main file= %s\n",tsoutptr->main_file);
      fprintf(stderr,"         nxts= %d\n",tshead.nx);
      fprintf(stderr,"         nzts= %d\n",tshead.nz);
      fprintf(stderr,"         ntts= %d\n",tshead.nt);
      fprintf(stderr,"         dtts= %d\n",outp->dtts);
      fprintf(stderr,"     xz-plane time slice output buffer is %.2f Mbytes\n",(3.0*(tsoutptr->np)*rpOLD->span*sizeof(float))/1.0e+06);
      fprintf(stderr,"     Output file is %.2f Mbytes\n\n",(3.0*tshead.nt*tshead.nx*sizeof(float) + sizeof(struct tsheader))/1.0e+06);
      fflush(stderr);

      if(di->enable_flag)
         strcpy(di->local_outputdir,outp->seisdir);
      }
   }
else
   tsoutptr->iflag = 0;

tsoutptr = &(outp->yzslice);
if(outp->ts_yz)
   {
   fprintf(stderr,"     yz-plane time slice:\n");
   fflush(stderr);

   tsoutptr->iflag = 1;

   iyleft = rpOLD->ny1 + 2;
   if(rpOLD->ny1 == 0) iyleft = 0;

   iyright = rpOLD->ny2 - 3;
   if(rpOLD->ny2 == rpOLD->globny) iyright = rpOLD->ny2 - 1;

   tshead_p.iyleft = -1;
   tshead_p.localny = 0;
   for(iy=0;iy<rpOLD->globny;iy++)
      {
      if(iy >= iyleft && iy <= iyright && iy%(outp->dyts) == 0)
         {
         if(tshead_p.iyleft < 0)
            tshead_p.iyleft = iy;

         tshead_p.iyright = iy;

         tshead_p.localny = tshead_p.localny + 1;
         }
      }

   tshead_p.ix0 = outp->ix_ts;
   tshead_p.iy0 = 0;
   tshead_p.iz0 = 0;
   tshead_p.it0 = 0;
   tshead_p.nx = 1;
   tshead_p.ny = (int)(((rpOLD->globny-1)/(outp->dyts))+1);
   tshead_p.nz = (int)(((rpOLD->nz-1)/(outp->dzts))+1);
   tshead_p.nx = tshead_p.nz;  /* to make merge_ts work 2011-11-09 */
   tshead_p.nt = (int)(((rpOLD->nt-1)/(outp->dtts))+1);
   tshead_p.dx = rpOLD->h;
   tshead_p.dy = rpOLD->h*outp->dyts;
   tshead_p.dz = rpOLD->h*outp->dxts;
   tshead_p.dt = rpOLD->dt*outp->dtts;
   tshead_p.modelrot = rpOLD->modelrot;
   tshead_p.modellat = rpOLD->modellat;
   tshead_p.modellon = rpOLD->modellon;

   tsoutptr->ix0 = outp->ix_ts;
   tsoutptr->iy0 = 0;
   tsoutptr->iz0 = 0;
   tsoutptr->it0 = 0;
   tsoutptr->idx = outp->ix_ts;
   tsoutptr->idy = outp->dyts;
   tsoutptr->idz = outp->dzts;
   tsoutptr->idt = outp->dtts;
   tsoutptr->nz = (int)(((rpOLD->nz-1)/(tsoutptr->idz))+1);
   tsoutptr->ny = tshead_p.localny;
   tsoutptr->np = (tsoutptr->nz)*(tsoutptr->ny);
   tsoutptr->ntout = tshead_p.nt;

   tsoutptr->iyleft = iyleft;
   tsoutptr->iyright = iyright;

   tsoutptr->buflen = 0;
   tsoutptr->nbufrite = 0;
   tsoutptr->flushbuf = 0;
   tsoutptr->head_off = sizeof(struct tsheader_proc);

   tsoutptr->vbuf = (float *) check_malloc (3*(tsoutptr->np)*(rpOLD->span)*sizeof(float));

   sprintf(tsoutptr->local_file,"%s/%s_yzts-%.4d.e3d",outp->seisdir,name,rpOLD->segmentId);
   sprintf(tsoutptr->main_file,"%s/%s_yzts-%.4d.e3d",di->main_dir,di->name,rpOLD->segmentId);

   if(rst->read_flag == 0)
      {
      tsoutptr->fdw = croptrfile(tsoutptr->local_file);

      rite(tsoutptr->fdw,&tshead_p,sizeof(struct tsheader_proc));

      initval = (float *) check_malloc (3*tshead_p.nz*sizeof(float));
      for(iy=0;iy<3*tshead_p.nz;iy++)
         initval[iy] = 0;

      for(it=0;it<tshead_p.nt*tshead_p.localny;it++)
         rite(tsoutptr->fdw,initval,3*tshead_p.nz*sizeof(float));

      free(initval);

      tsoutptr->cur_off = sizeof(struct tsheader_proc);
      }
   else
      {
      fcp_read_write(tsoutptr->main_file,tsoutptr->local_file);

      tsoutptr->fdw = opfile(tsoutptr->local_file);

      /* re-write header just in case nt has changed */
      rite(tsoutptr->fdw,&tshead_p,sizeof(struct tsheader_proc));

      if(tshead_p.localny > 0)
         {
         bblen = (off_t)(3*tsoutptr->ntout*tsoutptr->np*sizeof(float)) + (off_t)(sizeof(struct tsheader_proc));
         lseek(tsoutptr->fdw,0,SEEK_SET);
         off = lseek(tsoutptr->fdw,0,SEEK_END);

         npad = (bblen - off)/(off_t)(3*sizeof(float)*tsoutptr->np);
         if(npad > 0)
            {
            initval = (float *) check_malloc (3*tshead_p.nz*sizeof(float));
            for(i=0;i<3*tshead_p.nz;i++)
               initval[i] = 0;

            for(it=0;it<npad*tshead_p.localny;it++)
               rite(tsoutptr->fdw,initval,3*tshead_p.nz*sizeof(float));

            free(initval);
            }

         for(it=0;it<=rst->it;it=it+tsoutptr->idt)
            tsoutptr->nbufrite = tsoutptr->nbufrite + 1;

         bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*sizeof(float);
         tsoutptr->cur_off = sizeof(struct tsheader_proc) + bblen;
         }
      }

   lseek(tsoutptr->fdw,tsoutptr->cur_off,SEEK_SET);

   fprintf(stderr,"   local file= %s\n",tsoutptr->local_file);
   fprintf(stderr,"    main file= %s\n",tsoutptr->main_file);
   fprintf(stderr,"         nzts= %d\n",tshead_p.nz);
   fprintf(stderr,"         nyts= %d\n",tshead_p.localny);
   fprintf(stderr,"         ntts= %d\n",tshead_p.nt);
   fprintf(stderr,"         dtts= %d\n",outp->dtts);
   fprintf(stderr,"     yz-plane time slice output buffer is %.2f Mbytes\n",(3.0*(tsoutptr->np)*rpOLD->span*sizeof(float))/1.0e+06);
   fprintf(stderr,"     Output file is %.2f Mbytes\n\n",(3.0*tshead_p.nt*tshead_p.localny*tshead_p.nz*sizeof(float) + sizeof(struct tsheader_proc))/1.0e+06);
   fflush(stderr);

   if(di->enable_flag)
      strcpy(di->local_outputdir,outp->seisdir);
   }
else
   tsoutptr->iflag = 0;

sgtoutptr = &(outp->sgtpnts);
if(outp->sgtout)
   {
   fprintf(stderr,"**** Strain Green's tensor output:\n\n");
   fprintf(stderr,"     xmom=%13.5e\n",rpars->xmom);
   fprintf(stderr,"     ymom=%13.5e\n",rpars->ymom);
   fprintf(stderr,"     zmom=%13.5e\n\n",rpars->zmom);
   fflush(stderr);

   sgtmastptr = &(sgtoutptr->sgtmast);

   fpr = fopfile(outp->sgtcords,"r");
   fprintf(stderr,"     Output coordinate locations file: %s\n",outp->sgtcords);

   sgtoutptr->iflag = 1;
   sgtoutptr->buflen = 0;
   sgtoutptr->nbufrite = 0;
   sgtoutptr->flushbuf = 0;

   sgtmastptr->geoproj = rpars->geoproj;
   sgtmastptr->modellon = rpars->modellon;
   sgtmastptr->modellat = rpars->modellat;
   sgtmastptr->modelrot = rpars->modelrot;
   sgtmastptr->xshift = rpars->xshift;
   sgtmastptr->yshift = rpars->yshift;
   sgtmastptr->nt = 1 + (rpars->nt-1)/(sgtoutptr->tinc);

   fgets(string,1024,fpr);
   while(strncmp(string,"#",1) == 0)  /* skip comment lines */
      {
      fprintf(stderr,"     %s",string);
      fgets(string,1024,fpr);
      }
   fprintf(stderr,"\n");
   fflush(stderr);

   sscanf(string,"%d",&(sgtmastptr->globnp));

   sgtoutptr->sgtindx = (struct sgtindex *) check_malloc ((sgtmastptr->globnp)*sizeof(struct sgtindex));
   sgtoutptr->sgthead = (struct sgtheader *) check_malloc ((sgtmastptr->globnp)*sizeof(struct sgtheader));

   sgtindxptr = sgtoutptr->sgtindx;
   sgtheadptr = sgtoutptr->sgthead;
   iseis = 0;
   for(i=0;i<(sgtmastptr->globnp);i++)
      {
      fgets(string,1024,fpr);
      sscanf(string,"%d %d %d %Ld %f %f %f",&sgtindxptr[i].xsgt,&sgtindxptr[i].ysgt,&sgtindxptr[i].zsgt,&sgtindxptr[i].indx,&lonsgt,&latsgt,&depsgt);

      sgtindxptr[i].h = rpars->h;

      if(sgtindxptr[i].xsgt >= ni->ixminus && sgtindxptr[i].xsgt <= ni->ixplus &&
            sgtindxptr[i].ysgt >= ni->iyminus && sgtindxptr[i].ysgt <= ni->iyplus &&
	       sgtindxptr[i].zsgt >= ni->izminus && sgtindxptr[i].zsgt <= ni->izplus)
         {
         sgtheadptr[iseis].indx = sgtindxptr[i].indx;
         sgtheadptr[iseis].geoproj = rpars->geoproj;
         sgtheadptr[iseis].modellon = rpars->modellon;
         sgtheadptr[iseis].modellat = rpars->modellat;
         sgtheadptr[iseis].modelrot = rpars->modelrot;
         sgtheadptr[iseis].xshift = rpars->xshift;
         sgtheadptr[iseis].yshift = rpars->yshift;
         sgtheadptr[iseis].xazim = 90.0 + rpars->modelrot;

         sgtheadptr[iseis].nt = sgtmastptr->nt;
         sgtheadptr[iseis].dt = rpars->dt*sgtoutptr->tinc;
         sgtheadptr[iseis].h = rpars->h;

	 /* store source and SGT locations as global indices */
	 if(srcs->ix[0] < 0 || srcs->iy[0] < 0 || srcs->iz[0] < 0) /* source not in this node */
	    {
            sgtheadptr[iseis].xsrc = -srcs->ix[0];
            sgtheadptr[iseis].ysrc = -srcs->iy[0];
            sgtheadptr[iseis].zsrc = -srcs->iz[0];
	    }
	 else
	    {
            sgtheadptr[iseis].xsrc = srcs->ix[0] + ni->nx1;;
            sgtheadptr[iseis].ysrc = srcs->iy[0] + ni->ny1;
            sgtheadptr[iseis].zsrc = srcs->iz[0] + ni->nz1;;
	    }

         sgtheadptr[iseis].xsgt = sgtindxptr[i].xsgt;
         sgtheadptr[iseis].ysgt = sgtindxptr[i].ysgt;
         sgtheadptr[iseis].zsgt = sgtindxptr[i].zsgt;

         sgtheadptr[iseis].xmom = rpars->xmom;
         sgtheadptr[iseis].ymom = rpars->ymom;
         sgtheadptr[iseis].zmom = rpars->zmom;

         sgtheadptr[iseis].tst = -rpars->tdelay;

	 xx = sgtheadptr[iseis].xsrc*rpars->h;
	 yy = sgtheadptr[iseis].ysrc*rpars->h;

	 gcproj(&xx,&yy,&(sgtheadptr[iseis].src_lon),&(sgtheadptr[iseis].src_lat),&rpars->erad,&rpars->g0,&rpars->b0,rpars->amat,rpars->ainv,0);

         sgtheadptr[iseis].src_dep = sgtheadptr[iseis].zsrc*rpars->h;
	 if(rpars->freesurf)
            sgtheadptr[iseis].src_dep = sgtheadptr[iseis].src_dep - rpars->h;

	 xx = sgtheadptr[iseis].xsgt*rpars->h;
	 yy = sgtheadptr[iseis].ysgt*rpars->h;

	 gcproj(&xx,&yy,&(sgtheadptr[iseis].sgt_lon),&(sgtheadptr[iseis].sgt_lat),&rpars->erad,&rpars->g0,&rpars->b0,rpars->amat,rpars->ainv,0);

         sgtheadptr[iseis].sgt_dep = depsgt;

	 xsgt = sgtheadptr[iseis].xsgt - sgtheadptr[iseis].xsrc;
	 ysgt = sgtheadptr[iseis].ysgt - sgtheadptr[iseis].ysrc;
	 zsgt = sgtheadptr[iseis].zsgt - sgtheadptr[iseis].zsrc;

	 sgtheadptr[iseis].cdist = (rpars->h)*sqrt(xsgt*xsgt + ysgt*ysgt + zsgt*zsgt);

         iseis++;
         }
      }
   fclose(fpr);

   sgtmastptr->localnp = iseis;

   if(sgtmastptr->localnp)
      {
      sgtoutptr->sgthead = (struct sgtheader *) check_realloc (sgtoutptr->sgthead,(sgtmastptr->localnp)*sizeof(struct sgtheader));

      sgtoutptr->sgt = (float *) check_malloc (6*(sgtmastptr->localnp)*sizeof(float));

      sprintf(sgtoutptr->local_file,"%s/%s_sgt-%.5d.e3d",outp->seisdir,name,ni->segmentId);
      sprintf(sgtoutptr->main_file,"%s/%s_sgt-%.5d.e3d",di->main_dir,di->name,ni->segmentId);

      if(rst->read_flag == 0)
         sgtoutptr->fdw = croptrfile(sgtoutptr->local_file);
      else
         {
         fcp_read_write(sgtoutptr->main_file,sgtoutptr->local_file);

         sgtoutptr->fdw = opfile(sgtoutptr->local_file);
         }

      fprintf(stderr,"local file= %s\n",sgtoutptr->local_file);
      fprintf(stderr," main file= %s\n",sgtoutptr->main_file);
      fprintf(stderr,"     npnts= %d\n",sgtmastptr->localnp);
      fprintf(stderr,"     ntout= %d\n",sgtmastptr->nt);
      fprintf(stderr,"     dtout= %f\n",rpars->dt*sgtoutptr->tinc);
      fprintf(stderr,"     Strain Green's Tensor buffer is %.2f Mbytes\n",(6.0*(sgtmastptr->localnp)*sizeof(float))/1.0e+06);
      fprintf(stderr,"     Output file is %.2f Mbytes\n\n",(sizeof(struct sgtmaster) + (sgtmastptr->globnp)*sizeof(struct sgtindex) + (sgtmastptr->localnp)*(sizeof(struct sgtheader) + 6.0*(sgtmastptr->nt)*sizeof(float)))/1.0e+06);
      fflush(stderr);

   /* fill up file with header info and dummy time histories  */

      if(rst->read_flag == 0)
         {
         rite(sgtoutptr->fdw,sgtmastptr,sizeof(struct sgtmaster));
         rite(sgtoutptr->fdw,sgtindxptr,(sgtmastptr->globnp)*sizeof(struct sgtindex));

         sgtheadptr = sgtoutptr->sgthead;
         for(i=0;i<(sgtmastptr->localnp);i++)
            rite(sgtoutptr->fdw,&sgtheadptr[i],sizeof(struct sgtheader));

         initval = (float *) check_malloc (6*(sgtmastptr->nt)*sizeof(float));
         for(i=0;i<6*(sgtmastptr->nt);i++)
            initval[i] = 0;

         for(i=0;i<(sgtmastptr->localnp);i++)
            rite(sgtoutptr->fdw,initval,6*(sgtmastptr->nt)*sizeof(float));

         free(initval);

         lseek(sgtoutptr->fdw,0,SEEK_SET);
         sgtoutptr->cur_off = 0;
	 }
      else
	 {
	 rst_sgtmast = (struct sgtmaster *)check_malloc(sizeof(struct sgtmaster));
         reed(sgtoutptr->fdw,rst_sgtmast,sizeof(struct sgtmaster));

	 if(rst_sgtmast->globnp != sgtmastptr->globnp)
	    {
	    fprintf(stderr,"***** rst_sgtmast->globnp (%d) != sgtmastptr->globnp (%d), file:\n",rst_sgtmast->globnp,sgtmastptr->globnp);
	    fprintf(stderr,"      %s,\t\texiting ...\n",sgtoutptr->local_file);
	    fflush(stderr);
	    exit(-1);
	    }

	 if(rst_sgtmast->localnp != sgtmastptr->localnp)
	    {
	    fprintf(stderr,"***** rst_sgtmast->localnp (%d) != sgtmastptr->localnp (%d), file:\n",rst_sgtmast->localnp,sgtmastptr->localnp);
	    fprintf(stderr,"      %s,\t\texiting ...\n",sgtoutptr->local_file);
	    fflush(stderr);
	    exit(-1);
	    }

	 free(rst_sgtmast);

/* if we made it here, everything is OK, re-write headers just in case nt has changed */

         lseek(sgtoutptr->fdw,0,SEEK_SET);

         rite(sgtoutptr->fdw,sgtmastptr,sizeof(struct sgtmaster));
         rite(sgtoutptr->fdw,sgtindxptr,(sgtmastptr->globnp)*sizeof(struct sgtindex));

         sgtheadptr = sgtoutptr->sgthead;
         for(i=0;i<(sgtmastptr->localnp);i++)
            rite(sgtoutptr->fdw,&sgtheadptr[i],sizeof(struct sgtheader));

/* check if we need to pad for new nt */

         bblen = sizeof(struct sgtmaster) + (sgtmastptr->globnp)*sizeof(struct sgtindex)
                + (sgtmastptr->localnp)*sizeof(struct sgtheader)
	        + 6*(sgtmastptr->localnp)*(sgtmastptr->nt)*sizeof(float);

         lseek(sgtoutptr->fdw,0,SEEK_SET);
         off = lseek(sgtoutptr->fdw,0,SEEK_END);

	 npad = (bblen - off)/(off_t)(3*sizeof(float)*(sgtmastptr->localnp));
	 if(npad > 0)
	    {
            initval = (float *) check_malloc (3*npad*sizeof(float));
            for(i=0;i<3*npad;i++)
               initval[i] = 0;

            for(i=0;i<(sgtmastptr->localnp);i++)
               rite(sgtoutptr->fdw,initval,3*npad*sizeof(float));

	    free(initval);
	    }

	 for(it=0;it<=rst->it;it=it+sgtoutptr->tinc)
	    sgtoutptr->nbufrite = sgtoutptr->nbufrite + 1;

         sgtoutptr->cur_off = sizeof(struct sgtmaster)
	        + (sgtmastptr->globnp)*sizeof(struct sgtindex)
                + (sgtmastptr->localnp)*sizeof(struct sgtheader)
	        + 6*(sgtmastptr->localnp)*(sgtoutptr->nbufrite)*sizeof(float);

         lseek(sgtoutptr->fdw,sgtoutptr->cur_off,SEEK_SET);
	 }
      }
   else
      {
      sgtmastptr->localnp = 0;
      sgtoutptr->iflag = 0;
      free(sgtoutptr->sgthead);
      }

   if(di->enable_flag)
      strcpy(di->local_outputdir,outp->seisdir);
   }
else
   sgtoutptr->iflag = 0;
}

void check_rparsP3(char *file,struct runparamsP3 *rp0,struct runparamsP3 *rp1)
{
int fail_flag;
struct nodeinfo *ni0, *ni1;

ni0 = &(rp0->ni);
ni1 = &(rp1->ni);

fail_flag = 0;

if(ni0->nproc != ni1->nproc)
   fail_flag++;
if(ni0->nproc_x != ni1->nproc_x)
   fail_flag++;
if(ni0->nproc_y != ni1->nproc_y)
   fail_flag++;
if(ni0->nproc_z != ni1->nproc_z)
   fail_flag++;
if(ni0->min_nproc != ni1->min_nproc)
   fail_flag++;
if(ni0->segmentId != ni1->segmentId)
   fail_flag++;
if(ni0->minusId_x != ni1->minusId_x)
   fail_flag++;
if(ni0->plusId_x != ni1->plusId_x)
   fail_flag++;
if(ni0->minusId_y != ni1->minusId_y)
   fail_flag++;
if(ni0->plusId_y != ni1->plusId_y)
   fail_flag++;
if(ni0->minusId_z != ni1->minusId_z)
   fail_flag++;
if(ni0->plusId_z != ni1->plusId_z)
   fail_flag++;
if(ni0->procId_x != ni1->procId_x)
   fail_flag++;
if(ni0->procId_y != ni1->procId_y)
   fail_flag++;
if(ni0->procId_z != ni1->procId_z)
   fail_flag++;
if(ni0->globnx != ni1->globnx)
   fail_flag++;
if(ni0->nx1 != ni1->nx1)
   fail_flag++;
if(ni0->nx2 != ni1->nx2)
   fail_flag++;
if(ni0->ixminus != ni1->ixminus)
   fail_flag++;
if(ni0->ixplus != ni1->ixplus)
   fail_flag++;
if(ni0->loc_nx != ni1->loc_nx)
   fail_flag++;
if(ni0->globny != ni1->globny)
   fail_flag++;
if(ni0->ny1 != ni1->ny1)
   fail_flag++;
if(ni0->ny2 != ni1->ny2)
   fail_flag++;
if(ni0->iyminus != ni1->iyminus)
   fail_flag++;
if(ni0->iyplus != ni1->iyplus)
   fail_flag++;
if(ni0->loc_ny != ni1->loc_ny)
   fail_flag++;
if(ni0->globnz != ni1->globnz)
   fail_flag++;
if(ni0->nz1 != ni1->nz1)
   fail_flag++;
if(ni0->nz2 != ni1->nz2)
   fail_flag++;
if(ni0->izminus != ni1->izminus)
   fail_flag++;
if(ni0->izplus != ni1->izplus)
   fail_flag++;
if(ni0->loc_nz != ni1->loc_nz)
   fail_flag++;

if(rp0->nx != rp1->nx)
   fail_flag++;
if(rp0->ny != rp1->ny)
   fail_flag++;
if(rp0->nz != rp1->nz)
   fail_flag++;
if(rp0->freesurf != rp1->freesurf)
   fail_flag++;
if(rp0->geoproj != rp1->geoproj)
   fail_flag++;

/* these are floats, need better way to check them accurately 20111115
if(rp0->h != rp1->h)
   fail_flag++;
if(rp0->dt != rp1->dt)
   fail_flag++;
if(rp0->modelrot != rp1->modelrot)
   fail_flag++;
if(rp0->modellon != rp1->modellon)
   fail_flag++;
if(rp0->modellat != rp1->modellat)
   fail_flag++;
if(rp0->xmom != rp1->xmom)
   fail_flag++;
if(rp0->ymom != rp1->ymom)
   fail_flag++;
if(rp0->zmom != rp1->zmom)
   fail_flag++;
if(rp0->tdelay != rp1->tdelay)
   fail_flag++;
*/

if(fail_flag != 0)
   {
   fprintf(stderr,"Problem with input parameters in restart file:\n");
   fprintf(stderr,"   %s\n",file);
   fprintf(stderr,"     fail_flag= %d\n\n",fail_flag);
   print_rparsP3(rp0,rp1);
   fprintf(stderr,"\n    exiting ...\n");
   fflush(stderr);
   mpi_exit(-1);
   }
}

void get_restartP3(struct restartinfoP3 *rst,float *pvf,int nx,int ny,int nz)
{
off_t blen;

blen = (off_t)(N_WAVE_VARS*nx*ny*nz*sizeof(float));
reed(rst->fdr,pvf,blen);
close(rst->fdr);
}

void dump_restartP3(struct restartinfoP3 *rst,char *name,int sid,struct runparamsP3 *rpars,struct modelstorage *mds,float *pvf,int it)
{
float *pvf_ptr;
int nx, ny, nz;
int iy, rmflag;
char string[2048];
off_t blen;
int fdr, fdw;
char *cbuf;

nx = rpars->nx;
ny = rpars->ny;
nz = rpars->nz;

rst->it = it;

if(strcmp(rst->dumpfile,"NOT_STARTED_YET") == 0)
   {
   sprintf(rst->dumpfile,"%s/%s_rst-%.5d.e3d",rst->dir,name,sid);
   rst->fdw = croptrfile(rst->dumpfile);

   sprintf(rst->bkupfile,"%s.backup",rst->dumpfile);
   rmflag = 0;
   }
else
   { /* RWG 06/28/07 NEW WAY: use reed() and rite() pairs to copy files */
   fcp_read_write(rst->dumpfile,rst->bkupfile);

   rmflag = 1;
   }

fprintf(stderr,"**** Dumping wavefield to restart file:\n");
fprintf(stderr,"     %s\n",rst->dumpfile);
fflush(stderr);

lseek(rst->fdw,0,SEEK_SET);
rite(rst->fdw,&rst->it,sizeof(int));
rite(rst->fdw,rpars,sizeof(struct runparamsP3));

blen = (off_t)(N_WAVE_VARS*nx*ny*nz*sizeof(float));
rite(rst->fdw,pvf,blen);

/* RWG 06/28/07 NEW WAY: use unlink() to remove file */

if(rmflag != 0)
   unlink(rst->bkupfile);
}

void print_rparsP3(struct runparamsP3 *rp0,struct runparamsP3 *rp1)
{
struct nodeinfo *ni0, *ni1;

ni0 = &(rp0->ni);
ni1 = &(rp1->ni);

fprintf(stderr,"ni0->nproc= %d\tni1->nproc= %d\n",ni0->nproc,ni1->nproc);
fprintf(stderr,"ni0->nproc_x= %d\tni1->nproc_x= %d\n",ni0->nproc_x,ni1->nproc_x);
fprintf(stderr,"ni0->nproc_y= %d\tni1->nproc_y= %d\n",ni0->nproc_y,ni1->nproc_y);
fprintf(stderr,"ni0->nproc_z= %d\tni1->nproc_z= %d\n",ni0->nproc_z,ni1->nproc_z);
fprintf(stderr,"ni0->min_nproc= %d\tni1->min_nproc= %d\n",ni0->min_nproc,ni1->min_nproc);
fprintf(stderr,"ni0->segmentId= %d\tni1->segmentId= %d\n",ni0->segmentId,ni1->segmentId);
fprintf(stderr,"ni0->minusId_x= %d\tni1->minusId_x= %d\n",ni0->minusId_x,ni1->minusId_x);
fprintf(stderr,"ni0->plusId_x= %d\tni1->plusId_x= %d\n",ni0->plusId_x,ni1->plusId_x);
fprintf(stderr,"ni0->minusId_y= %d\tni1->minusId_y= %d\n",ni0->minusId_y,ni1->minusId_y);
fprintf(stderr,"ni0->plusId_y= %d\tni1->plusId_y= %d\n",ni0->plusId_y,ni1->plusId_y);
fprintf(stderr,"ni0->minusId_z= %d\tni1->minusId_z= %d\n",ni0->minusId_z,ni1->minusId_z);
fprintf(stderr,"ni0->plusId_z= %d\tni1->plusId_z= %d\n",ni0->plusId_z,ni1->plusId_z);
fprintf(stderr,"ni0->procId_x= %d\tni1->procId_x= %d\n",ni0->procId_x,ni1->procId_x);
fprintf(stderr,"ni0->procId_y= %d\tni1->procId_y= %d\n",ni0->procId_y,ni1->procId_y);
fprintf(stderr,"ni0->procId_z= %d\tni1->procId_z= %d\n",ni0->procId_z,ni1->procId_z);
fprintf(stderr,"ni0->globnx= %d\tni1->globnx= %d\n",ni0->globnx,ni1->globnx);
fprintf(stderr,"ni0->nx1= %d\tni1->nx1= %d\n",ni0->nx1,ni1->nx1);
fprintf(stderr,"ni0->nx2= %d\tni1->nx2= %d\n",ni0->nx2,ni1->nx2);
fprintf(stderr,"ni0->ixminus= %d\tni1->ixminus= %d\n",ni0->ixminus,ni1->ixminus);
fprintf(stderr,"ni0->ixplus= %d\tni1->ixplus= %d\n",ni0->ixplus,ni1->ixplus);
fprintf(stderr,"ni0->loc_nx= %d\tni1->loc_nx= %d\n",ni0->loc_nx,ni1->loc_nx);
fprintf(stderr,"ni0->globny= %d\tni1->globny= %d\n",ni0->globny,ni1->globny);
fprintf(stderr,"ni0->ny1= %d\tni1->ny1= %d\n",ni0->ny1,ni1->ny1);
fprintf(stderr,"ni0->ny2= %d\tni1->ny2= %d\n",ni0->ny2,ni1->ny2);
fprintf(stderr,"ni0->iyminus= %d\tni1->iyminus= %d\n",ni0->iyminus,ni1->iyminus);
fprintf(stderr,"ni0->iyplus= %d\tni1->iyplus= %d\n",ni0->iyplus,ni1->iyplus);
fprintf(stderr,"ni0->loc_ny= %d\tni1->loc_ny= %d\n",ni0->loc_ny,ni1->loc_ny);
fprintf(stderr,"ni0->globnz= %d\tni1->globnz= %d\n",ni0->globnz,ni1->globnz);
fprintf(stderr,"ni0->nz1= %d\tni1->nz1= %d\n",ni0->nz1,ni1->nz1);
fprintf(stderr,"ni0->nz2= %d\tni1->nz2= %d\n",ni0->nz2,ni1->nz2);
fprintf(stderr,"ni0->izminus= %d\tni1->izminus= %d\n",ni0->izminus,ni1->izminus);
fprintf(stderr,"ni0->izplus= %d\tni1->izplus= %d\n",ni0->izplus,ni1->izplus);
fprintf(stderr,"ni0->loc_nz= %d\tni1->loc_nz= %d\n",ni0->loc_nz,ni1->loc_nz);

fprintf(stderr,"rp0->nx= %d\trp1->nx= %d\n",rp0->nx,rp1->nx);
fprintf(stderr,"rp0->ny= %d\trp1->ny= %d\n",rp0->ny,rp1->ny);
fprintf(stderr,"rp0->nz= %d\trp1->nz= %d\n",rp0->nz,rp1->nz);
fprintf(stderr,"rp0->nt= %d\trp1->nt= %d\n",rp0->nt,rp1->nt);
fprintf(stderr,"rp0->freesurf= %d\trp1->freesurf= %d\n",rp0->freesurf,rp1->freesurf);
fprintf(stderr,"rp0->h= %13.5e\trp1->h= %13.5e\n",rp0->h,rp1->h);
fprintf(stderr,"rp0->dt= %13.5e\trp1->dt= %13.5e\n",rp0->dt,rp1->dt);
fprintf(stderr,"rp0->modelrot= %13.5e\trp1->modelrot= %13.5e\n",rp0->modelrot,rp1->modelrot);
fprintf(stderr,"rp0->modellon= %13.5e\trp1->modellon= %13.5e\n",rp0->modellon,rp1->modellon);
fprintf(stderr,"rp0->modellat= %13.5e\trp1->modellat= %13.5e\n",rp0->modellat,rp1->modellat);
fprintf(stderr,"rp0->xmom= %13.5e\trp1->xmom= %13.5e\n",rp0->xmom,rp1->xmom);
fprintf(stderr,"rp0->ymom= %13.5e\trp1->ymom= %13.5e\n",rp0->ymom,rp1->ymom);
fprintf(stderr,"rp0->zmom= %13.5e\trp1->zmom= %13.5e\n",rp0->zmom,rp1->zmom);
fprintf(stderr,"rp0->tdelay= %13.5e\trp1->tdelay= %13.5e\n",rp0->tdelay,rp1->tdelay);
}

void wroutP3(struct outputfields *outp,float *pvf,float *medf,int iy,int it,struct runparamsP3 *rpars)
{
struct tsoutputparams *tsoutptr;
struct seisoutputparams *soutptr;
struct seisheader *shptr;
struct sgtoutputparams *sgtoutptr;
struct sgtmaster *sgtmastptr;
struct sgtheader *sgtheadptr;
float *pin, *pout;
int i, j, np, ip;
off_t off, bblen, nrite;
int its, kp;
float den, num;
float one = 1.0;
float onep5 = 1.5;
float two = 2.0;
float three = 3.0;
float quart = 0.25;
float normf = 1.0e-05;

int nx, nz;
int glob_iy;
struct nodeinfo *ni;
int nseis_comps = 3;

/* change 3 components to 9 components for nseis output */
nseis_comps = 9;

ni = &(rpars->ni);
nx = rpars->nx;
nz = rpars->nz;

np = nx*nz;
glob_iy = ni->ny1 + iy;

if(outp->nseis)
   {
   soutptr = &(outp->spnts);

   if(soutptr->iflag)  /* all_in_one output */
      {
      if(soutptr->flushbuf == 1)
         {
         bblen = nseis_comps*(soutptr->np)*(soutptr->nbufrite)*sizeof(float);
         off = sizeof(int) + (soutptr->np)*sizeof(struct seisheader) + bblen
                        - soutptr->cur_off;
         lseek(soutptr->fdw,off,SEEK_CUR);

         bblen = nseis_comps*(soutptr->np)*(soutptr->buflen)*sizeof(float);
         nrite = rite(soutptr->fdw,soutptr->s,bblen);

         soutptr->cur_off = soutptr->cur_off + off + nrite;
         soutptr->nbufrite = soutptr->nbufrite + soutptr->buflen;
         soutptr->buflen = 0;
         }
      else
         {
         if(it%(soutptr->tinc) == 0)
            {
            shptr = soutptr->shead;
            for(j=0;j<(soutptr->np);j++)
               {
               if(iy == shptr[j].iy)
                  {
                  if(j == 0) /* increment buffer counter */
                     soutptr->buflen = soutptr->buflen + 1;

                  ip = shptr[j].iz*nx + shptr[j].ix;

                  its = it/(soutptr->tinc) - (soutptr->nbufrite);
                  kp = nseis_comps*(soutptr->np)*its + nseis_comps*j;

                  for(i=0;i<nseis_comps;i++)
                     soutptr->s[i+kp] = pvf[i*np+ip];
                  }
               }
            }
         }
      }
   }   

tsoutptr = &(outp->xyslice);
if(tsoutptr->iflag)
   {
   if(tsoutptr->flushbuf == 1 && tsoutptr->buflen)
      {
      bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*sizeof(float);
      off = bblen + tsoutptr->head_off - tsoutptr->cur_off;
      lseek(tsoutptr->fdw,off,SEEK_CUR);

      bblen = 3*(tsoutptr->np)*(tsoutptr->buflen)*sizeof(float);
      nrite = rite(tsoutptr->fdw,tsoutptr->vbuf,bblen);

      tsoutptr->cur_off = tsoutptr->cur_off + off + nrite;
      tsoutptr->nbufrite = tsoutptr->nbufrite + tsoutptr->buflen;
      tsoutptr->buflen = 0;
      }
   else
      {
      if(glob_iy >= tsoutptr->iyleft && glob_iy <= tsoutptr->iyright && glob_iy%(tsoutptr->idy) == 0 && it%(tsoutptr->idt) == 0 && it >= 0)
         {
	 its = it/(tsoutptr->idt) - (tsoutptr->nbufrite);
	 tsoutptr->buflen = 1 + its;

         for(i=0;i<3;i++)
            {
            pin = pvf + i*np + tsoutptr->idz*nx;
	    pout = tsoutptr->vbuf + (3*its + i)*(tsoutptr->np) + tsoutptr->nx*((glob_iy - tsoutptr->iyleft)/tsoutptr->idy);

            copy(pin,pout,tsoutptr->nx,tsoutptr->idx);
            }
         }
      }
   }

tsoutptr = &(outp->xzslice);
if(tsoutptr->iflag)
   {
   if(tsoutptr->flushbuf == 1 && tsoutptr->buflen)
      {
      bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*sizeof(float);
      off = bblen + tsoutptr->head_off - tsoutptr->cur_off;
      lseek(tsoutptr->fdw,off,SEEK_CUR);

      bblen = 3*(tsoutptr->np)*(tsoutptr->buflen)*sizeof(float);
      nrite = rite(tsoutptr->fdw,tsoutptr->vbuf,bblen);

      tsoutptr->cur_off = tsoutptr->cur_off + off + nrite;
      tsoutptr->nbufrite = tsoutptr->nbufrite + tsoutptr->buflen;
      tsoutptr->buflen = 0;
      }
   else
      {
      if(glob_iy == tsoutptr->idy && it%(tsoutptr->idt) == 0 && it >= 0)
         {
         its = it/(tsoutptr->idt) - (tsoutptr->nbufrite);
         tsoutptr->buflen = 1 + its;

         for(i=0;i<3;i++)
            {
            for(j=0;j<tsoutptr->nz;j++)
               {
               pin = pvf + i*np + (tsoutptr->iz0 + j*tsoutptr->idz)*nx;
               pout = tsoutptr->vbuf + (3*its + i)*(tsoutptr->np) + j*tsoutptr->nx;
               copy(pin,pout,tsoutptr->nx,tsoutptr->idx);
               }
            }
         }
      }
   }

tsoutptr = &(outp->yzslice);
if(tsoutptr->iflag)
   {
   if(tsoutptr->flushbuf == 1 && tsoutptr->buflen)
      {
      bblen = 3*(tsoutptr->np)*(tsoutptr->nbufrite)*sizeof(float);
      off = bblen + tsoutptr->head_off - tsoutptr->cur_off;
      lseek(tsoutptr->fdw,off,SEEK_CUR);

      bblen = 3*(tsoutptr->np)*(tsoutptr->buflen)*sizeof(float);
      nrite = rite(tsoutptr->fdw,tsoutptr->vbuf,bblen);

      tsoutptr->cur_off = tsoutptr->cur_off + off + nrite;
      tsoutptr->nbufrite = tsoutptr->nbufrite + tsoutptr->buflen;
      tsoutptr->buflen = 0;
      }
   else
      {
      if(glob_iy >= tsoutptr->iyleft && glob_iy <= tsoutptr->iyright && glob_iy%(tsoutptr->idy) == 0 && it%(tsoutptr->idt) == 0 && it >= 0)
         {
         its = it/(tsoutptr->idt) - (tsoutptr->nbufrite);
         tsoutptr->buflen = 1 + its;

         for(i=0;i<3;i++)
            {
            pin = pvf + i*np + tsoutptr->idx;
            pout = tsoutptr->vbuf + (3*its + i)*(tsoutptr->np) + tsoutptr->nz*((glob_iy - tsoutptr->iyleft)/tsoutptr->idy);

            copy(pin,pout,tsoutptr->nz,tsoutptr->idz*nx);
            }
         }
      }
   }

sgtoutptr = &(outp->sgtpnts);
if(sgtoutptr->iflag)
   {
   sgtmastptr = &(sgtoutptr->sgtmast);
   if(sgtoutptr->flushbuf == 1)
      {
      bblen = 6*(sgtmastptr->localnp)*(sgtoutptr->nbufrite)*sizeof(float);
      off = sizeof(struct sgtmaster) + (sgtmastptr->globnp)*sizeof(struct sgtindex)
             + (sgtmastptr->localnp)*sizeof(struct sgtheader) + bblen - sgtoutptr->cur_off;
      lseek(sgtoutptr->fdw,off,SEEK_CUR);

      bblen = 6*(sgtmastptr->localnp)*(sgtoutptr->buflen)*sizeof(float);
      nrite = rite(sgtoutptr->fdw,sgtoutptr->sgt,bblen);

      sgtoutptr->cur_off = sgtoutptr->cur_off + off + nrite;
      sgtoutptr->nbufrite = sgtoutptr->nbufrite + sgtoutptr->buflen;
      sgtoutptr->buflen = 0;
      }
   else
      {
      if(it%(sgtoutptr->tinc) == 0)
         {
         sgtheadptr = sgtoutptr->sgthead;
         for(j=0;j<(sgtmastptr->localnp);j++)
            {
            if(glob_iy == sgtheadptr[j].ysgt)  /* recall xsgt,ysgt,zsgt are global indices */
               {
	       if(j == 0) /* increment buffer counter */
                  sgtoutptr->buflen = sgtoutptr->buflen + 1;

	       ip = (sgtheadptr[j].zsgt-ni->nz1)*nx + (sgtheadptr[j].xsgt-ni->nx1);

	       num = two*(medf[11*np+ip] + medf[12*np+ip])/medf[11*np+ip];
	       den = normf*quart*medf[11*np+ip]/(medf[12*np+ip]*(onep5*medf[11*np+ip] + medf[12*np+ip]));

	       if(it == 0) /* set lambda,rigidity,density and rewrite header */
	          {
	          sgtheadptr[j].lam = (1.0e+10)*medf[11*np+ip];
	          sgtheadptr[j].mu = (1.0e+10)*medf[12*np+ip];
	          sgtheadptr[j].rho = three/(medf[5*np+ip] + medf[6*np+ip] + medf[7*np+ip]);

                  off = sizeof(struct sgtmaster) + (sgtmastptr->globnp)*sizeof(struct sgtindex)
		         + j*sizeof(struct sgtheader) - sgtoutptr->cur_off;
                  lseek(sgtoutptr->fdw,off,SEEK_CUR);

	          nrite = rite(sgtoutptr->fdw,&sgtheadptr[j],sizeof(struct sgtheader));
                  sgtoutptr->cur_off = sgtoutptr->cur_off + off + nrite;
	          }

	       its = it/(sgtoutptr->tinc) - (sgtoutptr->nbufrite);
               kp = 6*(sgtmastptr->localnp)*its + 6*j;

/*
            Note that the strains are NOT normalized by the factor
            of (h).  This factor would balance the (1/h) factor
	    of the body force used to specify the moment.  Since these
	    factors cancel, neither is applied.  See main.c
*/

	       /* mxx */
	       sgtoutptr->sgt[kp]   = den*(num*pvf[3*np+ip]
				        - pvf[4*np+ip] - pvf[5*np+ip]);
	       /* myy */
	       sgtoutptr->sgt[1+kp] = den*(num*pvf[4*np+ip]
				        - pvf[3*np+ip] - pvf[5*np+ip]);
	       /* mzz */
	       sgtoutptr->sgt[2+kp] = den*(num*pvf[5*np+ip]
				        - pvf[3*np+ip] - pvf[4*np+ip]);

	       /* mxy */
	       sgtoutptr->sgt[3+kp] = normf*pvf[6*np+ip]/medf[12*np+ip];
	       /* mxz */
	       sgtoutptr->sgt[4+kp] = normf*pvf[7*np+ip]/medf[12*np+ip];
	       /* myz */
	       sgtoutptr->sgt[5+kp] = normf*pvf[8*np+ip]/medf[12*np+ip];
               }
            }
         }
      }
   }
}
