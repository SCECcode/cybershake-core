#include "include.h"
#include "defs.h"
#include "structure.h"
#include "srf_structure.h"
#include "function.h"

void init_plane_srf(struct standrupformat *srf,float *elon,float *elat,int nx,int ny,float *fl,float *fw,float *dx,float *dy,float *stk,float *dip,float *dtop,float *sh,float *dh)
{
struct srf_prectsegments *prseg_ptr;
int ig;

sprintf(srf[0].type,"PLANE");

srf[0].srf_prect.nseg = 1;
srf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));
prseg_ptr = srf[0].srf_prect.prectseg;

prseg_ptr[0].elon = *elon;
prseg_ptr[0].elat = *elat;
prseg_ptr[0].nstk = nx;
prseg_ptr[0].ndip = ny;
prseg_ptr[0].flen = *fl;
prseg_ptr[0].fwid = *fw;
prseg_ptr[0].dlen = *dx;
prseg_ptr[0].dwid = *dy;
prseg_ptr[0].stk = *stk;
prseg_ptr[0].dip = *dip;
prseg_ptr[0].dtop = *dtop;
prseg_ptr[0].shyp = *sh;
prseg_ptr[0].dhyp = *dh;

srf[0].srf_apnts.np = 0;
for(ig=0;ig<srf[0].srf_prect.nseg;ig++)
   srf[0].srf_apnts.np = srf[0].srf_apnts.np + (prseg_ptr[ig].nstk)*(prseg_ptr[ig].ndip);

srf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));
}

void load_slip_srf(struct standrupformat *srf,struct stfpar *spar,struct pointsource *ps)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
float area;
float *stf;
int i, j, ip, ig, ntot;

float dmin = 4.0;
float dmax = 6.0;
float rtfac, tzero;
float rtfac0 = 1.0;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

for(ip=0;ip<apnts_ptr->np;ip++)
   {
   apval_ptr[ip].stf1 = (float *)check_malloc(spar->nt*sizeof(float));
   stf = apval_ptr[ip].stf1;

   apval_ptr[ip].dt = spar->dt;
   if(ps[ip].slip > MINSLIP)
      {
      if(ps[ip].dep >= dmax)
         rtfac = 1.0;
      else if(ps[ip].dep < dmax && ps[ip].dep > dmin)
         rtfac = 1.0 + rtfac0*(dmax-(ps[ip].dep))/(dmax-dmin);
      else
         rtfac = 1.0 + rtfac0;

      tzero = rtfac*spar->trise;
      apval_ptr[ip].nt1 = gen_2tri_stf(&(ps[ip].slip),&spar->trise,stf,spar->nt,&spar->dt,&ps[ip].dep);
      }
   else
      apval_ptr[ip].nt1 = 0;

   if(apval_ptr[ip].nt1)
      apval_ptr[ip].stf1 = (float *)check_realloc(apval_ptr[ip].stf1,(apval_ptr[ip].nt1)*sizeof(float));
   else
      {
      free(apval_ptr[ip].stf1);
      apval_ptr[ip].stf1 = NULL;
      }

   apval_ptr[ip].lon = ps[ip].lon;
   apval_ptr[ip].lat = ps[ip].lat;
   apval_ptr[ip].dep = ps[ip].dep;
   apval_ptr[ip].stk = ps[ip].stk;
   apval_ptr[ip].dip = ps[ip].dip;
   apval_ptr[ip].area = ps[ip].area;
   apval_ptr[ip].rake = ps[ip].rak;
   apval_ptr[ip].slip1 = ps[ip].slip;

   apval_ptr[ip].slip2 = 0.0;
   apval_ptr[ip].nt2 = 0;
   apval_ptr[ip].stf2 = NULL;
   apval_ptr[ip].slip3 = 0.0;
   apval_ptr[ip].nt3 = 0;
   apval_ptr[ip].stf3 = NULL;
   }
}

void load_rupt_srf(struct standrupformat *srf,struct pointsource *ps,float *sh,float *dh)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
int ip;

(srf->srf_prect).prectseg[0].shyp = *sh;
(srf->srf_prect).prectseg[0].dhyp = *dh;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

for(ip=0;ip<apnts_ptr->np;ip++)
   apval_ptr[ip].tinit = ps[ip].rupt;
}

void write_srf1(struct standrupformat *srf,char *file,int bflag)
{
FILE *fpw, *fopfile();
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float area;
float *stf;
int i, j, k, nt6, it, ip, ig, ntot;

char pword[32];
int fdw;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

if(bflag)
   {
   if(strcmp(file,"stdout") == 0)
      fdw = STDOUT_FILENO;
   else
      fdw = croptrfile(file);

   rite(fdw,srf->version,sizeof(srf->version));

   if(strcmp(srf->type,"PLANE") == 0)
      {
      rite(fdw,srf->type,sizeof(srf->type));
      rite(fdw,&(prect_ptr->nseg),sizeof(prect_ptr->nseg));
      rite(fdw,prseg_ptr,(prect_ptr->nseg)*sizeof(struct srf_prectsegments));
      }

   sprintf(pword,"POINTS");
   rite(fdw,pword,sizeof(pword));
   rite(fdw,&(apnts_ptr->np),sizeof(apnts_ptr->np));
   for(i=0;i<apnts_ptr->np;i++)
      {
      rite(fdw,&(apval_ptr[i].lon),sizeof(float));
      rite(fdw,&(apval_ptr[i].lat),sizeof(float));
      rite(fdw,&(apval_ptr[i].dep),sizeof(float));
      rite(fdw,&(apval_ptr[i].stk),sizeof(float));
      rite(fdw,&(apval_ptr[i].dip),sizeof(float));
      rite(fdw,&(apval_ptr[i].area),sizeof(float));
      rite(fdw,&(apval_ptr[i].tinit),sizeof(float));
      rite(fdw,&(apval_ptr[i].dt),sizeof(float));
      rite(fdw,&(apval_ptr[i].rake),sizeof(float));
      rite(fdw,&(apval_ptr[i].slip1),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt1),sizeof(int));
      rite(fdw,&(apval_ptr[i].slip2),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt2),sizeof(int));
      rite(fdw,&(apval_ptr[i].slip3),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt3),sizeof(int));

      rite(fdw,apval_ptr[i].stf1,(apval_ptr[i].nt1)*sizeof(float));
      rite(fdw,apval_ptr[i].stf2,(apval_ptr[i].nt2)*sizeof(float));
      rite(fdw,apval_ptr[i].stf3,(apval_ptr[i].nt3)*sizeof(float));
      }
   close(fdw);
   }
else
   {
   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      fpw = fopfile(file,"w");

   fprintf(fpw,"%s\n",srf->version);

   if(strcmp(srf->type,"PLANE") == 0)
      {
      fprintf(fpw,"%s %d\n",srf->type,prect_ptr->nseg);
      for(ig=0;ig<prect_ptr->nseg;ig++)
         {
         fprintf(fpw,"%10.4f %9.4f %5d %5d %8.2f %8.2f\n",prseg_ptr[ig].elon,
                                                        prseg_ptr[ig].elat,
                                                        prseg_ptr[ig].nstk,
                                                        prseg_ptr[ig].ndip,
                                                        prseg_ptr[ig].flen,
                                                        prseg_ptr[ig].fwid);
         fprintf(fpw,"%4.0f %4.0f %8.2f %8.2f %8.2f\n",prseg_ptr[ig].stk,
                                                    prseg_ptr[ig].dip,
                                                    prseg_ptr[ig].dtop,
                                                    prseg_ptr[ig].shyp,
                                                    prseg_ptr[ig].dhyp);
         }
      }

   fprintf(fpw,"POINTS %d\n",apnts_ptr->np);
   for(i=0;i<apnts_ptr->np;i++)
      {
      fprintf(fpw,"%10.4f %9.4f %9.4f %4.0f %4.0f %12.5e %10.4f %12.5e\n",
                                              apval_ptr[i].lon,
                                              apval_ptr[i].lat,
                                              apval_ptr[i].dep,
                                              apval_ptr[i].stk,
                                              apval_ptr[i].dip,
                                              apval_ptr[i].area,
                                              apval_ptr[i].tinit,
                                              apval_ptr[i].dt);
      fprintf(fpw,"%4.0f %8.2f %6d %8.2f %6d %8.2f %6d\n",
                                              apval_ptr[i].rake,
                                              apval_ptr[i].slip1,
                                              apval_ptr[i].nt1,
                                              apval_ptr[i].slip2,
                                              apval_ptr[i].nt2,
                                              apval_ptr[i].slip3,
                                              apval_ptr[i].nt3);

      stf = apval_ptr[i].stf1;
      nt6 = (apval_ptr[i].nt1)/6;
      for(k=0;k<nt6;k++)
         {
         for(j=0;j<6;j++)
            {
            it = 6*k + j;
            fprintf(fpw,"%13.5e",stf[it]);
            }
         fprintf(fpw,"\n");
         }

      if(6*nt6 != (apval_ptr[i].nt1))
         {
         for(j=6*nt6;j<(apval_ptr[i].nt1);j++)
            fprintf(fpw,"%13.5e",stf[j]);

         fprintf(fpw,"\n");
         }

      stf = apval_ptr[i].stf2;
      nt6 = (apval_ptr[i].nt2)/6;
      for(k=0;k<nt6;k++)
         {
         for(j=0;j<6;j++)
            {
            it = 6*k + j;
            fprintf(fpw,"%13.5e",stf[it]);
            }
         fprintf(fpw,"\n");
         }

      if(6*nt6 != (apval_ptr[i].nt2))
         {
         for(j=6*nt6;j<(apval_ptr[i].nt2);j++)
            fprintf(fpw,"%13.5e",stf[j]);

         fprintf(fpw,"\n");
         }

      stf = apval_ptr[i].stf3;
      nt6 = (apval_ptr[i].nt3)/6;
      for(k=0;k<nt6;k++)
         {
         for(j=0;j<6;j++)
            {
            it = 6*k + j;
            fprintf(fpw,"%13.5e",stf[it]);
            }
         fprintf(fpw,"\n");
         }

      if(6*nt6 != (apval_ptr[i].nt3))
         {
         for(j=6*nt6;j<(apval_ptr[i].nt3);j++)
            fprintf(fpw,"%13.5e",stf[j]);

         fprintf(fpw,"\n");
         }
      }
   fclose(fpw);
   }
}

#ifndef _V3_3_1
void write_srf2(struct standrupformat *srf,char *file,int bflag)
{
FILE *fpw, *_fopfile();
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float area;
float *stf;
int i, j, k, nt6, it, ig;

char *sptr, pword[32];
int fdw;

int ip, np_rite;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

if(bflag)
   {
   if(strcmp(file,"stdout") == 0)
      fdw = STDOUT_FILENO;
   else
      fdw = _croptrfile(file);

   _rite(fdw,srf->version,sizeof(srf->version));

   if(strcmp(srf->type,"PLANE") == 0)
      {
      _rite(fdw,srf->type,sizeof(srf->type));
      _rite(fdw,&(prect_ptr->nseg),sizeof(prect_ptr->nseg));
      _rite(fdw,prseg_ptr,(prect_ptr->nseg)*sizeof(struct srf_prectsegments));
      }

   sprintf(pword,"POINTS");
   _rite(fdw,pword,sizeof(pword));
   _rite(fdw,&(apnts_ptr->np),sizeof(apnts_ptr->np));
   for(i=0;i<apnts_ptr->np;i++)
      {
      _rite(fdw,&(apval_ptr[i].lon),sizeof(float));
      _rite(fdw,&(apval_ptr[i].lat),sizeof(float));
      _rite(fdw,&(apval_ptr[i].dep),sizeof(float));
      _rite(fdw,&(apval_ptr[i].stk),sizeof(float));
      _rite(fdw,&(apval_ptr[i].dip),sizeof(float));
      _rite(fdw,&(apval_ptr[i].area),sizeof(float));
      _rite(fdw,&(apval_ptr[i].tinit),sizeof(float));
      _rite(fdw,&(apval_ptr[i].dt),sizeof(float));
      _rite(fdw,&(apval_ptr[i].rake),sizeof(float));
      _rite(fdw,&(apval_ptr[i].slip1),sizeof(float));
      _rite(fdw,&(apval_ptr[i].nt1),sizeof(int));
      _rite(fdw,&(apval_ptr[i].slip2),sizeof(float));
      _rite(fdw,&(apval_ptr[i].nt2),sizeof(int));
      _rite(fdw,&(apval_ptr[i].slip3),sizeof(float));
      _rite(fdw,&(apval_ptr[i].nt3),sizeof(int));

      _rite(fdw,apval_ptr[i].stf1,(apval_ptr[i].nt1)*sizeof(float));
      _rite(fdw,apval_ptr[i].stf2,(apval_ptr[i].nt2)*sizeof(float));
      _rite(fdw,apval_ptr[i].stf3,(apval_ptr[i].nt3)*sizeof(float));
      }
   close(fdw);
   }
else
   {
   if(strcmp(file,"stdout") == 0) {
      fpw = stdout;
   } else {
      fpw = _fopfile(file,"w");
   }

   fprintf(fpw,"%s\n",srf->version);

   for(i=0;i<srf->srf_hcmnt.nline;i++)
      {
      sptr = (srf->srf_hcmnt.cbuf) + i*MAXLINE;
      fprintf(fpw,"%s",sptr);
      }

   if(strcmp(srf->type,"PLANE") == 0)
      {
      fprintf(fpw,"%s %d\n",srf->type,prect_ptr->nseg);
      for(ig=0;ig<prect_ptr->nseg;ig++)
         {
         fprintf(fpw,"%12.6f %11.6f %5d %5d %10.4f %10.4f\n",prseg_ptr[ig].elon,
                                                        prseg_ptr[ig].elat,
                                                        prseg_ptr[ig].nstk,
                                                        prseg_ptr[ig].ndip,
                                                        prseg_ptr[ig].flen,
                                                        prseg_ptr[ig].fwid);
         fprintf(fpw,"%4.0f %4.0f %10.4f %10.4f %10.4f\n",prseg_ptr[ig].stk,
                                                    prseg_ptr[ig].dip,
                                                    prseg_ptr[ig].dtop,
                                                    prseg_ptr[ig].shyp,
                                                    prseg_ptr[ig].dhyp);
         }
      }

   np_rite = 0;
   for(ig=0;ig<srf->nseg;ig++)
      {
      fprintf(fpw,"POINTS %d\n",srf->np_seg[ig]);
      for(ip=0;ip<srf->np_seg[ig];ip++)
         {
         i = ip + np_rite;

         fprintf(fpw,"%12.6f %11.6f %12.5e %4.0f %4.0f %12.5e %13.6e %12.5e %13.5e %13.5e\n",
                                              apval_ptr[i].lon,
                                              apval_ptr[i].lat,
                                              apval_ptr[i].dep,
                                              apval_ptr[i].stk,
                                              apval_ptr[i].dip,
                                              apval_ptr[i].area,
                                              apval_ptr[i].tinit,
                                              apval_ptr[i].dt,
                                              apval_ptr[i].vs,
                                              apval_ptr[i].den);

         fprintf(fpw,"%4.0f %10.4f %6d %10.4f %6d %10.4f %6d\n",
                                              apval_ptr[i].rake,
                                              apval_ptr[i].slip1,
                                             apval_ptr[i].nt1,
                                              apval_ptr[i].slip2,
                                              apval_ptr[i].nt2,
                                              apval_ptr[i].slip3,
                                              apval_ptr[i].nt3);

         stf = apval_ptr[i].stf1;
         nt6 = (apval_ptr[i].nt1)/6;
         for(k=0;k<nt6;k++)
            {
            for(j=0;j<6;j++)
               {
               it = 6*k + j;
               fprintf(fpw,"%13.5e",stf[it]);
               }
            fprintf(fpw,"\n");
            }

         if(6*nt6 != (apval_ptr[i].nt1))
            {
            for(j=6*nt6;j<(apval_ptr[i].nt1);j++)
               fprintf(fpw,"%13.5e",stf[j]);

            fprintf(fpw,"\n");
            }

         stf = apval_ptr[i].stf2;
         nt6 = (apval_ptr[i].nt2)/6;
         for(k=0;k<nt6;k++)
            {
            for(j=0;j<6;j++)
               {
               it = 6*k + j;
               fprintf(fpw,"%13.5e",stf[it]);
               }
            fprintf(fpw,"\n");
            }

         if(6*nt6 != (apval_ptr[i].nt2))
            {
            for(j=6*nt6;j<(apval_ptr[i].nt2);j++)
               fprintf(fpw,"%13.5e",stf[j]);

            fprintf(fpw,"\n");
            }

         stf = apval_ptr[i].stf3;
         nt6 = (apval_ptr[i].nt3)/6;
         for(k=0;k<nt6;k++)
            {
            for(j=0;j<6;j++)
               {
               it = 6*k + j;
               fprintf(fpw,"%13.5e",stf[it]);
               }
            fprintf(fpw,"\n");
            }

         if(6*nt6 != (apval_ptr[i].nt3))
            {
            for(j=6*nt6;j<(apval_ptr[i].nt3);j++)
               fprintf(fpw,"%13.5e",stf[j]);

            fprintf(fpw,"\n");
            }
         }
      np_rite = np_rite + srf->np_seg[ig];
      }

   fclose(fpw);
   }
}
#endif

void write_srf(struct standrupformat *srf,char *file,int bflag) {
	printf("Writing SRF to file %s, version %s.\n", file, srf->version);
	if (strcmp(srf->version, "1.0")==0) {
		write_srf1(srf, file, bflag);
	}
	#ifndef _V3_3_1
	 else if (strcmp(srf->version, "2.0")==0) {
		write_srf2(srf, file, bflag);
	}
	#endif
}

/*void free_srf_stf(struct standrupformat *srf)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
int i;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

for(i=0;i<apnts_ptr->np;i++)
   {
   free(apval_ptr[i].stf1);
   free(apval_ptr[i].stf2);
   free(apval_ptr[i].stf3);
   }
}*/

void read_srf(struct standrupformat *srf,char *file,int bflag)
{
FILE *fpr, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
char str[1024];

float *stf;
int i, j, k, nt6, it, ip, ig, ntot;

char pword[32];
int fdr;

printf("Read SRF %s\n", file);

if(bflag)
   {
   if(strcmp(file,"stdin") == 0)
      fdr = STDIN_FILENO;
   else
      fdr = opfile_ro(file);

   reed(fdr,srf->version,sizeof(srf->version));

   reed(fdr,pword,sizeof(pword));
   if(strcmp(pword,"PLANE") == 0)
      {
      sprintf(srf->type,"PLANE");

      reed(fdr,&(srf[0].srf_prect.nseg),sizeof(int));
      srf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));
      prseg_ptr = srf[0].srf_prect.prectseg;

      reed(fdr,prseg_ptr,(srf[0].srf_prect.nseg)*sizeof(struct srf_prectsegments));

      while(strncmp(pword,"POINTS",6) != 0)
         reed(fdr,pword,sizeof(pword));
      }

   if(strncmp(pword,"POINTS",6) == 0)
      {
      reed(fdr,&(srf[0].srf_apnts.np),sizeof(int));
      srf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

      apval_ptr = srf[0].srf_apnts.apntvals;

      for(i=0;i<srf[0].srf_apnts.np;i++)
         {
         reed(fdr,&(apval_ptr[i].lon),sizeof(float));
         reed(fdr,&(apval_ptr[i].lat),sizeof(float));
         reed(fdr,&(apval_ptr[i].dep),sizeof(float));
         reed(fdr,&(apval_ptr[i].stk),sizeof(float));
         reed(fdr,&(apval_ptr[i].dip),sizeof(float));
         reed(fdr,&(apval_ptr[i].area),sizeof(float));
         reed(fdr,&(apval_ptr[i].tinit),sizeof(float));
         reed(fdr,&(apval_ptr[i].dt),sizeof(float));
         reed(fdr,&(apval_ptr[i].rake),sizeof(float));
         reed(fdr,&(apval_ptr[i].slip1),sizeof(float));
         reed(fdr,&(apval_ptr[i].nt1),sizeof(int));
         reed(fdr,&(apval_ptr[i].slip2),sizeof(float));
         reed(fdr,&(apval_ptr[i].nt2),sizeof(int));
         reed(fdr,&(apval_ptr[i].slip3),sizeof(float));
         reed(fdr,&(apval_ptr[i].nt3),sizeof(int));

         apval_ptr[i].stf1 = (float *)check_malloc((apval_ptr[i].nt1)*sizeof(float));
         apval_ptr[i].stf2 = (float *)check_malloc((apval_ptr[i].nt2)*sizeof(float));
         apval_ptr[i].stf3 = (float *)check_malloc((apval_ptr[i].nt3)*sizeof(float));

         reed(fdr,apval_ptr[i].stf1,(apval_ptr[i].nt1)*sizeof(float));
         reed(fdr,apval_ptr[i].stf2,(apval_ptr[i].nt2)*sizeof(float));
         reed(fdr,apval_ptr[i].stf3,(apval_ptr[i].nt3)*sizeof(float));
         }
      }
   close(fdr);
   }
else
   {
   if(strcmp(file,"stdin") == 0)
      fpr = stdin;
   else
      fpr = fopfile(file,"r");

   fgets(str,1024,fpr);
   sscanf(str,"%s",&(srf[0].version));

   fgets(str,1024,fpr);
   sscanf(str,"%s",pword);

   if(strncmp(pword,"PLANE",5) == 0)
      {
      sscanf(str,"%s %d",srf[0].type,&(srf[0].srf_prect.nseg));

      srf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));
      prseg_ptr = srf[0].srf_prect.prectseg;

      for(ig=0;ig<srf[0].srf_prect.nseg;ig++)
         {
         fgets(str,1024,fpr);
         sscanf(str,"%f %f %d %d %f %f",&(prseg_ptr[ig].elon),
                                     &(prseg_ptr[ig].elat),
                                     &(prseg_ptr[ig].nstk),
                                     &(prseg_ptr[ig].ndip),
                                     &(prseg_ptr[ig].flen),
                                     &(prseg_ptr[ig].fwid));
         fgets(str,1024,fpr);
         sscanf(str,"%f %f %f %f %f",&(prseg_ptr[ig].stk),
                                  &(prseg_ptr[ig].dip),
                                  &(prseg_ptr[ig].dtop),
                                  &(prseg_ptr[ig].shyp),
                                  &(prseg_ptr[ig].dhyp));
         }

      while(strncmp(pword,"POINTS",6) != 0)
	 {
         fgets(str,1024,fpr);
         sscanf(str,"%s",pword);
	 }
      }

   if(strncmp(pword,"POINTS",6) == 0)
      {
      sscanf(str,"%*s %d",&(srf[0].srf_apnts.np));
      srf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

      apval_ptr = srf[0].srf_apnts.apntvals;

      for(i=0;i<srf[0].srf_apnts.np;i++)
         {
         fgets(str,1024,fpr);
         sscanf(str,"%f %f %f %f %f %f %f %f",&(apval_ptr[i].lon),
                                           &(apval_ptr[i].lat),
                                           &(apval_ptr[i].dep),
                                           &(apval_ptr[i].stk),
                                           &(apval_ptr[i].dip),
                                           &(apval_ptr[i].area),
                                           &(apval_ptr[i].tinit),
                                           &(apval_ptr[i].dt));
         fgets(str,1024,fpr);
         sscanf(str,"%f %f %d %f %d %f %d",&(apval_ptr[i].rake),
                                        &(apval_ptr[i].slip1),
                                        &(apval_ptr[i].nt1),
                                        &(apval_ptr[i].slip2),
                                        &(apval_ptr[i].nt2),
                                        &(apval_ptr[i].slip3),
                                        &(apval_ptr[i].nt3));

         apval_ptr[i].stf1 = (float *)check_malloc((apval_ptr[i].nt1)*sizeof(float));
         stf = apval_ptr[i].stf1;

         for(it=0;it<(apval_ptr[i].nt1);it++)
            fscanf(fpr,"%f",&stf[it]);

         apval_ptr[i].stf2 = (float *)check_malloc((apval_ptr[i].nt2)*sizeof(float));
         stf = apval_ptr[i].stf2;

         for(it=0;it<(apval_ptr[i].nt2);it++)
            fscanf(fpr,"%f",&stf[it]);

         apval_ptr[i].stf3 = (float *)check_malloc((apval_ptr[i].nt3)*sizeof(float));
         stf = apval_ptr[i].stf3;

         for(it=0;it<(apval_ptr[i].nt3);it++)
            fscanf(fpr,"%f",&stf[it]);

         /* get rouge newline character */
         if((apval_ptr[i].nt1) || (apval_ptr[i].nt2) || (apval_ptr[i].nt3))
            fgets(str,1024,fpr);
         }
      }
   fclose(fpr);
   }
}

int gen_2tri_stfOLD(float *slip,float *trise,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
int ip, it0, it1, it2;
float tr, amp, a0;
float sum;

float alpha = 0.1;      /* 1st triangle has pulse width = 2*alpha*trise */
float betadeep = 0.2;       /* 2nd triangle has amplitude = beta*A (z0>dmax)*/
float betashal = 0.5;       /* 2nd triangle has amplitude = beta*A (z0<dmin)*/
float beta, dbdd;

float dmin = 4.0;
float dmax = 8.0;
float rtfac;
float rtfac0 = 1.0;

dbdd = (betadeep - betashal)/(dmax-dmin);

if((*z0) >= dmax)
   beta = betadeep;
else if((*z0) < dmax && (*z0) > dmin)
   beta = betadeep - (dmax-(*z0))*dbdd;
else
   beta = betashal;

zapit(stf,nt);

if((*z0) >= dmax)
   rtfac = 1.0;
else if((*z0) < dmax && (*z0) > dmin)
   rtfac = 1.0 + rtfac0*(dmax-(*z0))/(dmax-dmin);
else
   rtfac = 1.0 + rtfac0;

tr = (*trise)*rtfac;

alpha = alpha*tr;

it0 = (int)((alpha)/(*dt) + 0.5);
if(it0 < 2)
   it0 = 2;
it1 = (int)((tr)/(*dt) + 0.5);
if(it1 < 4)
   it1 = 4;

it2 = (2 - beta)*it0;

a0 = 1.0;
amp = a0/(float)(it0);

for(it=0;it<it0;it++)
   stf[it] = it*amp;

for(it=it0;it<it2;it++)
   stf[it] = (2*it0-it)*amp;

amp = beta*a0/(float)(it1-it2);

for(it=it2;it<it1;it++)
   stf[it] = beta*a0 + (it2-it)*amp;

nstf = nt-1;
while(stf[nstf] == (float)(0.0) && nstf)
   nstf--;

if(nstf == 0)
   return(0);

if(nstf < nt-1)
   nstf = nstf + 2;;

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_brune_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
float t95, tend, sfac, tfac;
float sum;

zapit(stf,nt);

t95 = 1.745*exp(1.0)*(*t0);
tend = 3.0*t95;

nstf = (int)((tend)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

sfac = (*slip)/(*t0);
tfac = (*dt)/(*t0);
for(it=0;it<nstf;it++)
   stf[it] = sfac*(it*tfac)*exp(-it*tfac);

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

return(nstf);
}

int gen_2tri_stf(float *slip,float *trise,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
int ip, it0, it1, it2;
float tr, amp, a0;
float sum;
float alpha = 0.1;      /* 1st triangle has pulse width = 2*alpha*trise */
float betadeep = 0.2;       /* 2nd triangle has amplitude = beta*A (z0>dmax)*/
float betashal = 0.5;       /* 2nd triangle has amplitude = beta*A (z0<dmin)*/
float beta, dbdd;

float dmin = 4.0;
float dmax = 6.0;

dbdd = (betadeep - betashal)/(dmax-dmin);

if((*z0) >= dmax)
   beta = betadeep;
else if((*z0) < dmax && (*z0) > dmin)
   beta = betadeep - (dmax-(*z0))*dbdd;
else
   beta = betashal;

zapit(stf,nt);

tr = (*trise);
alpha = alpha*tr;

it0 = (int)((alpha)/(*dt) + 0.5);
if(it0 < 2)
   it0 = 2;
it1 = (int)((tr)/(*dt) + 0.5);
if(it1 < 4)
   it1 = 4;

it2 = (2 - beta)*it0;

a0 = 1.0;
amp = a0/(float)(it0);

for(it=0;it<it0;it++)
   stf[it] = it*amp;

for(it=it0;it<it2;it++)
   stf[it] = (2*it0-it)*amp;

amp = beta*a0/(float)(it1-it2);

for(it=it2;it<it1;it++)
   stf[it] = beta*a0 + (it2-it)*amp;

nstf = nt-1;
while(stf[nstf] == (float)(0.0) && nstf)
   nstf--;

if(nstf == 0)
   return(0);

if(nstf < nt-1)
   nstf = nstf + 2;;

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_esg2006_stf(float *slip,float *trise,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
int ip, it0, it1, it2;
float tr, amp, a0, tt;
float sum, arg1;
float alpha, beta, gamma;
float pi = 3.14159265;

zapit(stf,nt);

tr = (*trise);
alpha = 4.0/tr;
beta= 2.0*tr;
gamma = alpha/sqrt(pi);

nstf = (int)((2.0*beta)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

for(it=0;it<nstf;it++)
   {
   tt = it*(*dt);

   arg1 = alpha*(tt - beta);
   stf[it] = exp(-arg1*arg1);
   }

for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int gen_cos_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
float tfac, sum;
float pi = 3.14159265;
float one = 1.0;

zapit(stf,nt);

nstf = (int)((*t0)/(*dt));
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

tfac = 2.0*pi*(*dt)/(*t0);
for(it=0;it<nstf;it++)
   stf[it] = (one - cos(it*tfac));

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}
