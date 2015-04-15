#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

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
   if(ps[ip].slip*ps[ip].slip > MINSLIP*MINSLIP)
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

void write_srf(struct standrupformat *srf,char *file,int bflag)
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

void free_srf_stf(struct standrupformat *srf)
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
}

void read_srfX(struct standrupformat *srf,char *file,int bflag)
{
FILE *fpr, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
char str[1024];

float *stf;
int i, j, k, nt6, it, ip, ig, ntot;

char pword[32];
int fdr;

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
   /*
   sscanf(str,"%s",&(srf[0].version));
   */
   sscanf(str,"%s",srf[0].version);

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

      for(i=0;i<srf[0].srf_apnts.np;i++)
         {
         apval_ptr = &(srf[0].srf_apnts.apntvals[i]);

         fgets(str,1024,fpr);
         sscanf(str,"%f %f %f %f %f %f %f %f",&(apval_ptr->lon),
                                           &(apval_ptr->lat),
                                           &(apval_ptr->dep),
                                           &(apval_ptr->stk),
                                           &(apval_ptr->dip),
                                           &(apval_ptr->area),
                                           &(apval_ptr->tinit),
                                           &(apval_ptr->dt));
         fgets(str,1024,fpr);
         sscanf(str,"%f %f %d %f %d %f %d",&(apval_ptr->rake),
                                        &(apval_ptr->slip1),
                                        &(apval_ptr->nt1),
                                        &(apval_ptr->slip2),
                                        &(apval_ptr->nt2),
                                        &(apval_ptr->slip3),
                                        &(apval_ptr->nt3));

	 if(apval_ptr->nt1)
            apval_ptr->stf1 = (float *)check_malloc((apval_ptr->nt1)*sizeof(float));
	 else
            apval_ptr->stf1 = NULL;

         stf = apval_ptr->stf1;

         for(it=0;it<(apval_ptr->nt1);it++)
            fscanf(fpr,"%f",&stf[it]);

	 if(apval_ptr->nt2)
            apval_ptr->stf2 = (float *)check_malloc((apval_ptr->nt2)*sizeof(float));
	 else
            apval_ptr->stf2 = NULL;

         stf = apval_ptr->stf2;

         for(it=0;it<(apval_ptr->nt2);it++)
            fscanf(fpr,"%f",&stf[it]);

	 if(apval_ptr->nt3)
            apval_ptr->stf3 = (float *)check_malloc((apval_ptr->nt3)*sizeof(float));
	 else
            apval_ptr->stf3 = NULL;

         stf = apval_ptr->stf3;

         for(it=0;it<(apval_ptr->nt3);it++)
            fscanf(fpr,"%f",&stf[it]);

         /* get rouge newline character */
         if((apval_ptr->nt1) || (apval_ptr->nt2) || (apval_ptr->nt3))
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


int gen_ucsb_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
float tau, tau1, tau2, tau1x2, arg1, arg2;
float sum, t, alpha;
float pi = 3.141592654;

zapit(stf,nt);

tau = (*t0);
tau1 = 0.13*tau;
tau2 = tau - tau1;
tau1x2 = 2.0*tau1;


nstf = (int)((tau)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

for(it=0;it<nstf;it++)
   {
   t = it*(*dt);

   alpha = 0.0;
   if(t < tau1)
      {
      arg1 = pi*t/tau1;
      arg2 = 0.5*arg1;
      alpha = 0.7 - 0.7*cos(arg1) + 0.6*sin(arg2);
      }
   else if(t < tau1x2)
      {
      arg1 = pi*t/tau1;
      arg2 = pi*(t - tau1)/tau2;
      alpha = 1.0 - 0.7*cos(arg1) + 0.3*cos(arg2);
      }
   else if(t < tau)
      {
      arg1 = pi*(t - tau1)/tau2;
      alpha = 0.3 + 0.3*cos(arg1);
      }

   stf[it] = alpha;
   }


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

int gen_seki_stf(float *slip,float *t0,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
float tfac, sum, t, htan, ep, em, arg;

zapit(stf,nt);

nstf = (int)(1.5*(*t0)/(*dt));
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

tfac = 2.0/(*t0);
for(it=0;it<nstf;it++)
   {
   t = it*(*dt);
   arg = 2.0*(t*tfac - 1.0);
   ep = exp(arg);
   em = exp(-arg);
   htan = (ep - em)/(ep + em);

   stf[it] = (*slip)*tfac*(1.0 - htan*htan);
   }

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

void sum_srf(struct standrupformat *srf0,struct standrupformat *srf1,struct standrupformat *srf2,float *new_rake)
{
struct srf_prectsegments *prseg_ptr0, *prseg_ptr2;
struct srf_apointvalues *apval_ptr[3];
float *stf, *stf1, *stf2, *stf3;
float tdel, rdif;
double cosR, sinR;
int i, j, k, it, ip, ig, newnt, itdel, ir;

double rperd = 0.017453293;

if(srf0[0].srf_apnts.np != srf1[0].srf_apnts.np)
   {
   fprintf(stderr,"*** number of points in srf1 (%d) not equal to number of points in srf2 (%d), exiting...\n",srf0[0].srf_apnts.np,srf1[0].srf_apnts.np);
   exit(-1);
   }

/* 1st, copy all header info from srf0 to srf2 */

strcpy(srf2[0].version,srf0[0].version);

srf2[0].type[0] = '\0';
if(strncmp(srf0[0].type,"PLANE",5) == 0)
   {
   strcpy(srf2[0].type,srf0[0].type);
   srf2[0].srf_prect.nseg = srf0[0].srf_prect.nseg;

   srf2[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf2[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_ptr0 = srf0[0].srf_prect.prectseg;
   prseg_ptr2 = srf2[0].srf_prect.prectseg;

   for(ig=0;ig<srf2[0].srf_prect.nseg;ig++)
      {
      prseg_ptr2[ig].elon = prseg_ptr0[ig].elon;
      prseg_ptr2[ig].elat = prseg_ptr0[ig].elat;
      prseg_ptr2[ig].nstk = prseg_ptr0[ig].nstk;
      prseg_ptr2[ig].ndip = prseg_ptr0[ig].ndip;
      prseg_ptr2[ig].flen = prseg_ptr0[ig].flen;
      prseg_ptr2[ig].fwid = prseg_ptr0[ig].fwid;
      prseg_ptr2[ig].stk = prseg_ptr0[ig].stk;
      prseg_ptr2[ig].dip = prseg_ptr0[ig].dip;
      prseg_ptr2[ig].dtop = prseg_ptr0[ig].dtop;
      prseg_ptr2[ig].shyp = prseg_ptr0[ig].shyp;
      prseg_ptr2[ig].dhyp = prseg_ptr0[ig].dhyp;
      }
   }

srf2[0].srf_apnts.np = srf0[0].srf_apnts.np;
srf2[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

for(i=0;i<srf2[0].srf_apnts.np;i++)
   {
   apval_ptr[0] = &(srf0[0].srf_apnts.apntvals[i]);
   apval_ptr[1] = &(srf1[0].srf_apnts.apntvals[i]);
   apval_ptr[2] = &(srf2[0].srf_apnts.apntvals[i]);

   /* these should be all the same for srf0 and srf1 */
   apval_ptr[2]->lon = apval_ptr[0]->lon;
   apval_ptr[2]->lat = apval_ptr[0]->lat;
   apval_ptr[2]->dep = apval_ptr[0]->dep;
   apval_ptr[2]->stk = apval_ptr[0]->stk;
   apval_ptr[2]->dip = apval_ptr[0]->dip;
   apval_ptr[2]->area = apval_ptr[0]->area;
   apval_ptr[2]->dt = apval_ptr[0]->dt;
   
   /* find earliest initiation time, reset pointers so ptr0 is earliest and ptr2 is latest */

   if(apval_ptr[1]->tinit < apval_ptr[0]->tinit) /* need to reset pointers */
      {
      apval_ptr[0] = &(srf1[0].srf_apnts.apntvals[i]);
      apval_ptr[1] = &(srf0[0].srf_apnts.apntvals[i]);
      }
   apval_ptr[2]->tinit = apval_ptr[0]->tinit;

   /* resolve all slip to rake_u1=new_rake, rake_u2=new_rake+90 */

   apval_ptr[2]->rake = *new_rake;

   apval_ptr[2]->slip1 = 0.0;
   apval_ptr[2]->nt1 = 0;
   apval_ptr[2]->stf1 = NULL;

   apval_ptr[2]->slip2 = 0.0;
   apval_ptr[2]->nt2 = 0;
   apval_ptr[2]->stf2 = NULL;

   apval_ptr[2]->slip3 = 0.0;
   apval_ptr[2]->nt3 = 0;
   apval_ptr[2]->stf3 = NULL;

   /* loop over u1, u2, u3 for srf0 and srf1 to build up stfs */

   for(ir=0;ir<2;ir++)
      {
      rdif = apval_ptr[ir]->rake - apval_ptr[2]->rake;
      cosR = cos(rperd*rdif);
      sinR = sin(rperd*rdif);

      apval_ptr[2]->slip1 = apval_ptr[2]->slip1 +
                            apval_ptr[ir]->slip1*cosR - apval_ptr[ir]->slip2*sinR;

      apval_ptr[2]->slip2 = apval_ptr[2]->slip2 +
                            apval_ptr[ir]->slip1*sinR + apval_ptr[ir]->slip2*cosR;

      apval_ptr[2]->slip3 = apval_ptr[2]->slip3 + apval_ptr[ir]->slip3;

      tdel = apval_ptr[ir]->tinit - apval_ptr[2]->tinit;
      itdel = (int)(tdel/(apval_ptr[2]->dt) + 0.5);

      stf1 = apval_ptr[2]->stf1;
      stf2 = apval_ptr[2]->stf2;
      stf3 = apval_ptr[2]->stf3;

      if(apval_ptr[ir]->nt1)
	 {
	 newnt = itdel + apval_ptr[ir]->nt1;

	 if(newnt > apval_ptr[2]->nt1)
	    {
            apval_ptr[2]->stf1 = (float *)check_realloc(apval_ptr[2]->stf1,newnt*sizeof(float));
            stf1 = apval_ptr[2]->stf1;

            for(it=(apval_ptr[2]->nt1);it<newnt;it++)
	       stf1[it] = 0.0;

	    apval_ptr[2]->nt1 = newnt;
	    }

	 if(newnt > apval_ptr[2]->nt2)
	    {
            apval_ptr[2]->stf2 = (float *)check_realloc(apval_ptr[2]->stf2,newnt*sizeof(float));
            stf2 = apval_ptr[2]->stf2;

            for(it=(apval_ptr[2]->nt2);it<newnt;it++)
	       stf2[it] = 0.0;

	    apval_ptr[2]->nt2 = newnt;
	    }

         stf = apval_ptr[ir]->stf1;
         for(it=itdel;it<newnt;it++)
	    {
            stf1[it] = stf1[it] + stf[it-itdel]*cosR;
            stf2[it] = stf2[it] + stf[it-itdel]*sinR;
	    }
	 }

      if(apval_ptr[ir]->nt2)
	 {
	 newnt = itdel + apval_ptr[ir]->nt2;

	 if(newnt > apval_ptr[2]->nt1)
	    {
            apval_ptr[2]->stf1 = (float *)check_realloc(apval_ptr[2]->stf1,newnt*sizeof(float));
            stf1 = apval_ptr[2]->stf1;

            for(it=(apval_ptr[2]->nt1);it<newnt;it++)
	       stf1[it] = 0.0;

	    apval_ptr[2]->nt1 = newnt;
	    }

	 if(newnt > apval_ptr[2]->nt2)
	    {
            apval_ptr[2]->stf2 = (float *)check_realloc(apval_ptr[2]->stf2,newnt*sizeof(float));
            stf2 = apval_ptr[2]->stf2;

            for(it=(apval_ptr[2]->nt2);it<newnt;it++)
	       stf2[it] = 0.0;

	    apval_ptr[2]->nt2 = newnt;
	    }

         stf = apval_ptr[ir]->stf2;
         for(it=itdel;it<newnt;it++)
	    {
            stf1[it] = stf1[it] - stf[it-itdel]*sinR;
            stf2[it] = stf2[it] + stf[it-itdel]*cosR;
	    }
	 }

      if(apval_ptr[ir]->nt3)
	 {
	 newnt = itdel + apval_ptr[ir]->nt3;

	 if(newnt > apval_ptr[2]->nt3)
	    {
            apval_ptr[2]->stf3 = (float *)check_realloc(apval_ptr[2]->stf3,newnt*sizeof(float));
            stf3 = apval_ptr[2]->stf3;

            for(it=(apval_ptr[2]->nt3);it<newnt;it++)
	       stf3[it] = 0.0;

	    apval_ptr[2]->nt3 = newnt;
	    }

         stf = apval_ptr[ir]->stf3;
         for(it=itdel;it<newnt;it++)
            stf3[it] = stf3[it] + stf[it-itdel];
	 }
      }

   if(apval_ptr[2]->nt1 == 0)
      apval_ptr[2]->slip1 = 0.0;
   if(apval_ptr[2]->nt2 == 0)
      apval_ptr[2]->slip2 = 0.0;
   if(apval_ptr[2]->nt3 == 0)
      apval_ptr[2]->slip3 = 0.0;
   }
}

void join_srf(struct standrupformat *srf0,struct standrupformat *srf1,struct standrupformat *srf2)
{
struct srf_prectsegments *prseg_in, *prseg_out;
struct srf_apointvalues *apval_in, *apval_out;
float *stfin, *stfout;
int i, j, k, it, ip, ig;

/* 1st, copy all header info from srf0 to srf2 */

strcpy(srf2[0].version,srf0[0].version);

srf2[0].type[0] = '\0';
if(strncmp(srf0[0].type,"PLANE",5) == 0 && strncmp(srf1[0].type,"PLANE",5) == 0)
   {
   strcpy(srf2[0].type,srf0[0].type);

   srf2[0].srf_prect.nseg = srf0[0].srf_prect.nseg + srf1[0].srf_prect.nseg;
   srf2[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf2[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_out = srf2[0].srf_prect.prectseg;
   for(ig=0;ig<srf2[0].srf_prect.nseg;ig++)
      {
      if(ig < srf0[0].srf_prect.nseg)
         {
         k = ig;
         prseg_in = srf0[0].srf_prect.prectseg;
	 }
      else
         {
         k = ig - srf0[0].srf_prect.nseg;
         prseg_in = srf1[0].srf_prect.prectseg;
	 }

      prseg_out[ig].elon = prseg_in[k].elon;
      prseg_out[ig].elat = prseg_in[k].elat;
      prseg_out[ig].nstk = prseg_in[k].nstk;
      prseg_out[ig].ndip = prseg_in[k].ndip;
      prseg_out[ig].flen = prseg_in[k].flen;
      prseg_out[ig].fwid = prseg_in[k].fwid;
      prseg_out[ig].stk = prseg_in[k].stk;
      prseg_out[ig].dip = prseg_in[k].dip;
      prseg_out[ig].dtop = prseg_in[k].dtop;
      prseg_out[ig].shyp = prseg_in[k].shyp;
      prseg_out[ig].dhyp = prseg_in[k].dhyp;
      }
   }

srf2[0].srf_apnts.np = srf0[0].srf_apnts.np + srf1[0].srf_apnts.np;
srf2[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

for(i=0;i<srf2[0].srf_apnts.np;i++)
   {
   if(i < srf0[0].srf_apnts.np)
      {
      k = i;
      apval_in = &(srf0[0].srf_apnts.apntvals[k]);
      }
   else
      {
      k = i - srf0[0].srf_apnts.np;
      apval_in = &(srf1[0].srf_apnts.apntvals[k]);
      }

   apval_out = &(srf2[0].srf_apnts.apntvals[i]);

   apval_out->lon = apval_in->lon;
   apval_out->lat = apval_in->lat;
   apval_out->dep = apval_in->dep;
   apval_out->stk = apval_in->stk;
   apval_out->dip = apval_in->dip;
   apval_out->area = apval_in->area;
   apval_out->tinit = apval_in->tinit;
   apval_out->dt = apval_in->dt;
   apval_out->rake = apval_in->rake;

   apval_out->slip1 = apval_in->slip1;
   apval_out->nt1 = apval_in->nt1;
   apval_out->stf1 = NULL;

   if(apval_out->nt1)
      {
      apval_out->stf1 = (float *)check_realloc(apval_out->stf1,(apval_out->nt1)*sizeof(float));

      stfin = apval_in->stf1;
      stfout = apval_out->stf1;

      for(it=0;it<(apval_out->nt1);it++)
         stfout[it] = stfin[it];
      }

   apval_out->slip2 = apval_in->slip2;
   apval_out->nt2 = apval_in->nt2;
   apval_out->stf2 = NULL;

   if(apval_out->nt2)
      {
      apval_out->stf2 = (float *)check_realloc(apval_out->stf2,(apval_out->nt2)*sizeof(float));

      stfin = apval_in->stf2;
      stfout = apval_out->stf2;

      for(it=0;it<(apval_out->nt2);it++)
         stfout[it] = stfin[it];
      }

   apval_out->slip3 = apval_in->slip3;
   apval_out->nt3 = apval_in->nt3;
   apval_out->stf3 = NULL;

   if(apval_out->nt3)
      {
      apval_out->stf3 = (float *)check_realloc(apval_out->stf3,(apval_out->nt3)*sizeof(float));

      stfin = apval_in->stf3;
      stfout = apval_out->stf3;

      for(it=0;it<(apval_out->nt3);it++)
         stfout[it] = stfin[it];
      }

   }
}
