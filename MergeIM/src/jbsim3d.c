/*
Scott Callaghan - Modified to be called from a wrapper.
Changed all getpar arguments to parameters.
*/

#include "include.h"
#include "structure.h"
#include "function.h"

struct sgtindex *get_sgtpars(struct sgtfileparams *sgtfpar,
			     struct sgtmaster *sgtmast,
			     struct sgtindex *sgtindx);


float* jbsim3d_synth(char stat[], float slon, float slat, char rupmodfile[], char sgt_xfilename[], char sgt_yfilename[], int output_binary, int merge_output, int ntout, char seis_file[])
{
FILE *fopfile(), *fpr;
struct sgtfileparams sgtfilepar, sgtextract;
struct sgtparams *sgtparms;
struct sgtmaster sgtmast;
struct sgtindex *sgtindx;
struct sgtindex eqindx, statindx;
struct geoprojection geop;
float *gfmech;
float *stf, *seis, *subseis, *se, *sn, *sv;
float rt, scale;
float elon, elat, edep;
float vslip, *space;
float z0, strike, dip, rake;
int fdw, ip, maxmech, nstf, ntsum, maxnt, ig;
float mindt, maxdelta, fweight;

struct sgtheader *sgthead;
float *sgtbuf;

char string[256], outfile[128];
char outdir[256], sgtdir[256], sname[8];

struct standrupformat srf;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

struct mechparam mechpar;

long long *indx_master;
int nm, non_exact;

int extract_sgt = 0;

float tmom = 0.0;
float dtout = -1.0;
float slip_conv = 1.0;  /* input slip in cm on each subfault */
float tstart = 0.0;

float sdep = 0.0;

int apv_off = 0;
int inbin = 0;

int intmem = 0;
int memlen;

int ptol;
int print_tol = 25;

//char seis_file[256];

sname[0] = '\0';

sgtfilepar.xfile[0] = '\0';
sgtfilepar.yfile[0] = '\0';
sgtfilepar.zfile[0] = '\0';

sgtfilepar.xfdr = -1;
sgtfilepar.yfdr = -1;
sgtfilepar.zfdr = -1;

sgtextract.xfile[0] = '\0';
sgtextract.yfile[0] = '\0';
sgtextract.zfile[0] = '\0';

sgtextract.xfdr = -1;
sgtextract.yfdr = -1;
sgtextract.zfdr = -1;

sprintf(outdir,".");
sprintf(sgtdir,".");

/*setpar(ac, av);

mstpar("slat","f",&slat);
mstpar("slon","f",&slon);
getpar("outdir","s",outdir);
mstpar("stat","s",stat);

mstpar("rupmodfile","s",rupmodfile);
getpar("slip_conv","f",&slip_conv);
getpar("inbin","d",&inbin);

getpar("sgt_xfile","s",sgtfilepar.xfile);
getpar("sgt_yfile","s",sgtfilepar.yfile);
getpar("sgt_zfile","s",sgtfilepar.zfile);

getpar("extract_sgt","d",&extract_sgt);

getpar("outputBinary","d",&output_binary);*/

strcpy(sgtfilepar.xfile, sgt_xfilename);
strcpy(sgtfilepar.yfile, sgt_yfilename);

if (output_binary==1) { //if output as binary, merge the output files
  merge_output = 1;
}
//getpar("mergeOutput","d",&merge_output);
if (extract_sgt==0) {
  //mstpar("seis_file","s",seis_file);
}

if(sgtfilepar.xfile[0] == '\0' && sgtfilepar.yfile[0] == '\0' && sgtfilepar.zfile[0] == '\0')
   {
   fprintf(stderr,"*** need to specify at least one of sgt_xfile, sgt_yfile, or sgt_zfile; exiting ...\n");
   exit(-1);
   }

/*if(extract_sgt==1)
   {
   getpar("sgtdir","s",sgtdir);

   if(sgtfilepar.xfile[0] != '\0')
      mstpar("extract_sgt_xfile","s",sgtextract.xfile);
   if(sgtfilepar.yfile[0] != '\0')
      mstpar("extract_sgt_yfile","s",sgtextract.yfile);
   if(sgtfilepar.zfile[0] != '\0')
      mstpar("extract_sgt_zfile","s",sgtextract.zfile);
   }*/

/*
getpar("ntout","d",&ntout);
getpar("dtout","f",&dtout);
getpar("tstart","f",&tstart);
getpar("sname","s",sname);*/

//endpar();

printf("slon %f, slat %f, output_binary %d, merge_output %d, ntout %d\n", slon, slat, output_binary, merge_output, ntout);
fflush(stdout);
printf("rupmodfile %s\n", rupmodfile);
printf("inbin %s\n", inbin);

fflush(stdout);

read_srf(&srf,rupmodfile,inbin);
prect_ptr = &srf.srf_prect;
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &srf.srf_apnts;
apval_ptr = apnts_ptr->apntvals;

apv_off = 0;

sgtindx = get_sgtpars(&sgtfilepar,&sgtmast,sgtindx);
set_geoproj(&sgtmast,&geop);

eqindx.h = sgtindx[0].h;
statindx.h = sgtindx[0].h;

get_indx(&slon,&slat,&sdep,&statindx,&geop);

/* find sgt locations (indx) for all fault points */

fprintf(stdout,"Find SGTs for this rupture\n");

sgtparms = (struct sgtparams *) check_malloc ((srf.srf_apnts.np)*sizeof(struct sgtparams));

ptol = print_tol;
maxdelta = 0.0;
fweight = 1.0;
non_exact = 0;
for(ip=0;ip<srf.srf_apnts.np;ip++)
   {
   elon = apval_ptr[ip].lon;
   elat = apval_ptr[ip].lat;
   edep = apval_ptr[ip].dep;
   get_indx(&elon,&elat,&edep,&eqindx,&geop);

   find_sgt(&sgtparms[ip],&sgtmast,sgtindx,&eqindx,&statindx,&maxdelta,&fweight);

   if (sgtparms[ip].nsgt != 1)
	 non_exact++;

   if((float)(100.0*(float)(ip+1)/(float)(srf.srf_apnts.np)) >= ptol)
      {
      fprintf(stdout," %3d percent done (%d of %d)\n",ptol,ip,srf.srf_apnts.np);
      ptol = ptol + print_tol;
      }
   }

/* find unique and sorted master list of locations (indx) */

indx_master = (long long *) check_malloc (4*(srf.srf_apnts.np)*sizeof(long long));
get_master_list(sgtparms,srf.srf_apnts.np,indx_master,&nm);
indx_master = (long long *) check_realloc (indx_master,nm*sizeof(long long));

fprintf(stdout,"nm= %d non_exact= %d\n",nm,non_exact);
fprintf(stdout,"max_mindelta= %f weight= %f\n",sgtindx[0].h*sqrt(maxdelta),fweight);

if(extract_sgt==1)
   {
   fprintf(stdout,"Extracting SGTs for this rupture\n");

   sgt_subset(&sgtfilepar,&sgtextract,&sgtmast,sgtindx,nm,indx_master,sgtdir);
   exit(0);
   }

fprintf(stdout,"Constructing synthetic this rupture\n");

/* try to read all SGTs into memory */

memlen = sizeof(struct sgtmaster) + (sgtmast.globnp)*(sizeof(struct sgtindex) + sizeof(struct sgtheader) + 18*(sgtmast.nt)*sizeof(float));

fprintf(stdout,"Total memory for SGTs= %.2f Mb\n",memlen*1.0e-06);

sgthead = (struct sgtheader *) check_malloc ((sgtmast.globnp)*sizeof(struct sgtheader));
sgtbuf = (float *) check_malloc (18*(sgtmast.globnp)*(sgtmast.nt)*sizeof(float));

read_sgt(&sgtfilepar,&sgtmast,sgtindx,sgthead,sgtbuf);

maxnt = sgthead[0].nt;
mindt = sgthead[0].dt;

if(dtout < 0.0)
   dtout = mindt;

if(dtout < mindt)
   maxnt = (maxnt*mindt/dtout);

ntsum = 2;
while(ntsum < 4*maxnt)
   ntsum = ntsum*2;

if(ntout < 0)
   ntout = ntsum;

maxmech = 3;
mechpar.nmech = 1;
mechpar.flag[0] = U1FLAG;
mechpar.flag[1] = 0;
mechpar.flag[2] = 0;

gfmech = (float *) check_malloc (maxmech*12*ntsum*sizeof(float));
space = (float *) check_malloc (2*ntsum*sizeof(float));

seis = (float *) check_malloc (3*ntout*sizeof(float));
subseis = (float *) check_malloc (maxmech*3*ntout*sizeof(float));
stf = (float *) check_malloc (ntout*sizeof(float));

zapit(seis,3*ntout);

ptol = print_tol;
tmom = 0.0;
for(ip=0;ip<srf.srf_apnts.np;ip++)
   {
   zapit(subseis,maxmech*3*ntout);

   get_srfpars(&srf,apv_off,ip,&rt,&vslip,&mechpar.stk,&mechpar.dip,&mechpar.rak,&mechpar);
   scale = slip_conv*apval_ptr[ip].area;

   mech_sgt(gfmech,sgtbuf,sgthead,&sgtparms[ip],ntsum,mechpar,&scale);
   tmom = tmom + vslip*scale;

   sum_sgt(subseis,ntout,gfmech,&sgtparms[ip],sgthead,ntsum,&rt,&tstart,mechpar);
   srf_stf(&srf,apv_off,ip,seis,subseis,stf,ntout,&dtout,mechpar,space);

   if((float)(100.0*(float)(ip+1)/(float)(srf.srf_apnts.np)) >= ptol)
      {
      fprintf(stdout," %3d percent done (%d of %d)\n",ptol,ip,srf.srf_apnts.np);
      ptol = ptol + print_tol;
      }
   }

float* seis_return = seis;

sv = seis;
sn = seis + ntout;
se = seis + 2*ntout;

if(sname[0] == '\0') {
   strncpy(sname,stat,7);
   sname[7] = '\0';
}

char* last_slash = strrchr(seis_file, '/');
if (last_slash!=NULL) { //means there was a slash in the filename, might need to create a path
	char* path = malloc(sizeof(char)*strlen(seis_file));
	if (path==NULL) {
		printf("Error on malloc, exiting.\n");
		exit(2);
	}
	strncpy(path, seis_file, last_slash-seis_file);
	strcat(path, "\0");
	makedir(path);
}

if (merge_output==0) {
  if (sgtfilepar.xfile[0]!='\0') {
	write_seis(seis_file,sname,"000",sn,&dtout,ntout,&tstart,output_binary);
	seis_return = sn;
  }
  if (sgtfilepar.yfile[0]!='\0') {
	write_seis(seis_file,sname,"090",se,&dtout,ntout,&tstart,output_binary);
  }
  if (sgtfilepar.zfile[0]!='\0') {
	write_seis(seis_file,sname,"ver",sv,&dtout,ntout,&tstart,output_binary);
	seis_return = sv;
  }
} else { //merging output
  if (sgtfilepar.xfile[0]!='\0') {
	if (sgtfilepar.yfile[0]!='\0') {
	  if (sgtfilepar.zfile[0]!='\0') { //x,y,z
		write_seis(seis_file,sname,"grm",sv,&dtout,3*ntout,&tstart,output_binary);
	  } else { //x,y
		write_seis(seis_file,sname,"grm",sn,&dtout,2*ntout,&tstart,output_binary);
		seis_return = sn;
	  }
	}
	else if (sgtfilepar.zfile[0]!='\0') { //x,z
	  printf("Can't output merged X and Z components.\n");
	  exit(2);
	} else { //x
	  	write_seis(seis_file,sname,"grm",sn,&dtout,ntout,&tstart,output_binary);
	}
  } else if (sgtfilepar.yfile[0]!='\0') {
	if (sgtfilepar.zfile[0]!='\0') { //y,z
	 // write_seis(seis_file,sname,"grm",se,&dtout,2*ntout,&tstart,output_binary);
	  printf("Can't output merged Y and Z components.\n");
	} else { //y
	  write_seis(seis_file,sname,"grm",se,&dtout,ntout,&tstart,output_binary);
	  seis_return = se;
	}
  } else if (sgtfilepar.zfile[0]!='\0') { //z
	write_seis(seis_file,sname,"grm",sv,&dtout,ntout,&tstart,output_binary);
  }
}

fprintf(stdout,"Total moment= %13.5e\n",tmom);
printf("returned %ld: sv=%ld, sn=%ld, se=%ld\n", seis_return, sv, sn, se);
return seis_return;
}
