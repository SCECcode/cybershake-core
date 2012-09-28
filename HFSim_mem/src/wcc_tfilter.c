/**********************************************************************/
/*                                                                    */
/*           wcc_tfilter - written by RWG 03/92                       */
/*                                                                    */
/*           N-th order high-pass, low-pass and band-pass             */
/*           Butterworth filters for time series.                     */
/*                                                                    */
/**********************************************************************/

#include "include.h"
#include "structure.h"
#include "function.h"

#define 	PI 		3.14159265
#define         MAXFILES        500

char infilebuf[MAXFILES*128];
char outfilebuf[MAXFILES*128];
char stringbuf[256];

int wcc_tfilter (float** seis, int order, float fhi, float flo, int phase, int ncomps, int nt_in, float dt_in) {
struct statdata shead1;
char *infile[MAXFILES], *outfile[MAXFILES], *string, *readline(), filelist[128];
char str[512];
struct complex *q, *p, *tmpptr;
float *s1;
double trig_arg;
int j, it, i;
int nshft, inchar, outchar, nstat, nt6;
int rval;
float wplo, wphi, cosA, sinA;
float fnyq;
float are, aim;
float one = 1.0;
float two = 2.0;

int inbin = 0;
int outbin = 0;

FILE *fopfile(), *fmake_or_open(), *fpr, *fpw;
char outpath[256];

/*if(ac == 1)
   {
   printf("**** %s:\n",av[0]);
   printf("     'flo' is low pass frequency cutoff (pass for f < flo)\n");
   printf("     'fhi' is high pass frequency cutoff (pass for f > fhi)\n\n");
   printf("     This means that the passband will be:  'fhi' < f < 'flo'.\n\n");
   printf("     For a high-pass only filter, set fhi=f0, flo=1.0e+10;\n");
   printf("     for a low-pass only filter, set fhi=0.0, flo=f0;\n");
   printf("     where 'f0' is the cutoff frequency in either case.\n");
   exit(1);
   }


setpar(ac, av);
mstpar("filelist","s",filelist);
mstpar("outpath","s",outpath);
getpar("order","d",&order);
getpar("fhi","f",&fhi);
getpar("flo","f",&flo);
getpar("phase","d",&phase);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
endpar();*/

/*  read in input data filenames */
 
//fpr = fopfile(filelist,"r");

/*i = 0;                       
infile[0] = infilebuf;
outfile[0] = outfilebuf;
string = readline(fpr);
while(string != NULL)
   {
   nshft = 0;
   inchar = getname(string,infile[i],&nshft);
   outchar = getname(&string[nshft],outfile[i],&nshft);

   if(outfile[i][0] == '\0')
      {
      j = 0;                                              
      while(infile[i][j] != '\0')
         j++;
 
      outchar = 0;
      while(j >= 0 && infile[i][j] != '/')
	 {
         j--;
	 outchar++;
	 }
      j++;  
 
      strcpy(outfile[i],infile[i]+j);
      }

   i++;
   if(i == MAXFILES)
      break;
   else
      {
      infile[i] = infile[i-1] + inchar;
      outfile[i] = outfile[i-1] + outchar;
      }

   string = readline(fpr);
   }
nstat = i;
fclose(fpr);
fflush(stdout);
if(!nstat) {
  printf("No entries in file %s.\n");
  fflush(stdout);
  exit(0);
}*/

s1 = NULL;
p = NULL;
q = NULL;
shead1.nt = nt_in;
shead1.dt = dt_in;

int size_cx = sizeof(struct complex);
for(i=0;i<ncomps;i++)
   {
   //s1 = read_wccseis(infile[i],&shead1,s1,inbin);
   s1 = seis[i];

   p = (struct complex *) check_realloc(p,shead1.nt*size_cx);
   q = (struct complex *) check_realloc(q,shead1.nt*size_cx);

   for(it=0;it<shead1.nt;it++)
      {
      p[it].re = s1[it];
      p[it].im = 0.0;
      }

   czero(q,shead1.nt);
 
   fnyq = one/(two*shead1.dt);

   if(fhi < fnyq)
      wphi = tan(PI*fhi*shead1.dt);
   else
      wphi = 1.0e+10;

   if(flo < fnyq)
      wplo = tan(PI*flo*shead1.dt);

	   /*  forward pass  */

   for(j=0;j<order;j++)
      {
      trig_arg = 0.5*(two*j + one)*PI/order;
      cosA = cos(trig_arg);
      sinA = -sin(trig_arg);

      if(flo < fnyq)    /* low-pass filter */
         {
         are = wplo*sinA;
         aim = -wplo*cosA;
         lp_filter(q,p,shead1.nt,&are,&aim,1);
      
         tmpptr = p; p = q; q = tmpptr;
         }
   
      if(fhi >= 0.0)    /* high-pass filter */
         {
         are = wphi*sinA;
         aim = -wphi*cosA;
         hp_filter(q,p,shead1.nt,&are,&aim,1);
      
         tmpptr = p; p = q; q = tmpptr;
         }
      }
   
	      /*  reverse pass to obtain zero-phase response  */
   
   if(phase == 0)
      {
      for(j=0;j<order;j++)
         {
         trig_arg = 0.5*(two*j + one)*PI/order;
         cosA = cos(trig_arg);
         sinA = -sin(trig_arg);
   
         if(flo < fnyq)    /* low-pass filter */
            {
            are = wplo*sinA;
            aim = -wplo*cosA;
            lp_filter(q,p,shead1.nt,&are,&aim,-1);
         
            tmpptr = p; p = q; q = tmpptr;
            }
   
         if(fhi >= 0.0)    /* high-pass filter */
            {
            are = wphi*sinA;
            aim = -wphi*cosA;
            hp_filter(q,p,shead1.nt,&are,&aim,-1);
   
            tmpptr = p; p = q; q = tmpptr;
            }
         }
      }

   /* make sure output directory exists */
 
   //makedir(outpath);

   for(it=0;it<shead1.nt;it++)
      s1[it] = p[it].re;

   //set_fullpath(str,outpath,outfile[i]);
   //printf("Writing to %s.\n", str);
   //write_wccseis(str,&shead1,s1,outbin);
   free(p);
   free(q);
   p = NULL;
   q = NULL;
   }
   return(0);
}

czero(struct complex* p, int n) {
float zap = 0.0;
int i;

for(i=0;i<n;i++)
   p[i].re = p[i].im = zap;
}

hp_filter(q,p,n,alpha,beta,sgn)
struct complex *q, *p;
float *alpha, *beta;
int n, sgn;
{
float are, aim, bre, bim;
float tmpre, tmpim, denom;
float one = 1.0;
int i, n1, k;

tmpre = one - (*alpha);
tmpim = -(*beta);
denom = one/(tmpre*tmpre + tmpim*tmpim);
are = tmpre*denom;
aim = -tmpim*denom;

bre = one + (*alpha);
bim = (*beta);

n1 = n - 1;
if(sgn == 1)
   i = 1;
if(sgn == -1)
   i = n1 - 1;

k = i - sgn;

q[k].re = (are*p[k].re - aim*p[k].im);
q[k].im = (are*p[k].im + aim*p[k].re);

while(n1--)
   {
   tmpre = (p[i].re - p[k].re) + bre*q[k].re - bim*q[k].im;
   tmpim = (p[i].im - p[k].im) + bre*q[k].im + bim*q[k].re;

   q[i].re = are*tmpre - aim*tmpim;
   q[i].im = are*tmpim + aim*tmpre;

   i = i + sgn;
   k = i - sgn;
   }
}

lp_filter(q,p,n,alpha,beta,sgn)
struct complex *q, *p;
float *alpha, *beta;
int n, sgn;
{
float are, aim, bre, bim, cre, cim;
float tmpre, tmpim, denom;
float one = 1.0;
int i, n1, k;

tmpre = one - (*alpha);
tmpim = -(*beta);
denom = one/(tmpre*tmpre + tmpim*tmpim);
are = tmpre*denom;
aim = -tmpim*denom;

bre = one + (*alpha);
bim = (*beta);

cre = -(*alpha);
cim = -(*beta);

n1 = n - 1;
if(sgn == 1)
   i = 1;
if(sgn == -1)
   i = n1 - 1;

k = i - sgn;

tmpre = cre*p[k].re - cim*p[k].im;
tmpim = cre*p[k].im + cim*p[k].re;

q[k].re = are*tmpre - aim*tmpim;
q[k].im = are*tmpim + aim*tmpre;

while(n1--)
   {
   tmpre = cre*(p[i].re + p[k].re) - cim*(p[i].im + p[k].im)
      + bre*q[k].re - bim*q[k].im;
   tmpim = cre*(p[i].im + p[k].im) + cim*(p[i].re + p[k].re)
      + bre*q[k].im + bim*q[k].re;

   q[i].re = are*tmpre - aim*tmpim;
   q[i].im = are*tmpim + aim*tmpre;

   i = i + sgn;
   k = i - sgn;
   }
}

set_fullpath(fname,path,name)
char *path, *name, *fname;
{
int j;

j = 0;
while(path[j] != '\0')
   j++;
if(path[j-1] != '/')
   {  
   path[j] = '/';
   path[j+1] = '\0';
   }  

sprintf(fname,"%s%s",path,name);
}

makedir(path)
char *path;
{
char stmp[256], str[128];
int rtn, j; 
mode_t mode = 00777;

j = 0;   
while(path[j] != '\0')
   j++;  
 
j--;
while(path[j] == '/')
   j--;
path[j+1] = '\0';
path[j+2] = '\0';
         
j = 0;
while(path[j] != '\0')
   {
   while(path[j] != '/' && path[j] != '\0')
      j++;  
 
   if(j != 0)
      {
      strncpy(stmp,path,j);
      stmp[j] = '\0';
      rtn = mkdir(stmp,mode);
 
      if(rtn == -1)
         {
         if(errno != EEXIST)
            {
            sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
            perror(str);
            exit(-1);
            }
         }
      }
   j++;
   }
}

char *readline(fp)
FILE *fp;
{
char *ptr, c;

ptr = stringbuf;
while((c = getc(fp)) != EOF)
   {
   if(c == '\n')
      {
      *ptr = '\0';
      return(stringbuf);
      }
   *ptr++ = c;
   }
return(NULL);
}

getname(str,name,nshft)
char *str, *name;
int *nshft;
{
int nc, nc1;

nc1 = 0;
while(*str == ' ' || *str == '\t') /* remove blank space */
   {
   str++;
   nc1++;
   }
 
if(*str == '"')
   {
   str++;
   nc1++;
   nc = 0;
   while(str[nc] != '"')
      {
      if(str[nc] == '\0')
         break;
 
      nc++;
      }
 
   if(str[nc] == '"')
      nc1++;
   }
else  
   {
   nc = 0;
   while(str[nc] != '\0')
      {
      if(str[nc] == ' ' || str[nc] == '\t')
         break;
 
      nc++;
      }
   }

strncpy(name,str,nc);
name[nc] = '\0';
*nshft = *nshft + nc1 + nc;
nc++;
return(nc);
}
