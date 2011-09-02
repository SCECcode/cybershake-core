#include "include.h"
#include "structure.h"
#include "function.h"

#define         MAXLINE       256
#define         MAXMEM       1500   /* default max RAM in Mbytes */

main(int ac, char **av)
{
int myid, nproc, pnlen;
char procname[128];

FILE *fpr, *fopfile();
struct sgtmaster sgtmast;
struct sgtindex *sgtindx;
struct sgtheader114 *sgthead;
float *gfbuf, *gfptr;
int ig;
int i, nt, j, k, my_k, glob_k, nfiles, it, ip, gflen;
off_t off, head1_off, head2_off;
int localnp, local_cnt, fdr, fdw, nr;
char infile[512];
char stat[16];
char str[512];

int zero_outfile = 0;

float scale = 1.0;

int indx, indxst, indxnd, indx_inc;
int sflag = 1;
int xindx, yindx, zindx;
long long indx2get;

setpar(ac, av);
mstpar("infile","s",infile);
endpar();

fprintf(stderr,"processing file= %s\n",infile);
fflush(stderr);

fdr = opfile_ro(infile);

reed(fdr,&sgtmast,sizeof(struct sgtmaster));

fprintf(stderr,"   geoproj=  %12d\n",sgtmast.geoproj);
fprintf(stderr,"   modellon= %12.5f\n",sgtmast.modellon);
fprintf(stderr,"   modellat= %12.5f\n",sgtmast.modellat);
fprintf(stderr,"   modelrot= %12.5f\n",sgtmast.modelrot);
fprintf(stderr,"   xshift=   %12.5f\n",sgtmast.xshift);
fprintf(stderr,"   yshift=   %12.5f\n",sgtmast.yshift);
fprintf(stderr,"   globnp=   %12d\n",sgtmast.globnp);
fprintf(stderr,"   localnp=  %12d\n",sgtmast.localnp);
fprintf(stderr,"   nt=       %12d\n",sgtmast.nt);
fflush(stderr);

sgtindx = (struct sgtindex *) check_malloc ((sgtmast.globnp)*sizeof(struct sgtindex));
reed(fdr,sgtindx,(sgtmast.globnp)*sizeof(struct sgtindex));

/*
for(ip=0;ip<sgtmast.globnp;ip++)
   fprintf(stderr,"indx[%d]= %.12Ld\n",ip,sgtindx[ip].indx);
fflush(stderr);
*/

sgthead = (struct sgtheader114 *) check_malloc ((sgtmast.localnp)*sizeof(struct sgtheader114));
reed(fdr,sgthead,(sgtmast.localnp)*sizeof(struct sgtheader114));

for(ip=0;ip<sgtmast.localnp;ip++)
   fprintf(stderr,"header_indx[%d]= %.12Ld\n",ip,sgthead[ip].indx);
fflush(stderr);

gfbuf = (float *) check_malloc (6*(sgtmast.nt)*(sgtmast.localnp)*sizeof(float));
reed(fdr,gfbuf,6*(sgtmast.nt)*(sgtmast.localnp)*sizeof(float));

for(ip=0;ip<sgtmast.localnp;ip++)
   {
   fprintf(stderr,"ip=%8d\n",ip);
   for(ig=0;ig<6;ig++)
      {
      for(it=0;it<sgtmast.nt;it++)
         {
	 k = it*6*sgtmast.localnp + 6*ip + ig;

         if(isnan(gfbuf[k]))
            fprintf(stderr,"nan found: ig=%d it=%5d lam= %13.5e\n",ig,it,sgthead[ip].lam);
         }
      }
   }

/*
   fprintf(stderr,"   *** reading headers ... ");
   fflush(stderr);


   fprintf(stderr,"done\n");
   fflush(stderr);

      fprintf(stderr,"sgthead[0].indx= %.12Ld\n",sgthead[0].indx);
      fprintf(stderr,"sgthead[1].indx= %.12Ld\n",sgthead[1].indx);
      fflush(stderr);

   ip = 0;
   while(indx2get > sgthead[ip].indx && ip < (sgtmast.localnp)-1)
      ip++;

      fprintf(stderr,"ip= %d indx2get= %.12Ld sgthead[ip].indx= %.12Ld\n",ip,indx2get,sgthead[ip].indx);
      fprintf(stderr,"ip= %d indx2get= %.12Ld sgtindx[ip].indx= %.12Ld\n",ip,indx2get,sgtindx[ip].indx);
      fflush(stderr);

   if(indx2get == sgthead[ip].indx)
      {
      fprintf(stderr,"found indx= %.12Ld, getting SGT ...",indx2get);
      fflush(stderr);

      head1_off = (off_t)(sizeof(struct sgtmaster)) + (off_t)((sgtmast.globnp)*sizeof(struct sgtindex)) + (off_t)((sgtmast.localnp)*sizeof(struct sgtheader114));

      off = head1_off + ip*6*(sgtmast.nt)*sizeof(float);
      lseek(fdr,off,SEEK_SET);

      gfbuf = (float *) check_malloc (6*(sgtmast.nt)*sizeof(float));
      reed(fdr,gfbuf,6*(sgtmast.nt)*sizeof(float));
      close(fdr);

      sgtmast.localnp = 1;

      gfptr = gfbuf;
      write_seis("./",stat,stat,"gxx",gfptr,&sgthead[ip].dt,sgthead[ip].nt,&sgthead[ip].tst);
      gfptr = gfbuf + 1*(sgtmast.nt);
      write_seis("./",stat,stat,"gyy",gfptr,&sgthead[ip].dt,sgthead[ip].nt,&sgthead[ip].tst);
      gfptr = gfbuf + 2*(sgtmast.nt);
      write_seis("./",stat,stat,"gzz",gfptr,&sgthead[ip].dt,sgthead[ip].nt,&sgthead[ip].tst);
      gfptr = gfbuf + 3*(sgtmast.nt);
      write_seis("./",stat,stat,"gxy",gfptr,&sgthead[ip].dt,sgthead[ip].nt,&sgthead[ip].tst);
      gfptr = gfbuf + 4*(sgtmast.nt);
      write_seis("./",stat,stat,"gxz",gfptr,&sgthead[ip].dt,sgthead[ip].nt,&sgthead[ip].tst);
      gfptr = gfbuf + 5*(sgtmast.nt);
      write_seis("./",stat,stat,"gyz",gfptr,&sgthead[ip].dt,sgthead[ip].nt,&sgthead[ip].tst);

      fprintf(stderr,"done\n");
      fflush(stderr);
      }
*/

close(fdr);
exit(0);
}

void write_seis(char *dir,char *stat,char *sname,char *comp,float *st,float *dt,int nt,float *ts)
{
FILE *fopfile(), *fpw;
int i, j, nt6;
char outfile[256], header[128], stitle[128];

sprintf(outfile,"%s/%s.%s",dir,stat,comp);

fpw = fopfile(outfile,"w");

sprintf(stitle,"%-10s%3s %s\n",sname,comp,"TITLE");
fprintf(fpw,"%s",stitle);

sprintf(header,"%d %12.5e 0 0 %12.5e 0.0 0.0 0.0\n",nt,*dt,*ts);
fprintf(fpw,"%s",header);

nt6 = nt/6;
for(i=0;i<nt6;i++)
   {
   for(j=0;j<6;j++)
      fprintf(fpw,"%13.5e",st[6*i + j]);

   fprintf(fpw,"\n");
   }
if(6*nt6 != nt)
   {
   for(i=6*nt6;i<nt;i++)
      fprintf(fpw,"%13.5e",st[i]);

   fprintf(fpw,"\n");
   }
fclose(fpw);
}
