#include "include.h"
#include "structure.h"
#include "function.h"
#include <stdint.h>
#include <mpi.h>

#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         BLOCK_SIZE      1000000
#define         BNDPAD          3

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void get_filepar(char *str,char *file,int *nh,int *latfirst);
void merge_sort(int np, int* my_ixv, int* my_iyv, int* my_izv, long long* my_indx, int start, int end);

main(int ac,char **av)
{
FILE *fopfile(), *fp, *fpr;
float h, *rv, *dv;
float y2, r;
int *rinc, *dinc, nrad, ndep;
int ix, iy, iz, i, id, ir, ifound;
int nx, ny, nz, np, npz;
int xsrc, ysrc;
int *ixv, *iyv, *izv, *izlevel;
double invh, invh2;

int bcnt = 1;
int sflag = 1;
long long *indx, ll_int;

int nhead = 0;
int latfirst = 0;

int ixmin = -999;
int ixmax = -999;
int iymin = -999;
int iymax = -999;
int izstart = 2;
int izmax = -999;

char radiusfile[512], outfile[512], str[1024], nedfile[128];
float xp, yp, north, east, dep, modellon, modellat, modelrot;
float cosR, sinR;

char faultlist[512], infile[1024];
int ixp, iyp, izp;
float flon, flat, fdep, xx, yy, zz;

int geoproj = 1;
int xy2ll = 0;
int ll2xy = 1;
float kmlon, kmlat;

double g0, b0;
double amat[9], ainv[9];

double rperd = RPERD;
float ref_rad = ERAD;

//struct timeval tv;
//struct timeval new_tv;

MPI_Init(&ac, &av);

radiusfile[0] = '\0';
faultlist[0] = '\0';
nedfile[0] = '\0';

setpar(ac, av);

mstpar("nx","d",&nx);
mstpar("ny","d",&ny);
mstpar("nz","d",&nz);
mstpar("h","f",&h);
mstpar("xsrc","d",&xsrc);
mstpar("ysrc","d",&ysrc);

mstpar("modellon","f",&modellon);
mstpar("modellat","f",&modellat);
mstpar("modelrot","f",&modelrot);

mstpar("outfile","s",outfile);

getpar("radiusfile","s",radiusfile);
getpar("faultlist","s",faultlist);
getpar("nedfile","s",nedfile);

getpar("geoproj","d",&geoproj);

getpar("ixmin","d",&ixmin);
getpar("ixmax","d",&ixmax);
getpar("iymin","d",&iymin);
getpar("iymax","d",&iymax);
getpar("izstart","d",&izstart);
getpar("izmax","d",&izmax);

endpar();

int my_id = 0;
int num_procs = 1;
MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

gen_matrices(amat,ainv,&modelrot,&modellon,&modellat);

g0 = (double)(0.5*ny*h)/(double)(ref_rad);
b0 = (double)(0.5*nx*h)/(double)(ref_rad);

if(ixmin < 0)
   ixmin = 0;
if(ixmax < 0 || ixmax > nx)
   ixmax = nx;

if(iymin < 0)
   iymin = 0;
if(iymax < 0 || iymax > ny)
   iymax = ny;

if(izmax < 0 || izmax > nz-izstart)
   izmax = nz - izstart;

invh = 1.0/(double)(h);
invh2 = invh*invh;

ixv = (int *) check_malloc (bcnt*BLOCK_SIZE*sizeof(int));
iyv = (int *) check_malloc (bcnt*BLOCK_SIZE*sizeof(int));
izv = (int *) check_malloc (bcnt*BLOCK_SIZE*sizeof(int));

//Only rank 0 does adaptive mesh so we don't get duplicates
if(my_id==0) {
   //gettimeofday(&tv, NULL);

if(radiusfile[0] != '\0')
   {
   fp = fopfile(radiusfile,"r");

   fscanf(fp,"%d",&nrad);

   rv = (float *) check_malloc (nrad*sizeof(float));
   rinc = (int *) check_malloc (nrad*sizeof(int));

   for(i=0;i<nrad;i++)
      fscanf(fp,"%f",&rv[i]);
   for(i=0;i<nrad;i++)
      fscanf(fp,"%d",&rinc[i]);

   fscanf(fp,"%d",&ndep);

   dv = (float *) check_malloc (ndep*sizeof(float));
   dinc = (int *) check_malloc (ndep*sizeof(int));

   for(i=0;i<ndep;i++)
      fscanf(fp,"%f",&dv[i]);
   for(i=0;i<ndep;i++)
      fscanf(fp,"%d",&dinc[i]);

   fclose(fp);

/* renormalize values to save computations */
   for(i=0;i<nrad;i++)
      rv[i] = rv[i]*rv[i]*invh2;
   for(i=0;i<ndep;i++)
      dv[i] = dv[i]*invh;

   izlevel = (int *) check_malloc (nz*sizeof(int));
   npz = 0;
   for(iz=izstart;iz<izmax+izstart;iz++)
      {
      id = 0;
      while(dv[id] < (float)(iz-1.0) && id != ndep-1)
         id++;

      if((iz-izstart)%dinc[id] == 0)
         {
         izlevel[npz] = iz;
         npz++;
         }
      }

   np = 0;
   for(iy=iymin;iy<iymax;iy++)
      {
      y2 = (ysrc-iy);
      y2 = y2*y2;

      for(ix=ixmin;ix<ixmax;ix++)
         {
         r = (xsrc-ix)*(xsrc-ix) + y2;
         ir = 0;
         while(rv[ir] < r && ir != nrad-1)
            ir++;

         if(ix%rinc[ir] == 0 && iy%rinc[ir] == 0)
            {
	    for(iz=0;iz<npz;iz++)
	       {
	       np++;

               if(np > bcnt*BLOCK_SIZE)
                  {
                  bcnt++;

                  ixv = (int *) check_realloc (ixv,bcnt*BLOCK_SIZE*sizeof(int));
                  iyv = (int *) check_realloc (iyv,bcnt*BLOCK_SIZE*sizeof(int));
                  izv = (int *) check_realloc (izv,bcnt*BLOCK_SIZE*sizeof(int));
                  }

	       ixv[np-1] = ix;
	       iyv[np-1] = iy;
	       izv[np-1] = izlevel[iz];
	       }
	    }
         }
      }
   }

  //gettimeofday(&new_tv, NULL);
  //printf("%f sec for adaptive mesh gen.\n", (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));  
  }

int num_lines = 0;
int num_my_lines;
int avg_lines_per_proc;
char ** local_faultlist_data;
char* tmp;
int i0;
//Create fault files for each process
if (my_id==0) {
	//gettimeofday(&tv, NULL);
	char faultlist_data[10000][200];
	fp = fopfile(faultlist,"r");
	int counter = 0;
	char line[200];
	while(fgets(line, 200, fp)!=NULL) {
		strcpy(faultlist_data[counter], line);
		counter++;
	}
	fclose(fp);
	num_lines = counter;
	if (num_lines > 10000) {
		printf("# of lines is bigger than faultlist_data, exiting.\n");
		exit(1);
	}
	MPI_Bcast(&num_lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
	avg_lines_per_proc = num_lines/num_procs;
	char **data_to_send = check_malloc(sizeof(char*)*num_procs);
	int *lines_per_proc = check_malloc(sizeof(int)*num_procs);
	int i, j, k;
	for (i=0; i<num_procs; i++) {
		lines_per_proc[i] = avg_lines_per_proc;
		if (num_lines - avg_lines_per_proc*num_procs > i) {
			lines_per_proc[i]++;
		}
		data_to_send[i] = check_malloc(sizeof(char)*lines_per_proc[i]*200);
	}
	for (i=0; i<num_procs; i++) {
		k = 0;
		for(j=i; j<num_lines; j+=num_procs) {
			memcpy(data_to_send[i]+200*k, faultlist_data[j], 200);
			k++;
		}			
	}
	for (i=1; i<num_procs; i++) {
		MPI_Send(data_to_send[i], lines_per_proc[i]*200, MPI_CHAR, i, i, MPI_COMM_WORLD);
	}
	num_my_lines = lines_per_proc[0];
        local_faultlist_data = check_malloc(sizeof(char*)*num_my_lines);
        for (i0=0; i0<num_my_lines; i0++) {
                local_faultlist_data[i0] = check_malloc(sizeof(char)*200);
        }
	for (i=0; i<lines_per_proc[0]; i++) {
		strncpy(local_faultlist_data[i], faultlist_data[i], 200);
	}
	free(lines_per_proc);
	for (i=0; i<num_procs; i++) {
		free(data_to_send[i]);
	}
	free(data_to_send);
	//gettimeofday(&new_tv, NULL);
	//printf("%f sec for individual fault file creation.\n", my_id, (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));
} else {
	MPI_Bcast(&num_lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
	avg_lines_per_proc = num_lines/num_procs;
	num_my_lines = avg_lines_per_proc;
	if (num_lines - avg_lines_per_proc*num_procs > my_id) {
        	num_my_lines++;
        }
	local_faultlist_data = check_malloc(sizeof(char*)*num_my_lines);
	tmp = check_malloc(num_my_lines * 200 * sizeof(char));
	for (i0=0; i0<num_my_lines; i0++) {
		local_faultlist_data[i0] = check_malloc(sizeof(char)*200);
	}
	MPI_Status status;
	MPI_Recv(tmp, num_my_lines*200, MPI_CHAR, 0, my_id, MPI_COMM_WORLD, &status);
	int i;
	for (i=0; i<num_my_lines; i++) {
		strncpy(local_faultlist_data[i], tmp+200*i, 200);
	}
	free(tmp);
}

/*  now read in fault coords */
//for(i0=0; i0<num_my_lines; i0++) {
//	printf("Process %d, line %i: %s\n", my_id, i0, local_faultlist_data[i0]);
//}

//change to use buffer

//gettimeofday(&tv, NULL);
if(faultlist[0] != '\0')
   {
    for (i0=0; i0<num_my_lines; i0++) {
	strcpy(str, local_faultlist_data[i0]);
      get_filepar(str,infile,&nhead,&latfirst);

      fprintf(stderr,"%s nheader=%d latfirst=%d\n",infile,nhead,latfirst);
      fpr = fopfile(infile,"r");

      i = 0;
      while(i < nhead)
         {
         fgets(str,1024,fpr);
	 i++;
	 }

      while(fgets(str,1024,fpr) != NULL)
         {
	 if(latfirst)
            sscanf(str,"%f %f %f",&flat,&flon,&fdep);
	 else
            sscanf(str,"%f %f %f",&flon,&flat,&fdep);

         gcproj(&xx,&yy,&flon,&flat,&ref_rad,&g0,&b0,amat,ainv,ll2xy);

	 /*  go to closest point  */
         ixp = (int)((double)(xx)*invh + 0.5);
         iyp = (int)((double)(yy)*invh + 0.5);
         izp = (int)((double)(fdep)*invh + 1.5);

	 if(ixp >= BNDPAD && ixp < nx-BNDPAD && iyp >= BNDPAD && iyp < ny-BNDPAD && izp >= 1 && izp < nz-BNDPAD)
	    {
	    ifound = 0;
            for(ir=0;ir<np;ir++)
	       {
	       if(ixv[ir] == ixp && iyv[ir] == iyp && izv[ir] == izp)
	          {
	          ifound = 1;
	          break;
	          }
	       }

	    if(ifound == 0)
	       {
	       np++;

               if(np > bcnt*BLOCK_SIZE)
                  {
                  bcnt++;

                  ixv = (int *) check_realloc (ixv,bcnt*BLOCK_SIZE*sizeof(int));
                  iyv = (int *) check_realloc (iyv,bcnt*BLOCK_SIZE*sizeof(int));
                  izv = (int *) check_realloc (izv,bcnt*BLOCK_SIZE*sizeof(int));
                  }

	       ixv[np-1] = ixp;
	       iyv[np-1] = iyp;
	       izv[np-1] = izp;
	       }
	    }
         }

      fclose(fpr);
      }

   }

for(i0=0; i0<num_my_lines; i0++) {
        free(local_faultlist_data[i0]);
}
free(local_faultlist_data);


fprintf(stderr,"np= %d\n",np);

//gettimeofday(&new_tv, NULL);
//printf("proc %d took %f sec for reading in files.\n", my_id, (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));

indx = check_malloc(np*sizeof(long long));

for(ir=0;ir<np;ir++)
   indx[ir] = (long long)(ixv[ir])*(long long)(100000000) + (long long)(iyv[ir])*(long long)(10000) + (long long)(izv[ir]);


//printf("Starting local sort.\n");
//fflush(stdout);
//do local sort

//gettimeofday(&tv, NULL);
merge_sort(np, ixv, iyv, izv, indx, 0, np);
//gettimeofday(&new_tv, NULL);
//printf("proc %d took %f sec for local sort.\n", my_id, (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));

//printf("proc %d is done with local sort.\n", my_id);
//fflush(stdout);

//parallel merge on ixv, iyv, izv

//printf("%d has np %d\n", my_id, np);

//if (my_id==0) gettimeofday(&tv, NULL);

int iter=1;
int flag=0;
while (iter<num_procs && !flag) {
	if ((my_id/iter)%2==1) {
		//send np value 
		int partner_id = my_id-iter;
		MPI_Send(&np, 1, MPI_INT, partner_id, 10+iter, MPI_COMM_WORLD);
		MPI_Send(ixv, np, MPI_INT, partner_id, 20+iter, MPI_COMM_WORLD);
		MPI_Send(iyv, np, MPI_INT, partner_id, 30+iter, MPI_COMM_WORLD);
		MPI_Send(izv, np, MPI_INT, partner_id, 40+iter, MPI_COMM_WORLD);
		flag = 1;
	} else if ((my_id/iter)%2==0) {
		MPI_Status stat;
		//receive and merge
		//Make sure you have a receiving partner
		if (my_id+iter>=num_procs) {
			iter *= 2;
			continue;
		}
		int partner_id = my_id+iter;
		int partner_np;
		MPI_Recv(&partner_np, 1, MPI_INT, partner_id, 10+iter, MPI_COMM_WORLD, &stat);
		int* partner_ixv = check_malloc(sizeof(int)*partner_np);
		int* partner_iyv = check_malloc(sizeof(int)*partner_np);
		int* partner_izv = check_malloc(sizeof(int)*partner_np);
		long long* partner_indx = check_malloc(sizeof(long long)*partner_np);
		MPI_Recv(partner_ixv, partner_np, MPI_INT, partner_id, 20+iter, MPI_COMM_WORLD, &stat);
		MPI_Recv(partner_iyv, partner_np, MPI_INT, partner_id, 30+iter, MPI_COMM_WORLD, &stat);
		MPI_Recv(partner_izv, partner_np, MPI_INT, partner_id, 40+iter, MPI_COMM_WORLD, &stat);

		int* combined_ixv = check_malloc(sizeof(int) * (np + partner_np));
		int* combined_iyv = check_malloc(sizeof(int) * (np + partner_np));
		int* combined_izv = check_malloc(sizeof(int) * (np + partner_np));
		long long * combined_indx = check_malloc(sizeof(long long) * (np + partner_np));

		int combined_np = 0;		
	
		for (ir=0;ir<partner_np;ir++)
   			partner_indx[ir] = (long long)(partner_ixv[ir])*(long long)(100000000) + (long long)(partner_iyv[ir])*(long long)(10000) + (long long)(partner_izv[ir]);

		int j0 = 0;
		int i0 = 0;
		int ind = 0;
		int equal = 0;
		//printf("Proc %d here.\n", my_id);
		//printf("my np: %d, partner np: %d\n", np, partner_np);
		//fflush(stdout);
		while (i0<np || j0<partner_np) {
			if (ind>=np+partner_np) {
				printf("ind too large!\n");
			}
			if (i0==np) {
				for (j0; j0<partner_np; j0++) {
					combined_ixv[ind] = partner_ixv[j0];
					combined_iyv[ind] = partner_iyv[j0];
					combined_izv[ind] = partner_izv[j0];
					combined_indx[ind] = partner_indx[j0];
					ind++;
					combined_np++;
				}
			} else if (j0==partner_np) {
				for (i0; i0<np; i0++) {
					combined_ixv[ind] = ixv[i0];
					combined_iyv[ind] = iyv[i0];
					combined_izv[ind] = izv[i0];
					combined_indx[ind] = indx[i0];
					ind++;
					combined_np++;
				}
			} else {
				if (indx[i0]<partner_indx[j0]) {
                                        combined_ixv[ind] = ixv[i0];
                                        combined_iyv[ind] = iyv[i0];
                                        combined_izv[ind] = izv[i0];
                                        combined_indx[ind] = indx[i0];
					i0++;
				} else if (indx[i0]>partner_indx[j0]) {
                                        combined_ixv[ind] = partner_ixv[j0];
                                        combined_iyv[ind] = partner_iyv[j0];
                                        combined_izv[ind] = partner_izv[j0];
                                        combined_indx[ind] = partner_indx[j0];
					j0++;
				} else if (indx[i0]==partner_indx[j0]) {
					equal++;
					combined_ixv[ind] = ixv[i0];
                                        combined_iyv[ind] = iyv[i0];
                                        combined_izv[ind] = izv[i0];
                                        combined_indx[ind] = indx[i0];
                                        i0++;
					j0++;
				}
				ind++;
				combined_np++;
			}
		}
		np = combined_np;

		//printf("%d equal.\n", equal);

		combined_ixv = check_realloc(combined_ixv, sizeof(int)*np);
                combined_iyv = check_realloc(combined_iyv, sizeof(int)*np);
                combined_izv = check_realloc(combined_izv, sizeof(int)*np);
                combined_indx = check_realloc(combined_indx, sizeof(long long)*np);
		
		//printf("Starting frees.\n");
		//printf("partner_ixv %d, partner_iyv %d, partner_izv %d.\n", partner_ixv, partner_iyv, partner_izv);
		//fflush(stdout);

		free(partner_ixv);
		free(partner_iyv);
		free(partner_izv);
		free(partner_indx);
		//fflush(stdout);
		free(ixv);
		free(iyv);
		free(izv);
		free(indx);

		ixv = combined_ixv;
		iyv = combined_iyv;
		izv = combined_izv;
		indx = combined_indx;

		/*free(combined_ixv);
		free(combined_iyv);
		free(combined_izv);
		free(combined_indx);*/
	}
	iter *= 2;
}

if (my_id==0) {

//gettimeofday(&new_tv, NULL);
//printf("Took %f sec for cross-proc merge.\n", (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));

//gettimeofday(&tv, NULL);
if(nedfile[0] != '\0')
   {
   fp = fopfile(nedfile,"w");

   cosR = cos((90.0+modelrot)*rperd);
   sinR = sin((90.0+modelrot)*rperd);

   fprintf(fp,"%d\n",np);
   for(ir=0;ir<np;ir++)
      {
      xp = (xsrc - ixv[ir])*h;
      yp = (ysrc - iyv[ir])*h;

      north = xp*cosR - yp*sinR;
      east  = xp*sinR + yp*cosR;

      fprintf(fp,"%13.8f %13.8f %13.8f\n",north,east,sqrt(xp*xp + yp*yp));
      }

   fprintf(fp,"%d\n",npz);
   for(id=0;id<npz;id++)
      fprintf(fp,"%13.8f\n",(izv[id]-1.0)*h);

   fclose(fp);
   }

//gettimeofday(&new_tv, NULL);
//printf("Took %f sec to calc north/east/dist.\n", (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));
fflush(stdout);

if(geoproj == 0)
   sprintf(str,"localized flat earth on spherical ellipsoid (kmlon=%f, kmlat=%f)",kmlon,kmlat);
else if(geoproj == 1)
   sprintf(str,"great circle arcs");


//gettimeofday(&tv, NULL);
fp = fopfile(outfile,"w");
//char my_outfile[200];
//sprintf(my_outfile, "%s.%d", outfile, my_id);

//fp = fopfile(my_outfile, "w");

fprintf(fp,"# geoproj= %d (%s)\n",geoproj,str);
fprintf(fp,"# modellon= %12.5f modellat= %12.5f modelrot= %6.2f\n",modellon,modellat,modelrot);
fprintf(fp,"# xlen= %12.5f ylen= %12.5f\n",nx*h,ny*h);
fprintf(fp,"#\n");

fprintf(fp,"%d\n",np);
for(ir=0;ir<np;ir++)
   {
   xx = ixv[ir]*h;
   yy = iyv[ir]*h;
   zz = (izv[ir]-1)*h;   /* subtract 1 for free-surface */
   gcproj(&xx,&yy,&flon,&flat,&ref_rad,&g0,&b0,amat,ainv,xy2ll);
   fprintf(fp,"%6d %6d %6d %.12lld %11.5f %11.5f %12.5f\n",ixv[ir],iyv[ir],izv[ir],indx[ir],flon,flat,zz);
   }

fclose(fp);
//gettimeofday(&new_tv, NULL);
//printf("Took %f sec to write output file.\n", (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));
} //close if (my_id==0)

//printf("%d finalizing.\n", my_id);
//fflush(stdout);

MPI_Finalize();
}

void get_filepar(char *str,char *file,int *nh,int *latfirst)
{
int i;
char *sptr;

int lend = 0;

sscanf(str,"%s",file);
*nh = 0;
*latfirst = 0;

sptr = str;

/* 1st skip past 'file' */
while(sptr[0] != ' ' && sptr[0] != '\t')
   {
   if(sptr[0] == '\n' || sptr[0] == '\0')
      return;
   else
      sptr++;
   }

/* now, skip past white space */
while(sptr[0] == ' ' || sptr[0] == '\t')
   {
   if(sptr[0] == '\n' || sptr[0] == '\0')
      return;
   else
      sptr++;
   }

/* now, get 'nheaders' and/or 'latfirst', if they are given */
while(lend == 0)
   {
   i = 0;
   while(sptr[i] != '=')
      {
      if(sptr[i] == '\n' || sptr[i] == '\0')
         return;
      else
         i++;
      }

   if(strncmp(sptr,"nheader",7) == 0)
      sscanf(sptr+i+1,"%d",nh);

   if(strncmp(sptr,"latfirst",8) == 0)
      sscanf(sptr+i+1,"%d",latfirst);

   while(sptr[0] != ' ' && sptr[0] != '\t')
      {
      if(sptr[0] == '\n' || sptr[0] == '\0')
         {
         lend = 1;
         break;
         }
      else
         sptr++;
      }

   while(sptr[0] == ' ' || sptr[0] == '\t')
      {
      if(sptr[0] == '\n' || sptr[0] == '\0')
         {
         lend = 1;
         break;
         }
      else
         sptr++;
      }
   }
}

void merge_sort(int np, int* my_ixv, int* my_iyv, int* my_izv, long long* my_indx, int start, int end) {
	if (start<end-1) {
		int* local_ixv = check_malloc(np*sizeof(int));
		int* local_iyv = check_malloc(np*sizeof(int));
		int* local_izv = check_malloc(np*sizeof(int));
                long long* local_indx = check_malloc(np*sizeof(long long));
                int i;
                for (i=start; i<end; i++) {
			local_ixv[i] = my_ixv[i];
			local_iyv[i] = my_iyv[i];
			local_izv[i] = my_izv[i];
                        local_indx[i] = my_indx[i];
                }
                merge_sort(np, local_ixv, local_iyv, local_izv, local_indx, start, (start+end)/2);
                merge_sort(np, local_ixv, local_iyv, local_izv, local_indx, (start+end)/2, end);
                i = start;
                int j = (start+end)/2;
                int array_ind = i;
                while (array_ind<end) {
                        if (i==(start+end)/2) {
                                for (j; j<end; j++) {
					my_ixv[array_ind] = local_ixv[j];
					my_iyv[array_ind] = local_iyv[j];
					my_izv[array_ind] = local_izv[j];
                                        my_indx[array_ind] = local_indx[j];
                                        array_ind++;
                                }
                                break;
                        } else if (j==end) {
                                for (i; i<(start+end)/2; i++) {
					my_ixv[array_ind] = local_ixv[i];
                                        my_iyv[array_ind] = local_iyv[i];
                                        my_izv[array_ind] = local_izv[i];
                                        my_indx[array_ind] = local_indx[i];
                                        array_ind++;
                                }
                                break;
                        } else {
                                if (local_indx[i]<=local_indx[j]) {
                                        my_ixv[array_ind] = local_ixv[i];
                                        my_iyv[array_ind] = local_iyv[i];
                                        my_izv[array_ind] = local_izv[i];
                                        my_indx[array_ind] = local_indx[i];
                                        i++;
                                } else {
                                        my_ixv[array_ind] = local_ixv[j];
                                        my_iyv[array_ind] = local_iyv[j];
                                        my_izv[array_ind] = local_izv[j];
                                        my_indx[array_ind] = local_indx[j];
                                        j++;
                                }
                                array_ind++;
                        }
                }
                free(local_ixv);
		free(local_iyv);
		free(local_izv);
		free(local_indx);
        }
}

