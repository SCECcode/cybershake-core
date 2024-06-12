#include "include.h"
#include "structure.h"
#include "function.h"
#include <stdint.h>
#include <mpi.h>
#include <stddef.h>

#include "cfuhash.h"

#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         BLOCK_SIZE      5000000
#define         BNDPAD          3

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void get_filepar(char *str,char *file,int *nh,int *latfirst);

int my_id = -1;

struct entry {
	int x;
	int y;
	int z;
	long long index;
};

int compare_entry(const void* e1, const void* e2) {
	struct entry* entry1 = (struct entry*) e1;
	struct entry* entry2 = (struct entry*) e2;
	long long diff = entry1->index - entry2->index;
	if (diff>INT_MAX) {
		return INT_MAX;
	} else if (diff<INT_MIN) {
		return INT_MIN;
	} else {
		return (int)(diff);
	}
}

int main(int ac,char **av)
{
FILE *fp, *fpr;
float h, *rv, *dv;
float y2, r;
int *rinc, *dinc, nrad, ndep;
int ix, iy, iz, i, id, ir, ifound;
int nx, ny, nz, np=0, npz;
int xsrc, ysrc;
int *ixv, *iyv, *izv;
int *izlevel;
double invh, invh2;

int bcnt = 1;
int sflag = 1;
long long *indx, ll_int;

struct entry* points;

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

struct timeval tv;
struct timeval new_tv;


PI_Init(&ac, &av);

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

int num_procs = 1;

PI_Comm_rank(MPI_COMM_WORLD, &my_id);

PI_Comm_size(MPI_COMM_WORLD, &num_procs);

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

printf("%d starting adaptive mesh.\n", my_id);
fflush(stdout);

//Only rank 0 does adaptive mesh so we don't get duplicates
np = 0;

//hash table
cfuhash_table_t *hash = cfuhash_new_with_initial_size(10000000);
char hashkey[20];
//Need to allocate here because cfuhash doesn't make copies of the data, just the key
//Use hashvals as the list of unique points
struct entry* hashvals = check_malloc(sizeof(struct entry)*BLOCK_SIZE);
struct entry tmp_hashval;


if(my_id==0) {
   gettimeofday(&tv, NULL);

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

		  hashvals = check_realloc(hashvals, bcnt*BLOCK_SIZE*sizeof(struct entry));
                  }

		tmp_hashval.x = ix;
		tmp_hashval.y = iy;
		tmp_hashval.z = izlevel[iz];
		//Changed the way this is built to allot more digits per dimension
		tmp_hashval.index = (long long)tmp_hashval.x*(long long)1000000000000 + (long long)tmp_hashval.y*1000000 + (long long)tmp_hashval.z;
	        sprintf(hashkey, "%lld", tmp_hashval.index);

		if (cfuhash_get(hash, hashkey)==NULL) {
                        //Add it
			//np-1 because we already incremented it a few loops ago
                        hashvals[np-1].x = tmp_hashval.x;
			hashvals[np-1].y = tmp_hashval.y;
			hashvals[np-1].z = tmp_hashval.z;
			hashvals[np-1].index = tmp_hashval.index;
                        cfuhash_put(hash, hashkey, &(hashvals[np-1]));
                }
	       }
	    }
         }
      }
   }

  gettimeofday(&new_tv, NULL);
  printf("%f sec for adaptive mesh gen.\n", (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));  
  }

int num_lines = 0;
int num_my_lines;
int avg_lines_per_proc;
char ** local_faultlist_data;
char* tmp;
int i0;

printf("%d creating fault files.\n", my_id);
fflush(stdout);

//Create fault files for each process
int faultlist_chunk = 10000;
int faultlist_size = faultlist_chunk;
if (my_id==0) {
	gettimeofday(&tv, NULL);
	//char faultlist_data[10000][200];
	char** faultlist_data = check_malloc(sizeof(char*) * faultlist_size);
	for (i0=0; i0<faultlist_chunk; i0++) {
		faultlist_data[i0] = check_malloc(sizeof(char) * 200);
	}
	fp = fopfile(faultlist,"r");
	int counter = 0;
	char line[200];
	while(fgets(line, 200, fp)!=NULL) {
		strcpy(faultlist_data[counter], line);
		counter++;
		if (counter==faultlist_size) {
			faultlist_size += faultlist_chunk;
			printf("Expanding faultlist size to %d lines.\n", faultlist_size);
			faultlist_data = check_realloc(faultlist_data, sizeof(char*) * faultlist_size);
			for (i0=faultlist_size-faultlist_chunk; i0<faultlist_size; i0++) {
				faultlist_data[i0] = check_malloc(sizeof(char) * 200);
			}
		}
	}
	fclose(fp);
	num_lines = counter;
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
		strncpy(local_faultlist_data[i], faultlist_data[i*num_procs], 200);
	}
	free(lines_per_proc);
	for (i=0; i<num_procs; i++) {
		free(data_to_send[i]);
	}
	free(data_to_send);
	for (i=0; i<10000; i++) {
		free(faultlist_data[i]);
	}
	free(faultlist_data);
	gettimeofday(&new_tv, NULL);
	printf("%f sec for individual fault file creation.\n", my_id, (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));
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

//change to use buffer

gettimeofday(&tv, NULL);

if(faultlist[0] != '\0')
   {
   char outname[25];
   sprintf(outname, "core_%d.out", my_id);
   FILE* fp_out = fopen(outname, "w");
 
    for (i0=0; i0<num_my_lines; i0++) {
      strcpy(str, local_faultlist_data[i0]);
      get_filepar(str,infile,&nhead,&latfirst);
      printf("%d) Processing file %d of %d.\n", my_id, i0, num_my_lines);
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

	 tmp_hashval.x = ixp;
	 tmp_hashval.y = iyp;
	 tmp_hashval.z = izp;
	 tmp_hashval.index = (long long)ixp*(long long)1000000000000 + (long long)iyp*(long long)1000000 + izp;
	 sprintf(hashkey, "%lld", tmp_hashval.index); 

	 if(ixp >= BNDPAD && ixp < nx-BNDPAD && iyp >= BNDPAD && iyp < ny-BNDPAD && izp >= 1 && izp < nz-BNDPAD)
	    {
		if (cfuhash_get(hash, hashkey)==NULL) {
			//Add it
			hashvals[np].x = tmp_hashval.x;
			hashvals[np].y = tmp_hashval.y;
			hashvals[np].z = tmp_hashval.z;
			hashvals[np].index = tmp_hashval.index;
			cfuhash_put(hash, hashkey, &(hashvals[np]));
			np++;
			if (np > bcnt*BLOCK_SIZE) {
				bcnt++;
				printf("%d) Expanding hashvals.\n", my_id);
				fflush(stdout);
				hashvals = check_realloc(hashvals, sizeof(struct entry)*bcnt*BLOCK_SIZE);
			}
		}
	    }
         }

      fclose(fpr);
      }
	fflush(fp_out);
	fclose(fp_out);
   }

/*
char outname[25];
sprintf(outname, "core_%d.out", my_id);
FILE* fp_out = fopen(outname, "w");
for (i0=0; i<np; i++) {
	fprintf(fp_out, "%d %d %d %ld\n", hashvals[i].x, hashvals[i].y, hashvals[i].z, hashvals[i].index);
}
fflush(fp_out);
fclose(fp_out);
*/

//check # of entries to be sure it matches
size_t num_fault_points = cfuhash_num_entries(hash);

if (num_fault_points!=np) {
	fprintf(stderr, "Error: we think we inserted %d keys but we actually have %d.  Aborting.\n", np, num_fault_points);
	exit(2);
}

for(i0=0; i0<num_my_lines; i0++) {
        free(local_faultlist_data[i0]);
}
free(local_faultlist_data);


fprintf(stderr,"id=%d, np= %d\n",my_id, np);

gettimeofday(&new_tv, NULL);
printf("proc %d took %f sec for reading in files.\n", my_id, (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));

printf("%d has np %d\n", my_id, np);

fflush(stdout);

//parallel merge on ixv, iyv, izv

if (my_id==0) gettimeofday(&tv, NULL);

//Get total NP
int global_np = 0;
int* proc_nps = check_malloc(sizeof(int)*num_procs);
int* proc_disps = check_malloc(sizeof(int)*num_procs);
int rc = MPI_Gather(&np, 1, MPI_INT, proc_nps, 1, MPI_INT, 0, MPI_COMM_WORLD);
if (rc!=MPI_SUCCESS) {
	fprintf(stderr, "Error in MPI_Gather, aborting.\n");
	exit(-3);
}
if (my_id==0) {
	for (ir=0; ir<num_procs; ir++) {
		printf("Proc %d has np=%d\n", ir, proc_nps[ir]);
		proc_disps[ir] = global_np;
		global_np += proc_nps[ir];
		printf("Proc %d has offset=%d\n", ir, proc_disps[ir]);
	}
}

struct entry* recv_points = NULL;
if (my_id==0) {
	printf("Total of %d gathered elements.", global_np);
	fflush(stdout);
	printf("Allocating %ld bytes.\n", sizeof(struct entry)*(size_t)global_np);
	recv_points = check_malloc(sizeof(struct entry)*(size_t)global_np);
	//Gather points, put in hashmap, read back out
}

PI_Datatype entry_type;
int blen[5] = {1, 1, 1, 1, 1};
//Can't sent an 8-byte LONG via mpi, so use 2 ints
int x_offset, y_offset, z_offset, index_offset;
x_offset = offsetof(struct entry, x);
y_offset = offsetof(struct entry, y);
z_offset = offsetof(struct entry, z);
index_offset = offsetof(struct entry, index);

PI_Aint disps[5] = {x_offset, y_offset, z_offset, index_offset, index_offset+sizeof(int)};

PI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};

PI_Type_create_struct(5, blen, disps, types, &entry_type);

PI_Type_commit(&entry_type);

rc = MPI_Gatherv(hashvals, np, entry_type, recv_points, proc_nps, proc_disps, entry_type, 0, MPI_COMM_WORLD);
if (rc!=MPI_SUCCESS) {
	fprintf(stderr, "Error in MPI_Gatherv, aborting.\n");
	exit(-4);
}

//Add points in recv_points, starting with np, to hash table
if (my_id==0) {
	hashvals = check_realloc(hashvals, sizeof(struct entry)*global_np);
	for(ir=np; ir<global_np; ir++) {
		tmp_hashval.x = recv_points[ir].x;
		tmp_hashval.y = recv_points[ir].y;
		tmp_hashval.z = recv_points[ir].z;
		tmp_hashval.index = recv_points[ir].index;

		sprintf(hashkey, "%lld", tmp_hashval.index);
		if (cfuhash_get(hash, hashkey)==NULL) {
	        	//Add it
	                hashvals[np].x = tmp_hashval.x;
	                hashvals[np].y = tmp_hashval.y;
	                hashvals[np].z = tmp_hashval.z;
	                hashvals[np].index = tmp_hashval.index;
	                cfuhash_put(hash, hashkey, &(hashvals[np]));
			//printf("adding to hashvals: (%d, %d, %d, %ld)\n", hashvals[np].x, hashvals[np].y, hashvals[np].z, hashvals[np].index);
			np++;
	        }
	}
	hashvals = check_realloc(hashvals, sizeof(struct entry)*np);
	/*for(ir=0; ir<np; ir++) {
		printf("hashvals: (%d, %d, %d, %ld)\n", hashvals[ir].x, hashvals[ir].y, hashvals[ir].z, hashvals[ir].index);
	}*/
	/*keys = (char**)cfuhash_keys(hash, &num_fault_points, 0);
	if (num_fault_points!=np) {
		fprintf(stderr, "Error: expecting %d points but actually %d.\n", np, num_fault_points);
		exit(-1);
	}
	points = check_realloc(points, num_fault_points*sizeof(struct entry));
	for(ir=0; ir<np; ir++) {
	        struct entry* val = (struct entry*)cfuhash_get(hash, keys[ir]);
		printf("points: (%d, %d, %d)\n", val->x, val->y, val->z);
	        points[ir].index = val->index;
	        //populate ixv, iyv, izv
	        points[ir].x = val->x;
	        points[ir].y = val->y;
	        points[ir].z = val->z;
	}*/
	
	qsort(hashvals, np, sizeof(struct entry), compare_entry);
}

if (my_id==0) {

gettimeofday(&new_tv, NULL);
printf("Took %f sec for cross-proc merge.\n", (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));

gettimeofday(&tv, NULL);
if(nedfile[0] != '\0')
   {
   fp = fopfile(nedfile,"w");

   cosR = cos((90.0+modelrot)*rperd);
   sinR = sin((90.0+modelrot)*rperd);

   fprintf(fp,"%d\n",np);
   for(ir=0;ir<np;ir++)
      {
      xp = (xsrc - hashvals[ir].x*h);
      yp = (ysrc - hashvals[ir].y*h);

      north = xp*cosR - yp*sinR;
      east  = xp*sinR + yp*cosR;

      fprintf(fp,"%13.8f %13.8f %13.8f\n",north,east,sqrt(xp*xp + yp*yp));
      }

   fprintf(fp,"%d\n",npz);
   for(id=0;id<npz;id++) {
      fprintf(fp,"%13.8f\n",(hashvals[id].z-1.0)*h);
   }

   fclose(fp);
   }

gettimeofday(&new_tv, NULL);
printf("Took %f sec to calc north/east/dist.\n", (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));
fflush(stdout);

if(geoproj == 0)
   sprintf(str,"localized flat earth on spherical ellipsoid (kmlon=%f, kmlat=%f)",kmlon,kmlat);
else if(geoproj == 1)
   sprintf(str,"great circle arcs");


gettimeofday(&tv, NULL);
fp = fopfile(outfile,"w");

fprintf(fp,"# geoproj= %d (%s)\n",geoproj,str);
fprintf(fp,"# modellon= %12.5f modellat= %12.5f modelrot= %6.2f\n",modellon,modellat,modelrot);
fprintf(fp,"# xlen= %12.5f ylen= %12.5f\n",nx*h,ny*h);
fprintf(fp,"#\n");

fprintf(fp,"%d\n",np);
for(ir=0;ir<np;ir++)
   {
   xx = hashvals[ir].x*h;
   yy = hashvals[ir].y*h;
   zz = (hashvals[ir].z-1)*h;  /* subtract 1 for free-surface */
   gcproj(&xx,&yy,&flon,&flat,&ref_rad,&g0,&b0,amat,ainv,xy2ll);
   fprintf(fp,"%6d %6d %6d %.12lld %11.5f %11.5f %12.5f\n",hashvals[ir].x, hashvals[ir].y, hashvals[ir].z, hashvals[ir].index, flon, flat, zz);
   }

fclose(fp);
gettimeofday(&new_tv, NULL);
printf("Took %f sec to write output file.\n", (new_tv.tv_sec - tv.tv_sec + (new_tv.tv_usec - tv.tv_usec)/1000000.0));
} //close if (my_id==0)

free(recv_points);

cfuhash_destroy(hash);
free(proc_nps);
free(proc_disps);
free(hashvals);


PI_Finalize();
return 0;
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

