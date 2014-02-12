#include "mpi.h"
#include "include.h"
#include "function.h"
#include "func_mpi.h"
#include "ucvm.h"


/* Vp/Vs ratio */
#define MIN_V_RATIO 1.0
#define VX_NO_DATA -99998.0


int main(int ac,char **av)
{
/* MPI-RWG: MPI stuff and distributed computation variables */
int myid, nproc, pnlen;
char procname[128];

 FILE *fpr, *ppr;
float *fbuf;
float *mlon, *mlat, *mdep;

int fdw, i, j, k, nn, newn;
int nx, ny, nz, ix, iz, icnt;
int iy0, iy1, iy, local_ny;

char infile[512], outdir[512], fileroot[512], cordfile[512], depfile[512], outfile[512], modeldir[512], models[512];
 char logdir[512], logname[512], logfile[512], pointfile[512];
char str[512];

float meter2feet = 3.2808399;
float meter2km = 0.001;
int tomobg = 1;
int geot_layer = 1;
int max_iter = 20;

mpi_init(&ac,&av,&nproc,&myid,procname,&pnlen);

sprintf(logdir,".");
sprintf(logname,"vh_mpi");

setpar(ac,av);
mstpar("nx","d",&nx);
mstpar("ny","d",&ny);
mstpar("nz","d",&nz);
mstpar("cordfile","s",cordfile);
mstpar("depfile","s",depfile);
mstpar("modeldir","s",modeldir);

mstpar("outdir","s",outdir);
mstpar("fileroot","s",fileroot);

getpar("logdir","s",logdir);
getpar("logname","s",logname);

getpar("tomobg","d",&tomobg);
getpar("geot_layer","d",&geot_layer);
getpar("models","s",models);

endpar();

nn = nx*nz;

makedir(outdir);
makedir(logdir);

sprintf(logfile,"%s/%s-%.4d.log",logdir,logname,myid);
sprintf(pointfile,"%s/%s-%.4d.log",logdir,"points",myid);
freopen(logfile,"w",stderr);
ppr = fopen(pointfile, "w");

fprintf(stderr,"***** This is UCVM-mpi\n");
fprintf(stderr,"      Running on node= %s rank= %d\n\n",procname,myid);
fflush(stderr);

/*
    Find iy slices for this processor.
    Responsible for planes iy=iy0,...,iy1-1
    local_ny = iy1-iy0
*/

get_iylimits(nproc,myid,&iy0,&iy1,ny);
local_ny = iy1 - iy0;

mdep = (float *)check_malloc(nz*sizeof(float));

fpr = fopfile(depfile,"r");

for(iz=0;iz<nz;iz++)
   fscanf(fpr,"%f",&mdep[iz]);

fclose(fpr);

mlon = (float *)check_malloc(nx*local_ny*sizeof(float));
mlat = (float *)check_malloc(nx*local_ny*sizeof(float));

fpr = fopfile(cordfile,"r");

for(i=0;i<nx*iy0;i++)    /* scan thru file to line iy0*nx */
   fgets(str,512,fpr);

for(i=0;i<nx*local_ny;i++)   /* read-in (lat,lon) points for this processor */
   {
   fgets(str,512,fpr);
   sscanf(str,"%f %f",&mlon[i],&mlat[i]);
   }

fclose(fpr);

fprintf(stderr,"      Model slice=  %d (nx*nz) points\n",nn);
fprintf(stderr,"              iy0=  %d\n",iy0);
fprintf(stderr,"              iy1=  %d\n\n",iy1);
fprintf(stderr,"      Approximate memory (RAM)= %.2f Mb\n\n",(9*nn + 2*nx*local_ny + nz)*sizeof(float)*1.0e-06);
fflush(stderr);

fbuf = (float *) check_malloc(3*nn*sizeof(float));

sprintf(outfile,"%s/%s-%.4d.cvm",outdir,fileroot,myid);

fprintf(stderr,"      Output file= %s\n\n",outfile);
fflush(stderr);

fdw = croptrfile(outfile);

rite(fdw,&nx,sizeof(int));
rite(fdw,&ny,sizeof(int));
rite(fdw,&nz,sizeof(int));
rite(fdw,&iy0,sizeof(int));
rite(fdw,&local_ny,sizeof(int));

double rloc[3];
float* rtop = (float*)check_malloc(nx*local_ny*sizeof(float));
float surface_elev;
 float topo, mtop;

//epsilons taken from Patrick's code
float TOP_MODEL_EPSILON, ELEV_EPSILON = 0.01;
float CM_VOXEL_HEIGHT = 1000.0;
float LR_HR_VOXEL_HEIGHT = 100.0;

 //fprintf(stderr, "[%d] Sleeping for %d seconds\n", myid, myid*2);
 //fflush(stderr);
 //sleep(myid);

//run ucvm setup
 fprintf(stderr, "[%d] Setting up ucvm\n", myid);
 fflush(stderr);
 if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
   fprintf(stderr, "Failed to setup ucvm.\n");
   fflush(stderr);
   exit(-1);
 }

 //separate model names and add them
 if (strstr(models, ",")!=NULL) {
	 int ucvm_initialized = 0;
	 int ucvm_no_gtl_initialized = 0;
	 char* tok = strtok(models, ",");
	 while (tok!=NULL) {
		if (strcmp(tok, "cvmh")==0) {
			if (!ucvm_initialized) {
				if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
				   fprintf(stderr, "Failed to setup ucvm.\n");
				   fflush(stderr);
				   exit(-1);
				}
				ucvm_initialized = 1;
			}
			if (ucvm_no_gtl_initialized) {
				fprintf(stderr, "Can't use cvmh with a GTL after initializing UCVM without the GTL.");
				fflush(stderr);
				exit(-2);
			}
			if (ucvm_add_model(UCVM_MODEL_CVMH)!=UCVM_CODE_SUCCESS) {
			   fprintf(stderr, "Error retrieving CVM-H.\n");
			   fflush(stderr);
			   exit(-1);
			}
		} else if (strcmp(tok, "cvms")==0) {
			if (!ucvm_initialized) {
	                        if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
	                           fprintf(stderr, "Failed to setup ucvm.\n");
	                           fflush(stderr);
	                           exit(-1);
	                        }
				ucvm_initialized = 1;
			}
			if (ucvm_add_model(UCVM_MODEL_CVMS)!=UCVM_CODE_SUCCESS) {
			   fprintf(stderr, "Error retrieving CVM-S.\n");
	                   fflush(stderr);
			   exit(-2);
	                }
		} else if (strcmp(tok, "1d")==0) {
			if (!ucvm_initialized) {
	                        if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
	                           fprintf(stderr, "Failed to setup ucvm.\n");
	                           fflush(stderr);
	                           exit(-1);
	                        }
				ucvm_initialized = 1;
			}
			if (ucvm_add_model(UCVM_MODEL_1D) != UCVM_CODE_SUCCESS) {
			   fprintf(stderr, "Error retrieving 1D model.\n");
			   fflush(stderr);
			   exit(-3);
			}
		} else if (strcmp(tok, "cvmsi")==0) {
			if (!ucvm_initialized) {
	                        if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
	                           fprintf(stderr, "Failed to setup ucvm.\n");
	                           fflush(stderr);
	                           exit(-1);
				}
				ucvm_initialized = 1;
                        }
			if (ucvm_add_model(UCVM_MODEL_CVMSI) != UCVM_CODE_SUCCESS) {
	  		   fprintf(stderr, "Error retrieving CVM-SI.\n");
                           fflush(stderr);
                           exit(-4);
			}
		} else if (strcmp(tok, "hk")==0) {
			//Hadley-Kanamori
			if (!ucvm_initialized) {
				if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
	  			   fprintf(stderr, "Failed to setup ucvm.\n");
                                   fflush(stderr);
                                   exit(-1);
                                }
                                ucvm_initialized = 1;
			}
			if (ucvm_add_model(UCVM_MODEL_1D) != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving Hadley-Kanamori model.\n");
                           fflush(stderr);
                           exit(-3);
                        }
		} else if (strcmp(tok, "cvmh_nogtl")==0) {
			if (ucvm_initialized) {
				fprintf(stderr, "Error: need to put cvmh_nogtl first in list of models.");
				fflush(stderr);
				exit(-2);
			}
			ucvm_initialized = 1;
			ucvm_no_gtl_initialized = 1;	
			if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/cvmh_no_gtl.conf") != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Failed to setup ucvm.\n");
                           fflush(stderr);
                           exit(-1);
                        }
                        if (ucvm_add_model(UCVM_MODEL_CVMH) != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving CVM-H (no GTL).\n");
                           fflush(stderr);
                           exit(-4);
                        }
		} else {
			fprintf(stderr, "Model %s didn't match any known models, aborting.\n", tok);
			fflush(stderr);
			exit(-4);
		}
		tok = strtok(NULL, ",");
	}
 } else {
	if (strcmp(models, "cvmh")==0) {
            if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Failed to setup ucvm.\n");
                fflush(stderr);
                exit(-1);
             }
             if (ucvm_add_model(UCVM_MODEL_CVMH)!=UCVM_CODE_SUCCESS) {
             	fprintf(stderr, "Error retrieving CVM-H.\n");
             	fflush(stderr);
		exit(-1);
             }
        } else if (strcmp(models, "cvms")==0) {
            if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Failed to setup ucvm.\n");
                fflush(stderr);
                exit(-1);
             }
   	     if (ucvm_add_model(UCVM_MODEL_CVMS)!=UCVM_CODE_SUCCESS) {
         	    fprintf(stderr, "Error retrieving CVM-S.\n");
                    fflush(stderr);
		    exit(-2);
             }
        } else if (strcmp(models, "1d")==0) {
            if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Failed to setup ucvm.\n");
                fflush(stderr);
                exit(-1);
             }
             if (ucvm_add_model(UCVM_MODEL_1D)!=UCVM_CODE_SUCCESS) {
	             fprintf(stderr, "Error retrieving 1D model.\n");
                     fflush(stderr);
                     exit(-3);
             }
        } else if (strcmp(models, "cvmsi")==0) {
            if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Failed to setup ucvm.\n");
                fflush(stderr);
                exit(-1);
             }
             if (ucvm_add_model(UCVM_MODEL_CVMSI) != UCVM_CODE_SUCCESS) {
                     fprintf(stderr, "Error retrieving CVM-SI.\n");
                     fflush(stderr);
                     exit(-4);
             }
	} else if (strcmp(models, "hk")==0) {
		    if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
			fprintf(stderr, "Failed to setup ucvm.\n");
			fflush(stderr);
			exit(-1);
		     }
		     if (ucvm_add_model(UCVM_MODEL_1D) != UCVM_CODE_SUCCESS) {
			     fprintf(stderr, "Error retrieving 1D model.\n");
			     fflush(stderr);
			     exit(-4);
		     }
	} else if (strcmp(models, "cvmh_nogtl")==0) {
		     if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_13.9.0/conf/cvmh_no_gtl.conf") != UCVM_CODE_SUCCESS) {
                        fprintf(stderr, "Failed to setup ucvm.\n");
                        fflush(stderr);
                        exit(-1);
                     }
                     if (ucvm_add_model(UCVM_MODEL_CVMH) != UCVM_CODE_SUCCESS) {
                             fprintf(stderr, "Error retrieving 1D model.\n");
                             fflush(stderr);
                             exit(-4);
                     }
	} else {
 	       fprintf(stderr, "Model %s didn't match any known models, aborting.\n", models);
               fflush(stderr);
               exit(-3);
        }
 }

 // Query by depth
 if (ucvm_setparam(UCVM_PARAM_QUERY_MODE, UCVM_COORD_GEO_DEPTH)!=UCVM_CODE_SUCCESS) {
   fprintf(stderr, "Set query mode by depth failed.\n");
   fflush(stderr);
   exit(-2);
 }

 fprintf(stderr, "[%d] Generating CVM\n", myid);
 fflush(stderr);

 // output order: x,z,y
 int num_pts = nx;
 ucvm_point_t* pts = check_malloc(sizeof(ucvm_point_t)*num_pts);
 ucvm_data_t* props = check_malloc(sizeof(ucvm_data_t)*num_pts);

 for(j=0;j<local_ny;j++) {

   fprintf(stderr," %5d of %5d: generating CVM ... ",j+1,local_ny);
   fflush(stderr);

   for(iz=0;iz<nz;iz++) {

     for(ix=0; ix<nx; ix++) {
       i = ix + j*nx; // index into elevation map, latlon
       //k = ix + iz*nx; // index into output array

       pts[ix].coord[0] = mlon[i];
       pts[ix].coord[1] = mlat[i];
       pts[ix].coord[2] = mdep[iz];
     }       

     if (ucvm_query(num_pts, pts, props)!=UCVM_CODE_SUCCESS) {
	 fprintf(stderr, "Query UCVM failed.\n");
	 exit(-2);
     }

       /* perform sanity checks on the material properties */
     for (ix=0; ix<nx; ix++) { 
       k = ix + iz*nx;
       if (isnan(props[ix].cmb.vp) || isnan(props[ix].cmb.vs) || isnan(props[ix].cmb.rho) ||
	   isinf(props[ix].cmb.vp) || isinf(props[ix].cmb.vs) || isinf(props[ix].cmb.rho)) {
	 fprintf(stderr, "NaN/Inf detected at (%f, %f, %f)\n", 
		 pts[ix].coord[0], pts[ix].coord[1], pts[ix].coord[2]);
	 fflush(stderr);
	 exit(-1);
       }
       if ((props[ix].cmb.vp < 0.0) || (props[ix].cmb.vs < 0.0) || (props[ix].cmb.rho < 0.0)) {
	 fprintf(stderr, "Negative vals detected at (%f, %f, %f)\n", 
		 pts[ix].coord[0], pts[ix].coord[1], pts[ix].coord[2]);
	 fprintf(stderr, "vp=%f, vs=%f, rho=%lf\n", 
		 props[ix].cmb.vp, props[ix].cmb.vs, props[ix].cmb.rho);
	 fflush(stderr);
	 exit(-1);
       }

       if ((props[ix].cmb.vp / props[ix].cmb.vs) < MIN_V_RATIO) {
	 fprintf(stderr, 
		 "Vp/Vs < %f at (%f  %f  %f) (vp=%f,vs=%f)\n", 
		 MIN_V_RATIO,
		 pts[ix].coord[0], pts[ix].coord[1], pts[ix].coord[2], 
		 props[ix].cmb.vp, props[ix].cmb.vs);
	 fflush(stderr);
	 fprintf(ppr, 
		 "%f  %f  %f\n", pts[ix].coord[0], pts[ix].coord[1], pts[ix].coord[2]);
	 fflush(ppr);
	 //exit(-1);
       }

       fbuf[k] = meter2km*props[ix].cmb.vp;
       fbuf[k+nn] = meter2km*props[ix].cmb.vs;
       fbuf[k+2*nn] = meter2km*props[ix].cmb.rho;
     }
   }

   fprintf(stderr,"writing output ... ");
   fflush(stderr);

   rite(fdw,fbuf,3*nn*sizeof(float));

   fprintf(stderr,"DONE\n");
   fflush(stderr);
 }

free(pts);
free(props);

fprintf(stderr, "[%d] Closing output file\n", myid);
fflush(stderr);

close(fdw);
fclose(ppr);

fprintf(stderr, "[%d] Successful exit\n", myid);
fflush(stderr);

mpi_exit(0);
}
