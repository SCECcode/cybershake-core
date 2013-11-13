#include "mpi.h"
#include "include.h"
#include "function.h"
#include "func_mpi.h"
#include "ucvm.h"


/* Vp/Vs ratio */
#define MIN_V_RATIO 1.0
#define VX_NO_DATA -99998.0
#define START 0
#define END 1
#define RWG 0
#define AWP 1

int main(int ac,char **av)
{
/* MPI-RWG: MPI stuff and distributed computation variables */
int my_id, nproc, pnlen;
char procname[128];

 FILE *fpr, *ppr;
float *mlon, *mlat, *mdep;
float *vp_buf, *vs_buf, *rho_buf, *buf;
vp_buf = vs_buf = rho_buf = buf = NULL;

int fdw, i, j, k, nn, newn, s;
int nx, ny, nz, ix, iz, icnt;
int format;
//My processor coordinates
int my_nx, my_ny, my_nz;
//My endpoint bounds
int starting_stripe, ending_stripe;
int local_np;

char infile[512], cordfile[512], depfile[512], outfile[512], modeldir[512], models[512];
char logdir[512], logname[512], logfile[512], pointfile[512];
char str[512], format_name[16];

float meter2feet = 3.2808399;
float meter2km = 0.001;
int tomobg = 1;
int geot_layer = 1;
int max_iter = 20;

//In m/s (RWG expects km/s)
float min_vp, min_vs, min_rho;

mpi_init(&ac,&av,&nproc,&my_id,procname,&pnlen);

sprintf(logdir,".");
sprintf(logname,"v_mpi");

setpar(ac,av);
mstpar("nx","d",&nx);
mstpar("ny","d",&ny);
mstpar("nz","d",&nz);

mstpar("cordfile","s",cordfile);
mstpar("depfile","s",depfile);
mstpar("modeldir","s",modeldir);

mstpar("outfile","s",outfile);

//These should be in m/s
mstpar("min_vp","f",&min_vp);
mstpar("min_vs","f",&min_vs);
mstpar("min_rho","f",&min_rho);

fprintf(stderr, "min_vp=%f, min_vs=%f, min_rho=%f\n", min_vp, min_vs, min_rho);

getpar("logdir","s",logdir);
getpar("logname","s",logname);

getpar("tomobg","d",&tomobg);
getpar("geot_layer","d",&geot_layer);
getpar("models","s",models);

strcpy(format_name, "rwg");
getpar("format","s",format_name);

endpar();

//For debugging
sprintf(logfile,"%s/%s-%.4d.log", logdir, logname, my_id);
freopen(logfile,"w",stderr);

if (strcmp(format_name, "rwg")==0) {
	format = RWG;
} else if (strcmp(format_name, "awp")==0) {
	format = AWP;
} else {
	fprintf(stderr, "Format %s is not recognized.", format_name);
	exit(2);
}

makedir(logdir);

fprintf(stderr,"***** This is UCVM-mpi\n");
fprintf(stderr,"      Running on node= %s rank= %d\n\n",procname,my_id);
fflush(stderr);

//Define my coordinates, depending on format
if (format==RWG) {
	//fast x, z, y	
	//Assume we have few enough processors that everyone gets an X-stripe
	int num_x_stripes = ny*nz;
	int x_stripes_per_proc = (ny*nz)/nproc;
	if (x_stripes_per_proc==0) {
		fprintf(stderr, "%d) Trying to use too many processors.  Use no more than %d processors.\n", my_id, num_x_stripes);
		mpi_exit(2);
	} else if (nproc*x_stripes_per_proc!=ny*nz) {
		fprintf(stderr, "%d) Error - number of processors must be a factor of ny*nz.  %d must be a factor of %d.\n", my_id, nproc, num_x_stripes);
		mpi_exit(3);
	}
	
	//Determine my starting, ending y and z values
	starting_stripe = my_id*x_stripes_per_proc;
	ending_stripe = (my_id+1)*x_stripes_per_proc;
	if (ending_stripe>num_x_stripes) {
		ending_stripe = num_x_stripes;
	}
	local_np = nx * (ending_stripe - starting_stripe);

} else {
	//fast y, x, z
	//Everyone gets a y-stripe here
	int num_y_stripes = nx*nz;
	int y_stripes_per_proc = (nx*nz)/nproc;
        if (y_stripes_per_proc==0) {
                fprintf(stderr, "%d) Trying to use too many processors.  Use no more than %d processors.\n", my_id, num_y_stripes);
                mpi_exit(2);
        } else if (nproc*y_stripes_per_proc!=nx*nz) {
                fprintf(stderr, "%d) Error - number of processors must be a factor of nx*nz.  %d must be a factor of %d.\n", my_id, nproc, num_y_stripes);
		mpi_exit(3);
        }

	starting_stripe = my_id*y_stripes_per_proc;
        ending_stripe = (my_id+1)*y_stripes_per_proc;
        if (ending_stripe>num_y_stripes) {
                ending_stripe = num_y_stripes;
        }

	local_np = ny * (ending_stripe - starting_stripe);
}

mlon = check_malloc(nx*ny*sizeof(float));
mlat = check_malloc(nx*ny*sizeof(float));
mdep = check_malloc(nz*sizeof(float));

if (my_id==0) {
	//Read in the coords and the depths and broadcast
	fpr = fopfile(depfile,"r");
	for(iz=0;iz<nz;iz++) {
	   fscanf(fpr,"%f",&mdep[iz]);
	}
	fclose(fpr);

	fpr = fopfile(cordfile,"r");
	for(i=0;i<nx*ny;i++)   /* read-in (lat,lon) points */
	   {
	   fgets(str,512,fpr);
	   sscanf(str,"%f %f",&mlon[i],&mlat[i]);
	   }
	fclose(fpr);
}

//Broadcast depths, lats, lons
int error = MPI_Bcast(mdep, nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
if (error!=MPI_SUCCESS) {
	char string[256];
        int err_len;
        MPI_Error_string(error, string, &err_len);
        fprintf(stderr, "Process %d: error in mdep bcast\n", my_id);
        fprintf(stderr, "Error message: %s\n", string);
        MPI_Finalize();
        exit(1);
}
error = MPI_Bcast(mlat, nx*ny, MPI_FLOAT, 0, MPI_COMM_WORLD);
if (error!=MPI_SUCCESS) {
        char string[256];
        int err_len;
        MPI_Error_string(error, string, &err_len);
        fprintf(stderr, "Process %d: error in mlat bcast\n", my_id);
        fprintf(stderr, "Error message: %s\n", string);
        MPI_Finalize();
        exit(1);
}
error = MPI_Bcast(mlon, nx*ny, MPI_FLOAT, 0, MPI_COMM_WORLD);
if (error!=MPI_SUCCESS) {
        char string[256];
        int err_len;
        MPI_Error_string(error, string, &err_len);
        fprintf(stderr, "Process %d: error in mlon bcast\n", my_id);
        fprintf(stderr, "Error message: %s\n", string);
        MPI_Finalize();
        exit(1);
}

fprintf(stderr,"%d)      Model slice=  %d points\n", my_id, local_np);
fprintf(stderr,"%d)        Stripes: %d -> %d\n", my_id, starting_stripe, ending_stripe);
//input arrays, plus UCVM arrays and output array
fprintf(stderr,"%d)      Approximate memory (RAM)= %.2f Mb\n\n",my_id,((2*nx*ny+nz)*sizeof(float) + (sizeof(ucvm_point_t)+sizeof(ucvm_data_t)+3*sizeof(float))*local_np)*1.0e-06);
fflush(stderr);

if (format==RWG) {
	vp_buf = check_malloc(local_np*sizeof(float));
	vs_buf = check_malloc(local_np*sizeof(float));
	rho_buf = check_malloc(local_np*sizeof(float));
} else if (format==AWP) {
	buf = check_malloc(3*local_np*sizeof(float));
}

//epsilons taken from Patrick's code
float TOP_MODEL_EPSILON, ELEV_EPSILON = 0.01;
float CM_VOXEL_HEIGHT = 1000.0;
float LR_HR_VOXEL_HEIGHT = 100.0;

//run ucvm setup
 fprintf(stderr, "[%d] Setting up ucvm\n", my_id);
 fflush(stderr);
 if (ucvm_init("/work/00940/tera3d/CyberShake/software/UCVM/ucvm_12.2.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
   fprintf(stderr, "Failed to setup ucvm.\n");
   fflush(stderr);
   exit(-1);
 }

// Query by depth
 if (ucvm_setparam(UCVM_PARAM_QUERY_MODE, UCVM_COORD_GEO_DEPTH)!=UCVM_CODE_SUCCESS) {
   fprintf(stderr, "Set query mode by depth failed.\n");
   fflush(stderr);
   exit(-2);
 }
 //Add CVM-H models
 //separate model names and add them
 if (strstr(models, ",")!=NULL) {
	 char* tok = strtok(models, ",");
	 while (tok!=NULL) {
		if (strcmp(tok, "cvmh")==0) {
			if (ucvm_add_model(UCVM_MODEL_CVMH)!=UCVM_CODE_SUCCESS) {
			   fprintf(stderr, "Error retrieving CVM-H.\n");
			   fflush(stderr);
			   exit(-1);
			}
		} else if (strcmp(tok, "cvms")==0) {
			if (ucvm_add_model(UCVM_MODEL_CVMS)!=UCVM_CODE_SUCCESS) {
			   fprintf(stderr, "Error retrieving CVM-S.\n");
	                   fflush(stderr);
			   exit(-2);
	                }
		} else if (strcmp(tok, "1d")==0) {
			if (ucvm_add_model(UCVM_MODEL_1D) != UCVM_CODE_SUCCESS) {
			   fprintf(stderr, "Error retrieving 1D model.\n");
			   fflush(stderr);
			   exit(-3);
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
             if (ucvm_add_model(UCVM_MODEL_CVMH)!=UCVM_CODE_SUCCESS) {
             	fprintf(stderr, "Error retrieving CVM-H.\n");
             	fflush(stderr);
		exit(-1);
             }
        } else if (strcmp(models, "cvms")==0) {
   	     if (ucvm_add_model(UCVM_MODEL_CVMS)!=UCVM_CODE_SUCCESS) {
         	    fprintf(stderr, "Error retrieving CVM-S.\n");
                    fflush(stderr);
		    exit(-2);
             }
        } else if (strcmp(models, "1d")==0) {
             if (ucvm_add_model(UCVM_MODEL_1D)!=UCVM_CODE_SUCCESS) {
	             fprintf(stderr, "Error retrieving 1D model.\n");
                     fflush(stderr);
                     exit(-3);
             }
	} else {
 	       fprintf(stderr, "Model %s didn't match any known models, aborting.\n", models);
               fflush(stderr);
               exit(-3);
        }
 }


 fprintf(stderr, "[%d] Generating CVM\n", my_id);
 fflush(stderr);

 // output order: x,z,y
 int num_pts = local_np;
 ucvm_point_t* pts = check_malloc(sizeof(ucvm_point_t)*num_pts);
 ucvm_data_t* props = check_malloc(sizeof(ucvm_data_t)*num_pts);

 //put the if statement out here, for easier optimization
 if (format==RWG) {
	//fast x, z, y
	for(s=starting_stripe; s<ending_stripe; s++) {
		//Get z and y values from stripe
		int z_ind = s % nz;
		int y_ind = s / nz;
		for(i=0; i<nx; i++) {
			//offset in coordinate list (fast x, y)
			int input_offset = y_ind*nx + i;
			int output_offset = (s-starting_stripe)*nx+i;
			pts[output_offset].coord[0] = mlon[input_offset];
                        pts[output_offset].coord[1] = mlat[input_offset];
                        pts[output_offset].coord[2] = mdep[z_ind];
		}
	}

 } else {
	//AWP format
	//output is fast y, x, z
	for(s=starting_stripe; s<ending_stripe; s++) {
		int x_ind = s % nx;
		int z_ind = s / nx;
		for(j=0; j<ny; j++) {
			int input_offset = j*nx + x_ind;
			int output_offset = (s-starting_stripe)*ny + j;
			pts[output_offset].coord[0] = mlon[input_offset];
                        pts[output_offset].coord[1] = mlat[input_offset];
                        pts[output_offset].coord[2] = mdep[z_ind];
		}
	}
 }

 if (ucvm_query(local_np, pts, props)!=UCVM_CODE_SUCCESS) {
 	fprintf(stderr, "Query UCVM failed.\n");
        exit(-2);
 }

 for (i=0; i<nx; i++) {
	fprintf(stderr, "(%f, %f, %f)->surf=%f, vs30=%f, depth=%f, cmb.vp=%f, cmb.vs=%f, cmb.rho=%f\n", pts[i].coord[0], pts[i].coord[1], pts[i].coord[2], props[i].surf, props[i].vs30, props[i].depth, props[i].cmb.vp, props[i].cmb.vs, props[i].cmb.rho);
 }

 /* perform sanity checks on the material properties */     

 fprintf(stderr, "Enforcing min_vp=%f, min_vs=%f, min_rho=%f.\n", min_vp, min_vs, min_rho);

 for (i=0; i<local_np; i++) { 
	//Check for nan, inf
       if (isnan(props[i].cmb.vp) || isnan(props[i].cmb.vs) || isnan(props[i].cmb.rho) ||
	   isinf(props[i].cmb.vp) || isinf(props[i].cmb.vs) || isinf(props[i].cmb.rho)) {
	 fprintf(stderr, "NaN/Inf detected at (%f, %f, %f)\n", 
		 pts[i].coord[0], pts[i].coord[1], pts[i].coord[2]);
	 fflush(stderr);
	 exit(-1);
       }
	//Check for negative values
       if ((props[i].cmb.vp < 0.0) || (props[i].cmb.vs < 0.0) || (props[i].cmb.rho < 0.0)) {
	 fprintf(stderr, "Negative vals detected at (%f, %f, %f)\n", 
		 pts[i].coord[0], pts[i].coord[1], pts[i].coord[2]);
	 fprintf(stderr, "vp=%f, vs=%f, rho=%lf\n", 
		 props[i].cmb.vp, props[i].cmb.vs, props[i].cmb.rho);
	 fflush(stderr);
	 exit(-1);
       }

	//Check for min Vp, Vs, Rho
	if (props[i].cmb.vp<min_vp) {
		props[i].cmb.vp = min_vp;
	}
	if (props[i].cmb.vs<min_vs) {
		props[i].cmb.vs = min_vs;
	}
	if (props[i].cmb.rho<min_rho) {
		props[i].cmb.rho = min_rho;
	}
	
	//Check poisson ratio
	float fac = 1.45;
	if (props[i].cmb.vp/props[i].cmb.vs < fac) {
		fprintf(stderr, "Adjusting index %d - changing vs from %f to %f.\n", i, props[i].cmb.vs, props[i].cmb.vp/fac);
		props[i].cmb.vs = props[i].cmb.vp/fac;
	}

	if (format==RWG) {
	        vp_buf[i] = meter2km*props[i].cmb.vp;
                vs_buf[i] = meter2km*props[i].cmb.vs;
        	rho_buf[i] = meter2km*props[i].cmb.rho;
        } else if (format==AWP) {
		buf[3*i] = props[i].cmb.vp;
                buf[3*i+1] = props[i].cmb.vs;
                buf[3*i+2] = props[i].cmb.rho;
	}

   }

  //Now, everyone opens and writes to the files
 if (format==RWG) {
	MPI_Info info;
	MPI_File fp_out;
	MPI_Offset offset;
	MPI_Status status;
	char filename[512];
	char suffix[] = {'p', 's', 'd'};
	float* bufs[] = {vp_buf, vs_buf, rho_buf};
	for (i=0; i<3; i++) {
		sprintf(filename, "%s.%c", outfile, suffix[i]);
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_out);
		//output is fast x, z, y, but contiguous
		offset = (long)starting_stripe * nx * sizeof(float);
		fprintf(stderr, "%d) writing at offset %ld\n", my_id, offset);
		int rc = MPI_File_write_at_all(fp_out, offset, bufs[i], local_np, MPI_FLOAT, &status);
		if (rc!=MPI_SUCCESS) {
			char error_string[256];
			int len_err_string;
			MPI_Error_string(rc, error_string, &len_err_string);
			fprintf(stderr, "%d) Error writing to file %s: %s\n", my_id, filename, error_string);
			MPI_Finalize();
			exit(3);
		}
		//Close
		MPI_File_close(&fp_out);
	}
 } else if (format==AWP) {
	MPI_Info info;
        MPI_File fp_out;
        MPI_Offset offset;
        MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_out);
	//output is fast y, x, z, all 3 values
	offset = (long)starting_stripe * ny * 3 * sizeof(float);
	int rc = MPI_File_write_at_all(fp_out, offset, buf, 3*local_np, MPI_FLOAT, &status);
        if (rc!=MPI_SUCCESS) {
        	char error_string[256];
                int len_err_string;
                MPI_Error_string(rc, error_string, &len_err_string);
                fprintf(stderr, "%d) Error writing to file %s: %s\n", my_id, outfile, error_string);
                mpi_exit(3);
        }
	MPI_File_close(&fp_out);
 }

free(mlat);
free(mlon);
free(mdep);

free(pts);
free(props);

free(vs_buf);
free(vp_buf);
free(rho_buf);
free(buf);

fflush(stderr);

mpi_exit(0);
}