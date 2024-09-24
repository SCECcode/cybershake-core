#define _GNU_SOURCE
#include "mpi.h"
#include "include.h"
#include "function.h"
#include "func_mpi.h"
#include "ucvm.h"
#include "ucvm_interp.h"
#include "ucvm_model_elygtl.h"

/* Vp/Vs ratio */
#define MIN_V_RATIO 1.0
#define VX_NO_DATA -99998.0
#define START 0
#define END 1
#define RWG 0
#define AWP 1
#define AWP_Z 2

int main(int ac,char **av)
{
/* MPI-RWG: MPI stuff and distributed computation variables */
int my_id, nproc, pnlen;
char procname[128];

 FILE *fpr, *ppr;
float *mlon, *mlat, *mdep;
float *vp_buf, *vs_buf, *rho_buf, *buf;
vp_buf = vs_buf = rho_buf = buf = NULL;

char* config_filename = "/lustre/orion/proj-shared/geo156/CyberShake/software/UCVM/ucvm_22.7.0_withSFCVM/conf/ucvm.conf";

int fdw, i, j, k, nn, newn, s;
int nx, ny, nz, ix, iz, icnt;
int format;
//My processor coordinates
int my_nx, my_ny, my_nz;
//My endpoint bounds
int starting_stripe, ending_stripe;
int pts_per_stripe;
int local_np;

char infile[512], cordfile[512], depfile[512], outfile[512], modeldir[512], models[512];
char logdir[512], logname[512], logfile[512], pointfile[512];
char str[512], format_name[16];

float meter2feet = 3.2808399;
float meter2km = 0.001;
int tomobg = 1;
int geot_layer = 1;
int max_iter = 20;

//Apply a constant mantle below some depth
int const_mantle = 0;
float const_mantle_depth = 45000.0;

//Depth, in meters, to use when querying UCVM to populate the surface point. So 25 means use a depth of 25m for the surface point, instead of 0.
float surface_cvm_depth = 0.0;

//In m/s (RWG expects km/s)
float min_vp, min_vs, min_rho;

//When using AWP-Z mode, apply load-balancing
int load_balancing = 0;

int z_query_mode = UCVM_COORD_GEO_DEPTH;

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

getpar("const_mantle", "d", &const_mantle);
getpar("const_mantle_depth", "f", &const_mantle_depth);

//When populating the surface points, use this depth in meters instead
getpar("surface_cvm_depth", "f", &surface_cvm_depth);

//Option for using the Ely taper:
//'none': don't use (default)
//'all': always use
//'ifless': use which ever is smaller, down to transition depth
char ely_taper[64];
sprintf(ely_taper, "none");
getpar("ely_taper", "s", ely_taper);
//Transition depth for taper
float ely_transition_value = 700.0;
getpar("ely_transition_depth", "f", &ely_transition_value);
//List of comma-separated models to apply taper to, default is none
char ely_taper_models[256];
ely_taper_models[0]='\0';
getpar("ely_taper_models", "s", ely_taper_models);

char z_query_mode_string[64];
getpar("z_query_mode", "s", z_query_mode_string);
if (strcmp(z_query_mode_string, "depth")==0) {
	z_query_mode = UCVM_COORD_GEO_DEPTH;
} else if (strcmp(z_query_mode_string, "elevation")==0) {
	z_query_mode = UCVM_COORD_GEO_ELEV;
}

//If using query by elevation, invert transition value
if (z_query_mode==UCVM_COORD_GEO_ELEV) {
	ely_transition_value *= -1.0;
}

getpar("load_balancing", "d", &load_balancing);

endpar();

//For debugging
sprintf(logfile,"%s/%s-%.4d.log", logdir, logname, my_id);
freopen(logfile,"w",stderr);

if (strcmp(format_name, "rwg")==0) {
	format = RWG;
} else if (strcmp(format_name, "awp")==0) {
	format = AWP;
} else if (strcmp(format_name, "awpz")==0) {
    format = AWP_Z;
} else {
	fprintf(stderr, "Format %s is not recognized.", format_name);
	exit(2);
}

makedir(logdir);

fprintf(stderr,"***** This is UCVM-mpi\n");
fprintf(stderr,"      Running on node= %s core=%d rank= %d\n\n",procname,sched_getcpu(), my_id);
fflush(stderr);

//Define my coordinates, depending on format
if (format==RWG) {
	//fast x, z, y	
	//Assume we have few enough processors that everyone gets an X-stripe
	int num_x_stripes = ny*nz;
	float x_stripes_per_proc = ((float)(ny*nz))/((float)nproc);
	//int x_stripes_per_proc = (ny*nz)/nproc;
	if (x_stripes_per_proc==0) {
		fprintf(stderr, "%d) Trying to use too many processors.  Use no more than %d processors.\n", my_id, num_x_stripes);
		mpi_exit(2);
	}
	//Modifying to support non-even divisions
	// else if (nproc*x_stripes_per_proc!=ny*nz) {
	//	fprintf(stderr, "%d) Error - number of processors must be a factor of ny*nz.  %d must be a factor of %d.\n", my_id, nproc, num_x_stripes);
	//	mpi_exit(3);
	//}
	
	//Determine my starting, ending y and z values
	starting_stripe = (int)(my_id*x_stripes_per_proc);
	ending_stripe = (int)((my_id+1)*x_stripes_per_proc);
	if (ending_stripe>num_x_stripes) {
		ending_stripe = num_x_stripes;
	}
	if (my_id==nproc-1) {
		ending_stripe = num_x_stripes;
	}
	pts_per_stripe = nx;
	local_np = pts_per_stripe * (ending_stripe - starting_stripe);

} else if (format==AWP_Z) {
	//output is fast y, x, z
	//But we let Z vary fastest in the query, then do some memory shuffling
	int num_z_stripes = nx*ny;
	float z_stripes_per_proc = ((float)(nx*ny))/((float)nproc);
	if (z_stripes_per_proc==0) {
		fprintf(stderr, "%d) Trying to use too many processors.  Use no more than %d processors.\n", my_id, num_z_stripes);
	    mpi_exit(2); 
	}
	starting_stripe = (int)(my_id*z_stripes_per_proc);
	ending_stripe = (int)((my_id+1)*z_stripes_per_proc);
	if (ending_stripe>num_z_stripes) {
		ending_stripe = num_z_stripes;
	}
	if (my_id==nproc-1) {
		ending_stripe = num_z_stripes;
	}
	pts_per_stripe = nz;
	local_np = pts_per_stripe * (ending_stripe - starting_stripe);
} else {
	//fast y, x, z
	//Everyone gets a y-stripe here
	int num_y_stripes = nx*nz;
	float y_stripes_per_proc = ((float)(nx*nz))/((float)nproc);
	//int y_stripes_per_proc = (nx*nz)/nproc;
        if (y_stripes_per_proc==0) {
                fprintf(stderr, "%d) Trying to use too many processors.  Use no more than %d processors.\n", my_id, num_y_stripes);
                mpi_exit(2);
        }
	//Modifying to support non-even divisions
	//else if (nproc*y_stripes_per_proc!=nx*nz) {
        //        fprintf(stderr, "%d) Error - number of processors must be a factor of nx*nz.  %d must be a factor of %d.\n", my_id, nproc, num_y_stripes);
	//	mpi_exit(3);
        //}

	starting_stripe = (int)(my_id*y_stripes_per_proc);
        ending_stripe = (int)((my_id+1)*y_stripes_per_proc);
        if (ending_stripe>num_y_stripes) {
                ending_stripe = num_y_stripes;
        }
	if (my_id==nproc-1) {
		ending_stripe = num_y_stripes;
	}
	pts_per_stripe = ny;
	local_np = pts_per_stripe * (ending_stripe - starting_stripe);
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

if (format==AWP_Z && load_balancing==1) {
	if (my_id==0) {
		//Reassign stripes based on cost
		//Stripes are Z-columns, with fast y, x 
		float* col_cost = check_malloc(sizeof(float)*nx*ny);
		float lat, lon;
		//Col index is fast y, x
		int col_index;
		//Coordinate list is fast x, y
		int coord_index;
		double total_cost = 0.0;
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) {
				col_index = i*ny + j;
				coord_index = j*nx + i;
				lat = mlat[coord_index];
				lon = mlon[coord_index];
				//Check if in ocean
				if (lat < -1.4922 * lon - 145.3135) {
					 col_cost[col_index] = 10.0;
				//Check if outside of model to the east (1D)
				} else if (lat > -1.1052 * lon - 94.7573) {
					col_cost[col_index] = 0.1;
				//Check if outside of model to the south
				} else if (lat < -0.56894 * lon + 104.3751) {
					col_cost[col_index] = 0.1;
				} else {
					col_cost[col_index] = 1.0;
				}
				total_cost += col_cost[col_index];
			}
		}
		printf("Total cost = %f, per_proc_cost = %f\n", total_cost, total_cost/nproc);
		//Now, divide them up
		//Array has list of stripe cutoffs
		//[0, first one of P1, first one of P2, ...., N]
		int* stripe_array = check_malloc(sizeof(int)*(nproc+1));
		stripe_array[0] = 0;
		float per_proc_cost = total_cost/nproc;
		col_index = 0;
		int stripe_index = 1;
		double cur_cost = 0.0;
		while (col_index<nx*ny) {
			cur_cost += col_cost[col_index];
			if (cur_cost>=stripe_index*per_proc_cost) {
				stripe_array[stripe_index] = col_index;
				stripe_index++;
			}
			col_index++;
		}
		stripe_array[nproc] = nx*ny;
		for (i=0; i<nproc; i++) {
			printf("%d) proc %d: [%d-%d)\n", my_id, i, stripe_array[i], stripe_array[i+1]);
		}
		fflush(stdout);
		stripe_array[nproc] = nx*ny;
		int error = MPI_Bcast(stripe_array, nproc+1, MPI_INT, 0, MPI_COMM_WORLD);
		if (error!=MPI_SUCCESS) {
			char string[256];
			int err_len;
			MPI_Error_string(error, string, &err_len);
        	fprintf(stderr, "Process %d: error in stripe_array bcast\n", my_id);
        	fprintf(stderr, "Error message: %s\n", string);
        	MPI_Finalize();
        	exit(1);
		}
		starting_stripe = stripe_array[my_id];
		ending_stripe = stripe_array[my_id+1];
		local_np = pts_per_stripe * (ending_stripe - starting_stripe);
		free(col_cost);
		free(stripe_array);
	} else {
		int* stripe_array = check_malloc(sizeof(int)*(nproc+1));
		int error = MPI_Bcast(stripe_array, nproc+1, MPI_INT, 0, MPI_COMM_WORLD);
        if (error!=MPI_SUCCESS) {
            char string[256];
            int err_len;
            MPI_Error_string(error, string, &err_len);
            fprintf(stderr, "Process %d: error in stripe_array bcast\n", my_id);
            fprintf(stderr, "Error message: %s\n", string);
            MPI_Finalize();
            exit(1);
        }
        starting_stripe = stripe_array[my_id];
        ending_stripe = stripe_array[my_id+1];
        local_np = pts_per_stripe * (ending_stripe - starting_stripe);
        free(stripe_array);
	}
	printf("%d) pts_per_stripe = %d, starting_stripe = %d, ending_stripe = %d\n", my_id, pts_per_stripe, starting_stripe, ending_stripe);
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
fprintf(stderr,"%d)        Stripes: %d -> %d, %d points per stripe\n", my_id, starting_stripe, ending_stripe, pts_per_stripe);
//input arrays, plus UCVM arrays and output array
fprintf(stderr,"%d)      Approximate memory (RAM)= %.2f Mb\n\n",my_id,((2*nx*ny+nz)*sizeof(float) + (sizeof(ucvm_point_t)+3*sizeof(float))*local_np+pts_per_stripe*sizeof(ucvm_data_t))*1.0e-06);
fflush(stderr);

if (format==RWG) {
	vp_buf = check_malloc(local_np*sizeof(float));
	vs_buf = check_malloc(local_np*sizeof(float));
	rho_buf = check_malloc(local_np*sizeof(float));
} else if (format==AWP) {
	buf = check_malloc(3*local_np*sizeof(float));
} else if (format==AWP_Z) {
	buf = check_malloc(3*local_np*sizeof(float));
}

//epsilons taken from Patrick's code
float TOP_MODEL_EPSILON, ELEV_EPSILON = 0.01;
float CM_VOXEL_HEIGHT = 1000.0;
float LR_HR_VOXEL_HEIGHT = 100.0;

//run ucvm setup
 fprintf(stderr, "[%d] Setting up ucvm\n", my_id);
 fflush(stderr);

 //Add CVM-H models
 //separate model names and add them
 if (strstr(models, ",")!=NULL) {
	 int ucvm_initialized = 0;
	 int ucvm_no_gtl_initialized = 0;
	 char* save;
	 char* tok = strtok_r(models, ",", &save);
	 while (tok!=NULL) {
		if (strcmp(tok, "cvmh")==0) {
			if (ucvm_no_gtl_initialized) {
				fprintf(stderr, "Error:  Only cvmh or cvmh_nogtl may be used.\n");
				exit(-2);
			}
			if (!ucvm_initialized) {
		                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
		                   fprintf(stderr, "Failed to setup ucvm.\n");
		                   fflush(stderr);
		                   exit(-1);
	                        }
				ucvm_initialized = 1;
			}
			if (ucvm_add_model(UCVM_MODEL_CVMH)!=UCVM_CODE_SUCCESS) {
			   fprintf(stderr, "Error retrieving CVM-H.\n");
			   fflush(stderr);
			   exit(-1);
			}
		} else if (strcmp(tok, "cvms")==0) {
                        if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
                                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
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
                        if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
                                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
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
                        if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
                                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
                                   fprintf(stderr, "Failed to setup ucvm.\n");
                                   fflush(stderr);
                                   exit(-1);
                                }
                                ucvm_initialized = 1;
                        }
                        if (ucvm_add_model(UCVM_MODEL_CVMSI) != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving CVM-S4.26 model.\n");
                           fflush(stderr);
                           exit(-3);
                        }
                } else if (strcmp(tok, "scec1d")==0) {
                        if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
                                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
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
		} else if (strcmp(tok, "usgs")==0) {
			if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
                                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
                                   fprintf(stderr, "Failed to setup ucvm.\n");
                                   fflush(stderr);
                                   exit(-1);
                                }
                                ucvm_initialized = 1;
                        }
			if (ucvm_add_model(UCVM_MODEL_CENCAL) != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving USGS Bay Area model.\n");
                           fflush(stderr);
                           exit(-3);
                        } 
		} else if (strcmp(tok, "cvms5")==0) {
                        if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
                                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
                                   fprintf(stderr, "Failed to setup ucvm.\n");
                                   fflush(stderr);
                                   exit(-1);
                                }
                                ucvm_initialized = 1;
                        }
                        if (ucvm_add_model("cvms5") != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving CVM-S5 model.\n");
                           fflush(stderr);
                           exit(-3);
                        } else {
			   printf("Added CVM-S5.\n");
			}
                } else if (strcmp(tok, "cvmh_nogtl")==0) {
			if (ucvm_initialized) {
				fprintf(stderr, "Error: cvmh_nogtl must be first in list of models.\n");
				exit(-2);
			}
                        if (!ucvm_no_gtl_initialized) {
                                if (ucvm_init("/lustre/atlas/proj-shared/geo112/CyberShake/software/UCVM/ucvm-15.10.0/conf/cvmh_no_gtl.conf") != UCVM_CODE_SUCCESS) {
                                   fprintf(stderr, "Failed to setup ucvm.\n");
                                   fflush(stderr);
                                   exit(-1);
                                }
                                ucvm_no_gtl_initialized = 1;
                        }
                        if (ucvm_add_model(UCVM_MODEL_CVMH) != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving CVM-H (no GTL) model.\n"); 
                           fflush(stderr);
                           exit(-3);
                        }
		} else if (strcmp(tok, "bbp1d")==0) {
			if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
                                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
                                   fprintf(stderr, "Failed to setup ucvm.\n");
                                   fflush(stderr);
                                   exit(-1);
                                }
                                ucvm_initialized = 1;
                        }
                        if (ucvm_add_model(UCVM_MODEL_BBP1D) != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving BBP 1D model.\n");
                           fflush(stderr);
                         }
		} else if (strcmp(tok, "cca1d")==0) {
                        if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
                                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
                                   fprintf(stderr, "Failed to setup ucvm.\n");
                                   fflush(stderr);
                                   exit(-1);
                                }
                                ucvm_initialized = 1;
                        }
                        if (ucvm_add_model(UCVM_MODEL_BBP1D) != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving CCA 1D model.\n");
                           fflush(stderr);
                         }
                } else if (strcmp(tok, "cca")==0) {
                        if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
                                if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
                                   fprintf(stderr, "Failed to setup ucvm.\n");
                                   fflush(stderr);
                                   exit(-1);
                                }
                                ucvm_initialized = 1;
                        }
                        if (ucvm_add_model("cca") != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving CCA model.\n");
                           fflush(stderr);
                         }
		} else if (strcmp(tok, "sfcvm")==0) {
			if (!ucvm_initialized && !ucvm_no_gtl_initialized) {
				if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
                	fprintf(stderr, "Failed to setup ucvm.\n");
                    fflush(stderr);
                    exit(-1);
                }
                ucvm_initialized = 1;
			}
			if (ucvm_add_model("sfcvm") != UCVM_CODE_SUCCESS) {
				fprintf(stderr, "Error retrieving SFCVM model.\n");
				fflush(stdout);
			}
		} else {
			fprintf(stderr, "Model %s didn't match any known models, aborting.\n", tok);
			fflush(stderr);
			exit(-4);
		}
		tok = strtok_r(NULL, ",", &save);
	}
 } else {
	//If cvmh_nogtl, use the alternative config file
	if (strcmp(models, "cvmh_nogtl")==0) {
	     if (ucvm_init("/lustre/atlas/proj-shared/geo112/CyberShake/software/UCVM/ucvm-15.10.0/conf/cvmh_no_gtl.conf") != UCVM_CODE_SUCCESS) {
		fprintf(stderr, "Failed to setup ucvm.\n");
		fflush(stderr);
		exit(-1);
		}
	} else {
	     if (ucvm_init(config_filename) != UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Failed to setup ucvm.\n");
                fflush(stderr);
                exit(-1);
                }
	}
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
	} else if (strcmp(models, "cvmsi")==0) {
             if (ucvm_add_model(UCVM_MODEL_CVMSI)!=UCVM_CODE_SUCCESS) {
                     fprintf(stderr, "Error retrieving CVM-S4.26 model.\n");
                     fflush(stderr);
                     exit(-3);
             }
	} else if (strcmp(models, "cvms5")==0) {
             if (ucvm_add_model("cvms5") != UCVM_CODE_SUCCESS) {
             	fprintf(stderr, "Error retrieving CVM-S5 model.\n");
                fflush(stderr);
                exit(-3);
             }
	} else if (strcmp(models, "scec1d")==0) {
             if (ucvm_add_model(UCVM_MODEL_1D)!=UCVM_CODE_SUCCESS) {
                     fprintf(stderr, "Error retrieving 1D model.\n");
                     fflush(stderr);
                     exit(-3);
             }
	} else if (strcmp(models, "cvmh_nogtl")==0) {
             if (ucvm_add_model(UCVM_MODEL_CVMH)!=UCVM_CODE_SUCCESS) {
                     fprintf(stderr, "Error retrieving CVM-H (no GTL) model.\n");
                     fflush(stderr);
                     exit(-3);
             }
	} else if (strcmp(models, "bbp1d")==0) {
	     if (ucvm_add_model(UCVM_MODEL_BBP1D)!=UCVM_CODE_SUCCESS) {
		     fprintf(stderr, "Error retrieving BBP 1D model.\n");
		     fflush(stderr);
		     exit(-3);
	     }
	} else if (strcmp(models, "cca1d")==0) {
             if (ucvm_add_model(UCVM_MODEL_BBP1D)!=UCVM_CODE_SUCCESS) {
                     fprintf(stderr, "Error retrieving CCA 1D model.\n");
                     fflush(stderr);
                     exit(-3);
             }
	} else if (strcmp(models, "cca")==0) {
		printf("Preparing to load CCA model.\n");
	     if (ucvm_add_model("cca")!=UCVM_CODE_SUCCESS) {
		     fprintf(stderr, "Error retrieving CCA model.\n");
		     fflush(stderr);
		     exit(-3);
	     }
	} else if (strcmp(models, "sfcvm")==0) {
            if (ucvm_add_model("sfcvm") != UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Error retrieving SFCVM model.\n");
                fflush(stdout);
            }
	} else {
 	       fprintf(stderr, "Model %s didn't match any known models, aborting.\n", models);
		   fflush(stderr);
		   exit(-3);
        }
 }

 //Figure out which if any models are including the taper
 char taper_models_list[10][512];
 int num_taper_models = 0;
 if (strlen(ely_taper_models)>0 && strcmp(ely_taper_models, "none")!=0) {  
	//Only supported with awpz format for now
	if (format!=AWP_Z) {
		printf("Using the taper with only some models is only supported when using the AWP_Z format, aborting.\n");
		exit(-4);
	}
    char* save;
    char* tok = strtok_r(ely_taper_models, ",", &save);
    while (tok!=NULL) {
		if (my_id==0) {
			printf("Adding %s to models list.", tok);
			fflush(stdout);
		}
 		strcpy(taper_models_list[num_taper_models], tok);
		num_taper_models++;
		tok = strtok_r(NULL, ",", &save);
	}
  } else {
	if (my_id==0) {
		printf("Model list contains %s.\n", ely_taper_models);
		fflush(stdout);
	}
	if (strcmp(ely_taper_models, "all")==0) {
		strcpy(taper_models_list[num_taper_models], "all");
		num_taper_models = 1;
	} else if (strcmp(ely_taper_models, "none")!=0) {
		strcpy(taper_models_list[num_taper_models], ely_taper_models);
		num_taper_models = 1;
	}
	//If ely_taper_models=="none", don't do anything
  }

 // Query by z-type selected earlier
 if (ucvm_setparam(UCVM_MODEL_PARAM_QUERY_MODE, z_query_mode)!=UCVM_CODE_SUCCESS) {
   fprintf(stderr, "Set query mode by depth failed.\n");
   fflush(stderr);
   exit(-2);
 }

 fprintf(stderr, "[%d] Generating CVM\n", my_id);
 fflush(stderr);

 // output order: x,z,y
 int num_pts = local_np;
 ucvm_point_t* pts = check_malloc(sizeof(ucvm_point_t)*(size_t)num_pts);
 // To conserve memory, only allocate one stripe's worth of ucvm_data_t, but num_pts worth of velocity info
 ucvm_data_t* tmp_props = check_malloc(sizeof(ucvm_data_t)*pts_per_stripe);
 // props - final results after Ely taper processing stored here
 ucvm_data_t* props = NULL;
 ucvm_data_t* combine_props = NULL;
 ucvm_data_t* ely_props = NULL;
 
 

 //put the if statement out here, for easier optimization
 if (format==RWG) {
	//fast x, z, y
	for(s=starting_stripe; s<ending_stripe; s++) {
		//Get z and y values from stripe
		int z_ind = s % nz;
		int y_ind = s / nz;
		//Check to see if we should use const mantle
		float z_value = mdep[z_ind];
		if (const_mantle && z_value>const_mantle_depth) {
			z_value = const_mantle_depth;
		}
		//Check to see if we are at surface and should use a different value
		if (z_ind==0) {
			z_value = surface_cvm_depth;
		}
        //If doing query by elevation, need to invert depths to get elevations
        if (z_query_mode==UCVM_COORD_GEO_ELEV) {
            z_value = -1.0*z_value;
        }
		for(i=0; i<nx; i++) {
			//offset in coordinate list (fast x, y)
			int input_offset = y_ind*nx + i;
			int output_offset = (s-starting_stripe)*nx+i;
			pts[output_offset].coord[0] = mlon[input_offset];
                        pts[output_offset].coord[1] = mlat[input_offset];
                        pts[output_offset].coord[2] = z_value;
		}
	}
 } else if (format==AWP_Z) {
	//Fast Z, y, x
	float z_value;
	for(s=starting_stripe; s<ending_stripe; s++) {
		int y_ind = s % ny;
		int x_ind = s / ny;
		for (k=0; k<nz; k++) {
			//offset in coordinate list (fast x, y)
			int input_offset = y_ind*nx + x_ind;
			int output_offset = (s-starting_stripe)*nz + k;
			pts[output_offset].coord[0] = mlon[input_offset];
            pts[output_offset].coord[1] = mlat[input_offset];
			z_value = mdep[k];
			if (k==0) {
				z_value = surface_cvm_depth;
			}
			if (const_mantle && z_value>const_mantle_depth) {
				z_value = const_mantle_depth;
			}
		    //If doing query by elevation, need to invert depths to get elevations
        	if (z_query_mode==UCVM_COORD_GEO_ELEV) {
        	    z_value = -1.0*z_value;
        	}
			pts[output_offset].coord[2] = z_value;
		}
	}
 } else {
	//AWP format
	//output is fast y, x, z
	for(s=starting_stripe; s<ending_stripe; s++) {
		int x_ind = s % nx;
		int z_ind = s / nx;
		float z_value = mdep[z_ind];
        if (const_mantle && z_value>const_mantle_depth) {
                z_value = const_mantle_depth;
        }
		if (z_ind==0) {
			z_value = surface_cvm_depth;
		}
		//If doing query by elevation, need to invert depths to get elevations
		if (z_query_mode==UCVM_COORD_GEO_ELEV) {
			z_value = -1.0*z_value;
		}
		for(j=0; j<ny; j++) {
			int input_offset = j*nx + x_ind;
			int output_offset = (s-starting_stripe)*ny + j;
			pts[output_offset].coord[0] = mlon[input_offset];
            pts[output_offset].coord[1] = mlat[input_offset];
            pts[output_offset].coord[2] = z_value;
		}
	}
 }

 int ely_init = 0;

 for (s=0; s<(ending_stripe-starting_stripe); s++) {
	fprintf(stderr, "[%d] Stripe %d of %d\n", my_id, s+1, (ending_stripe-starting_stripe));
	if (s%10==0) {
		fflush(stderr);
	}
	struct rusage my_rusage;
	getrusage(RUSAGE_SELF, &my_rusage);
	//fprintf(stderr, "[%d] Using %ld kbytes.\n", my_id, my_rusage.ru_maxrss);
	//fflush(stderr);
	//fprintf(stderr, "pts_per_stripe = %d, pts[0]=(%f, %f, %f), @tmp_props=%x\n", pts_per_stripe, (pts+s*pts_per_stripe)->coord[0], (pts+s*pts_per_stripe)->coord[1], (pts+s*pts_per_stripe)->coord[2], tmp_props);
	//fflush(stderr);
	if (ucvm_query(pts_per_stripe, pts+s*pts_per_stripe, tmp_props)!=UCVM_CODE_SUCCESS) {
 		fprintf(stderr, "Query UCVM failed.\n");
    	exit(-2);
	}

	//Consider Ely taper
	if (strcmp(ely_taper, "none")==0) {
		//No taper here
		props = tmp_props;
	} else {
		//Do need to consider taper
		if (format==AWP_Z) {
			//Did this point come from a model we're tapering?
			//Since we're running z-columns, we only have to check the surface point
			int flag = 0;
			if (strcmp(taper_models_list[0], "all")==0) {
				flag = 1;
			} else {
				char pt_label[64];
				//printf("%d) Model ID: %d\n", my_id, tmp_props[0].crust.source);
				//fflush(stdout);
				ucvm_model_label(tmp_props[0].crust.source, pt_label, 64);
				for (i=0; i<num_taper_models; i++) {
					if (strcmp(pt_label, taper_models_list[i])==0) {
						flag = 1;
						//fprintf(stderr, "%d) Applying taper to model name %s.\n", my_id, pt_label);
						//fflush(stderr);
						break;
					}
				}
			}
			if (flag==1) {
				//fprintf(stderr,"Applying taper to column at (%f, %f)\n", pts[s*pts_per_stripe].coord[0], pts[s*pts_per_stripe].coord[1]);
				//fflush(stderr);
				//Need to apply taper to this column
				//Check to see how deep we need to go
				int cutoff_index = -1;
				for (i=0; i<pts_per_stripe; i++) {
					//Changed to fabs() because z-coordinate could be negative if querying by elevation
					if (fabs(pts[s*pts_per_stripe+i].coord[2])>fabs(ely_transition_value)) {
						cutoff_index = i;
						break;
					}
				}
				if (cutoff_index>0) {
					if (ely_init==0) {
						ucvm_modelconf_t conf;
	                    ucvm_elygtl_model_init(UCVM_MAX_MODELS-1, &conf);
	                    ely_init = 1;
					}
				}
				//First, populate ely_props with ely taper data, down to and including transition_depth
				if (ely_props==NULL) {
                    ely_props = check_malloc(sizeof(ucvm_data_t)*cutoff_index);
                }
				memcpy(ely_props, tmp_props, sizeof(ucvm_data_t)*cutoff_index);
				double* original_depths = check_malloc(sizeof(double)*cutoff_index);
				//Change the z-depth to the transition depth.  This is because we need to re-query the crustal model at the transition depth for use with the taper.
				for (i=0; i<cutoff_index; i++) {
					original_depths[i] = pts[s*pts_per_stripe+i].coord[2];
                    pts[s*pts_per_stripe+i].coord[2] = ely_transition_value;
                }
				//Query crustal model again, now that we have the transition depth set
                if (ucvm_query(cutoff_index, pts+s*pts_per_stripe, ely_props)!=UCVM_CODE_SUCCESS) {
					printf("Error querying UCVM.\n");
                    exit(-3);
            	}
                //Set domain to UCVM_DOMAIN_INTERP for elygtl query, and depth back for interpolation
                for (i=0; i<cutoff_index; i++) {
                    ely_props[i].domain = UCVM_DOMAIN_INTERP;
                    ely_props[i].depth = original_depths[i];
                }
				free(original_depths);
                //Set the floors to something very low, so that we can handle the vp/vs ratio scaling in here instead
                double ucvm_taper_floors[3] = {100.0, 500.0, 500.0};
                ucvm_setfloor(ucvm_taper_floors);

				ucvm_elygtl_model_query(UCVM_MAX_MODELS-1, z_query_mode, cutoff_index, pts+s*pts_per_stripe, ely_props);
				//Now run interpolator

                for (i=0; i<cutoff_index; i++) {
                    ucvm_interp_ely(0.0, ely_transition_value, z_query_mode, pts+s*pts_per_stripe+i, ely_props+i);
                }
				//Decide which to use based on method
                if (strcmp(ely_taper, "all")==0) {
                    //Always use the ely_taper value
                    props = ely_props;
                } else if (strcmp(ely_taper, "ifless")==0) {
                    //Select the one which is smallest.  We already know this stripe is above or at the transition depth.
                    if (combine_props==NULL) {
                        combine_props = check_malloc(sizeof(ucvm_data_t)*pts_per_stripe);
                    }
                    for (i=0; i<cutoff_index; i++) {
                        if (tmp_props[i].cmb.vs<ely_props[i].cmb.vs) {
                            memcpy(combine_props+i, tmp_props+i, sizeof(ucvm_data_t));
                        } else {
                            memcpy(combine_props+i, ely_props+i, sizeof(ucvm_data_t));
                        }
                    }
					memcpy(combine_props+cutoff_index, tmp_props+cutoff_index, sizeof(ucvm_data_t)*(pts_per_stripe-cutoff_index));
                    props = combine_props;
                } else {
                    fprintf(stderr, "Don't recognize Ely taper type '%s', aborting.\n", ely_taper);
                    exit(4);
                }
			} else {
				//Otherwise, don't need to apply taper
				props = tmp_props;
			}
		} else {
			//With other formats, taper is applied to all models
			//Check depth of this stripe to see if we even need to worry about it
			//Under both RWG and AWP, all points in a stripe have the same Z-depth; can just check the first one
			if (fabs(pts[s*pts_per_stripe].coord[2])<=fabs(ely_transition_value)) {
				if (ely_init==0) {
					ucvm_modelconf_t conf;
					ucvm_elygtl_model_init(UCVM_MAX_MODELS-1, &conf);
					ely_init = 1;
				}
				//First, populate ely_props with ely taper data, down to and including transition_depth
				if (ely_props==NULL) {
					ely_props = check_malloc(sizeof(ucvm_data_t)*pts_per_stripe);
				}
				memcpy(ely_props, tmp_props, sizeof(ucvm_data_t)*pts_per_stripe);
				double original_depth = pts[s*pts_per_stripe].coord[2];
				//Change the z-depth to the transition depth.  This is because we need to re-query the crustal model at the transition depth for use with the taper.
				for (i=0; i<pts_per_stripe; i++) {
					pts[s*pts_per_stripe+i].coord[2] = ely_transition_value;
				}
				//Query crustal model again, now that we have the transition depth set
				if (ucvm_query(pts_per_stripe, pts+s*pts_per_stripe, ely_props)!=UCVM_CODE_SUCCESS) {
					printf("Error querying UCVM.\n");
					exit(-3);
			}
				//Set domain to UCVM_DOMAIN_INTERP for elygtl query, and depth back for interpolation
				for (i=0; i<pts_per_stripe; i++) {
					ely_props[i].domain = UCVM_DOMAIN_INTERP;
					ely_props[i].depth = original_depth;
				}
				//printf("Before elygtl query, pt (%f, %f, %f), vs30=%lf, crust.vs=%f\n", pts[s*pts_per_stripe].coord[0], pts[s*pts_per_stripe].coord[1], pts[s*pts_per_stripe].coord[2], ely_props[0].vs30, ely_props[0].crust.vs);
				ucvm_elygtl_model_query(UCVM_MAX_MODELS-1, z_query_mode, pts_per_stripe, pts+s*pts_per_stripe, ely_props);
				//printf("pt (%f, %f, %f), vs30=%lf, GTL Vs=%f, Vs=%f\n", pts[s*pts_per_stripe].coord[0], pts[s*pts_per_stripe].coord[1], pts[s*pts_per_stripe].coord[2], ely_props[0].vs30, ely_props[0].gtl.vs, ely_props[0].cmb.vs);
				//Now run interpolator
				for (i=0; i<pts_per_stripe; i++) {
	                ucvm_interp_ely(0.0, ely_transition_value, z_query_mode, pts+s*pts_per_stripe+i, ely_props+i);
	            }
				//printf("After interpolator, pt (%f, %f, %f), vs30=%lf, GTL Vs=%f, Vs=%f\n", pts[s*pts_per_stripe].coord[0], pts[s*pts_per_stripe].coord[1], pts[s*pts_per_stripe].coord[2], ely_props[0].vs30, ely_props[0].gtl.vs, ely_props[0].cmb.vs);
	
				//Decide which to use based on method
				if (strcmp(ely_taper, "all")==0) {
					//Always use the ely_taper value
					props = ely_props;
				} else if (strcmp(ely_taper, "ifless")==0) {
					//Select the one which is smallest.  We already know this stripe is above or at the transition depth.
					if (combine_props==NULL) {
						combine_props = check_malloc(sizeof(ucvm_data_t)*pts_per_stripe);
					}
					for (i=0; i<pts_per_stripe; i++) {
						if (tmp_props[i].cmb.vs<ely_props[i].cmb.vs) {
							memcpy(combine_props+i, tmp_props+i, sizeof(ucvm_data_t));
						} else {
							memcpy(combine_props+i, ely_props+i, sizeof(ucvm_data_t));
						}
					}
					props = combine_props;
				} else {
					fprintf(stderr, "Don't recognize Ely taper type '%s', aborting.\n", ely_taper);
					exit(4); 
				}
			} else {
				//This stripe is too deep to be affected by the Ely taper
				props = tmp_props;
			}
		}
	}

	 /* perform sanity checks on the material properties */     

 	//fprintf(stderr, "Enforcing min_vp=%f, min_vs=%f, min_rho=%f.\n", min_vp, min_vs, min_rho);

 	for (i=0; i<pts_per_stripe; i++) { 
		//Check for nan, inf
		if (isnan(props[i].cmb.vp) || isnan(props[i].cmb.vs) || isnan(props[i].cmb.rho) ||
		   isinf(props[i].cmb.vp) || isinf(props[i].cmb.vs) || isinf(props[i].cmb.rho)) {
			 fprintf(stderr, "NaN/Inf detected at (%f, %f, %f)\n", 
			 pts[s*pts_per_stripe+i].coord[0], pts[s*pts_per_stripe+i].coord[1], pts[s*pts_per_stripe+i].coord[2]);
	 		 fflush(stderr);
	 		 exit(-1);
		}

		//Check for negative values
	       if ((props[i].cmb.vp < 0.0) || (props[i].cmb.vs < 0.0) || (props[i].cmb.rho < 0.0)) {
		 fprintf(stderr, "Negative vals detected at (%f, %f, %f)\n", 
			 pts[s*pts_per_stripe+i].coord[0], pts[s*pts_per_stripe+i].coord[1], pts[s*pts_per_stripe+i].coord[2]);
		 fprintf(stderr, "vp=%f, vs=%f, rho=%lf\n", 
			 props[i].cmb.vp, props[i].cmb.vs, props[i].cmb.rho);
		 fflush(stderr);
		 exit(-1);
	       }

		//Check for min Vp, Vs, Rho
		if (props[i].cmb.vs<min_vs) {
			//If we're doing z-slices, we can use the updated Vp/Vs ratio approach to avoid too-large Vp values, developed for Study 24.8
			float vpvs_ratio=props[i].cmb.vp/props[i].cmb.vs;
			if (format==AWP_Z) {
				//If this is the surface point, we'll calculate the ratio differently
				//if (pts[s*pts_per_stripe+i].coord[2]==0.0 || pts[s*pts_per_stripe+i].coord[2]==surface_cvm_depth) {
				if (props[i].depth==0.0 || props[i].depth==surface_cvm_depth) {
					//Calculate Vp/Vs ratio at the next point, which is one grid point down
					//printf("%d) Point index %d (%lf, %lf, %lf) has Vp=%lf, Vs=%lf, ratio=%f.  The next point is (%lf, %lf, %lf) with Vp=%lf, Vs=%lf, ratio=%f.\n", my_id, i, pts[s*pts_per_stripe+i].coord[0], pts[s*pts_per_stripe+i].coord[1], pts[s*pts_per_stripe+i].coord[2], props[i].cmb.vp, props[i].cmb.vs, vpvs_ratio, pts[s*pts_per_stripe+i+1].coord[0], pts[s*pts_per_stripe+i+1].coord[1], pts[s*pts_per_stripe+i+1].coord[2], props[i+1].cmb.vp, props[i+1].cmb.vs, props[i+1].cmb.vp/props[i+1].cmb.vs);
					float next_vpvs_ratio = props[i+1].cmb.vp/props[i+1].cmb.vs;
					if (next_vpvs_ratio>4.0) {
						next_vpvs_ratio = 4.0;
					}
					vpvs_ratio = next_vpvs_ratio;
				}
			}
	        props[i].cmb.vs = min_vs;
			props[i].cmb.vp = min_vs*vpvs_ratio;
        }
		if (props[i].cmb.vp<min_vp) {
			props[i].cmb.vp = min_vp;
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
		        vp_buf[s*pts_per_stripe+i] = meter2km*props[i].cmb.vp;
	                vs_buf[s*pts_per_stripe+i] = meter2km*props[i].cmb.vs;
	        	rho_buf[s*pts_per_stripe+i] = meter2km*props[i].cmb.rho;
		} else if (format==AWP || format==AWP_Z) {
			buf[3*(s*pts_per_stripe+i)] = props[i].cmb.vp;
	                buf[3*(s*pts_per_stripe+i)+1] = props[i].cmb.vs;
	                buf[3*(s*pts_per_stripe+i)+2] = props[i].cmb.rho;
		}
	
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
		//fprintf(stderr, "%d) writing at offset %ld\n", my_id, offset);
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
 } else if (format==AWP_Z) {
	MPI_Info info;
	MPI_File fp_out;
	MPI_Offset offset;
	MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_out);
	//Stripes are fast z, y, x, but output is fast y, x, z
	//Copy output into fast y, x, z array with constant Z
	//Multiple writes, since not contiguous in the file
	float* tmp_buf = check_malloc(3*sizeof(float)*(ending_stripe-starting_stripe));
	for (k=0; k<nz; k++) {
		//Copy values into tmp_buf
		for(j=0; j<ending_stripe-starting_stripe; j++) {
			tmp_buf[3*j] = buf[3*(j*nz + k)];
			tmp_buf[3*j + 1] = buf[3*(j*nz + k) + 1];
			tmp_buf[3*j + 2] = buf[3*(j*nz + k) + 2];
		}
		offset = ((long)k * nx * ny + starting_stripe) * 3 * sizeof(float);
		printf("%d) Offset is at %lld bytes.\n", my_id, offset);
		int rc = MPI_File_write_at_all(fp_out, offset, tmp_buf, 3*(ending_stripe-starting_stripe), MPI_FLOAT, &status);
		if (rc!=MPI_SUCCESS) {
	          char error_string[256];
              int len_err_string;
              MPI_Error_string(rc, error_string, &len_err_string);
              fprintf(stderr, "%d) Error writing to file %s: %s\n", my_id, outfile, error_string);
              mpi_exit(3);
        } else {
        	fprintf(stderr, "[%d] File write %d of %d completed.\n", my_id, i, nz);
    	}
	}
    MPI_File_close(&fp_out);
	free(tmp_buf);
 } else if (format==AWP) {
	fprintf(stderr, "[%d] Writing to file %s.\n", my_id, outfile);
	MPI_Info info;
        MPI_File fp_out;
        MPI_Offset offset;
        MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_out);
	//output is fast y, x, z, all 3 values
	offset = (long)starting_stripe * ny * 3 * sizeof(float);
	printf("%d) Offset is %lld bytes.\n", my_id, offset);
	int rc = MPI_File_write_at_all(fp_out, offset, buf, 3*local_np, MPI_FLOAT, &status);
        if (rc!=MPI_SUCCESS) {
        	char error_string[256];
                int len_err_string;
                MPI_Error_string(rc, error_string, &len_err_string);
                fprintf(stderr, "%d) Error writing to file %s: %s\n", my_id, outfile, error_string);
                mpi_exit(3);
        } else {
		fprintf(stderr, "[%d] Writing to file completed.\n", my_id);
	}
	MPI_File_close(&fp_out);
 }

free(mlat);
free(mlon);
free(mdep);

free(pts);
//Don't need to free props, it always points somewhere else
free(tmp_props);
if (ely_props!=NULL) {
	free(ely_props);
}
if (combine_props!=NULL) {
	free(combine_props);
}

free(vs_buf);
free(vp_buf);
free(rho_buf);
free(buf);

fprintf(stderr, "[%d] Completed.\n", my_id);
fflush(stderr);
fflush(stdout);

mpi_exit(0);
}
