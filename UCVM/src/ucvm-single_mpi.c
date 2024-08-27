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

//char config_filename[] = "/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0_01302019/conf/ucvm.conf";

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
		                if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
                                if (ucvm_init("/lustre/atlas/proj-shared/geo112/ucvm_18_5/UCVMC/install/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
                                if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
                                if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
                                if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
                                if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
                                if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
                                if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
                                if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
                                if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
	     if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
 ucvm_data_t* props = check_malloc(sizeof(ucvm_data_t)*pts_per_stripe);

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
		for(i=0; i<nx; i++) {
			//offset in coordinate list (fast x, y)
			int input_offset = y_ind*nx + i;
			int output_offset = (s-starting_stripe)*nx+i;
			pts[output_offset].coord[0] = mlon[input_offset];
                        pts[output_offset].coord[1] = mlat[input_offset];
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
		for(j=0; j<ny; j++) {
			int input_offset = j*nx + x_ind;
			int output_offset = (s-starting_stripe)*ny + j;
			pts[output_offset].coord[0] = mlon[input_offset];
                        pts[output_offset].coord[1] = mlat[input_offset];
                        pts[output_offset].coord[2] = z_value;
		}
	}
 }
 
 for (s=0; s<(ending_stripe-starting_stripe); s++) {
	//if (ucvm_query(local_np, pts, props)!=UCVM_CODE_SUCCESS) {
	fprintf(stderr, "[%d] Stripe %d of %d\n", my_id, s+1, (ending_stripe-starting_stripe));
	struct rusage my_rusage;
	getrusage(RUSAGE_SELF, &my_rusage);
	fprintf(stderr, "[%d] Using %ld kbytes.\n", my_id, my_rusage.ru_maxrss);
	if (ucvm_query(pts_per_stripe, pts+s*pts_per_stripe, props)!=UCVM_CODE_SUCCESS) {
 		fprintf(stderr, "Query UCVM failed.\n");
        	exit(-2);
	}
	/*if (my_id==0) {
		//Print the first point
		printf("%d) Point 0 i(%f, %f, %f) has properties vp=%lf, vs=%lf, rho=%lf, crust vp=%lf, crust vs=%lf, crust rho=%lf\n", my_id, pts[0].coord[0], pts[0].coord[1], pts[0].coord[2], props[0].cmb.vp, props[0].cmb.vs, props[0].cmb.rho, props[0].crust.vp, props[0].crust.vs, props[0].crust.rho);
	}*/

	 /* perform sanity checks on the material properties */     

 	//fprintf(stderr, "Enforcing min_vp=%f, min_vs=%f, min_rho=%f.\n", min_vp, min_vs, min_rho);

 	for (i=0; i<pts_per_stripe; i++) { 
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
		        vp_buf[s*pts_per_stripe+i] = meter2km*props[i].cmb.vp;
	                vs_buf[s*pts_per_stripe+i] = meter2km*props[i].cmb.vs;
	        	rho_buf[s*pts_per_stripe+i] = meter2km*props[i].cmb.rho;
	        } else if (format==AWP) {
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
free(props);

free(vs_buf);
free(vp_buf);
free(rho_buf);
free(buf);

fprintf(stderr, "[%d] Completed.\n", my_id);
fflush(stderr);
fflush(stdout);

mpi_exit(0);
}
