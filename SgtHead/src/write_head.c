#include "include.h"
#include "function.h"

/*
* This program generates a separate header file containing the data needed for emod3d
* to index a set of SGTs generated using AWP-ODC-SGT.
*
*/

void write_header(char* modelbox, char* coordfile, char* fdloc, char* gridout, float spacing, int nt, float dt, int sgt_tinc, char* component, float moment, float flo, char* awp_media_file, char* outfile, int mu_corr_flag);

int main(int argc, char** argv) {
	if (argc<14) {
		printf("Usage: %s <modelbox_file> <coordfile> <fdloc_file> <gridout_file> <spacing> <nt> <dt> <time decimation> <component> <moment> <max_freq> <awp_media_file> <outfile> [-c]\n", argv[0]);
		exit(1);
	}
	char* modelbox = argv[1];
	char* coordfile = argv[2];
	char* fdloc = argv[3];
	char* gridout = argv[4];
	float spacing = atof(argv[5]);
	int nt = atoi(argv[6]);
	float dt = atof(argv[7]);
	int sgt_tinc = atoi(argv[8]);
	char* component = argv[9];
	float moment = atof(argv[10]);
	float flo = atof(argv[11]);
	char* awp_media_file = argv[12];
	char* outfile = argv[13];
	int mu_corr_flag = 0;
	if (argc==15) {
		if (strcmp(argv[14], "-c")==0) {
			printf("Using unrelaxed/corrected mu.\n");
			mu_corr_flag = 1;
		}
	}
	write_header(modelbox, coordfile, fdloc, gridout, spacing, nt, dt, sgt_tinc, component, moment, flo, awp_media_file, outfile, mu_corr_flag);
}

/* Header copies format of headers in RWG SGTs
*
* 1  * sgtmast - projection, modelbox, np, nt
* np * sgtindex - point index, grid location, spacing
* np * sgtheader - index, modelbox, nt, dt,  site loc, params
* Good potential for optimizing by normalizing header information.
*/

void write_header(char* modelbox, char* coordfile, char* fdloc, char* gridout, float spacing, int nt, float dt, int sgt_tinc, char* component, float moment, float flo, char* awp_media_file, char* outfile, int mu_corr_flag) {
	FILE* modelbox_fp, *coordfile_fp, *fdloc_fp, *gridout_fp;
	int fd_out, media_fp;
	struct sgtmaster sgtmast;
	struct sgtindex* sgtindx;
	struct sgtheader* sgthead;
	//used to construct geocord matrices
	struct runparams rpars;
	char buffer[1024];
	int i, j, k;
	float lonsgt, latsgt, depsgt;
	int xsrc, ysrc;
	float xx, yy;
	float erad = ERAD;
	double g0;
	float xsgt, ysgt, zsgt;
	float* vel_data;
	float lam, mu;
	int velocity_index;
	float convert_factor = 1.0e+10;
	float scale_factor = 0.001; //to compensate for m/s in AWP velocity file

	//Set magic values for 0.5 Hz
	float fh = 25.0;
	float fp = 0.5;
	float qs = 25.0;

	//set defaults
	rpars.geoproj = 1;
	rpars.xshift = -1.0e+15;
	rpars.yshift = -1.0e+15;
	
	modelbox_fp = fopfile(modelbox, "r");
	//5th line has mlon, mlat, mrot information
	fgets(buffer, 1024, modelbox_fp);
        fgets(buffer, 1024, modelbox_fp);
        fgets(buffer, 1024, modelbox_fp);
        fgets(buffer, 1024, modelbox_fp);
        fgets(buffer, 1024, modelbox_fp);
	char* tok = strtok(buffer, "= ");
	tok = strtok(NULL, "= ");
	rpars.modellon = atof(tok);
        tok = strtok(NULL, "= ");
        tok = strtok(NULL, "= ");
        rpars.modellat = atof(tok);
        tok = strtok(NULL, "= ");
        tok = strtok(NULL, "= ");
	rpars.modelrot = atof(tok);
	fclose(modelbox_fp);

	//Get source location
	fdloc_fp = fopfile(fdloc, "r");
	fscanf(fdloc_fp, "%d %d", &xsrc, &ysrc);
	fclose(fdloc_fp);

	//Get nx, ny, nz
	gridout_fp = fopfile(gridout, "r");
	fgets(buffer, 1024, gridout_fp);
	fscanf(gridout_fp, "nx=%d", &(rpars.nx));
	for (i=0; i<rpars.nx+2; i++) {
		fgets(buffer, 1024, gridout_fp);
	}
	fscanf(gridout_fp, "ny=%d", &(rpars.globny));
        for (i=0; i<rpars.globny+2; i++) {
                fgets(buffer, 1024, gridout_fp);
        }
	fscanf(gridout_fp, "nz=%d", &(rpars.nz));
	fclose(gridout_fp);

	//Construct rpars to generate geocord matrices
	//printf("nx: %d, ny: %d, nz: %d\n", rpars.nx, rpars.globny, rpars.nz);
	rpars.h = spacing;
	rpars.ny1 = 0;
	rpars.ny2 = rpars.globny-1;
	
	init_modelc(&rpars);
	//printf("yshift: %f\n", rpars.yshift);

        sgtmast.geoproj = 1;
	sgtmast.modellon = rpars.modellon;
	sgtmast.modellat = rpars.modellat;
	sgtmast.modelrot = rpars.modelrot;
        sgtmast.xshift = rpars.xshift;
        sgtmast.yshift = rpars.yshift;
	
	//open velocity file
	media_fp = opfile_ro(awp_media_file);

	coordfile_fp = fopfile(coordfile, "r");
	fgets(buffer, 1024, coordfile_fp);
	while(strncmp(buffer,"#",1) == 0) {
		fgets(buffer, 1024, coordfile_fp);
	}
	//Number of points
	sscanf(buffer,"%d",&(sgtmast.globnp));
	sgtmast.localnp = sgtmast.globnp;
	sgtmast.nt = nt/sgt_tinc;

	//printf("globnp=%d\n", sgtmast.globnp);
	sgtindx = check_malloc(sgtmast.globnp * sizeof(struct sgtindex));
	sgthead = check_malloc(sgtmast.globnp * sizeof(struct sgtheader));

	vel_data = check_malloc(3*rpars.nx*rpars.globny*sizeof(float)*3);	
	struct z_structure** z_list = check_malloc(rpars.nz*sizeof(struct z_structure*));
	for (i=0; i<rpars.nz; i++) {
		z_list[i] = NULL;
	}

	struct z_structure* zst;

	//read in entire coordfile
	for (i=0; i<sgtmast.globnp; i++) {
		//fprintf(stderr, "cordfile point %d\n", i);
		fgets(buffer, 1024, coordfile_fp);
                sscanf(buffer, "%d %d %d %Ld %f %f %f",&sgtindx[i].xsgt,&sgtindx[i].ysgt,&sgtindx[i].zsgt,&sgtindx[i].indx,&lonsgt,&latsgt,&sgthead[i].sgt_dep);
		zst = check_malloc(sizeof(struct z_structure));
		zst->sgtindx = &(sgtindx[i]);
		zst->sgthdr = &(sgthead[i]);
		if (z_list[sgtindx[i].zsgt]!=NULL) {
			zst->next = z_list[sgtindx[i].zsgt];
		} else {
			zst->next = NULL;
		}
		z_list[sgtindx[i].zsgt] = zst;
	}

	fprintf(stderr, "Cordfile complete.\n");
	fflush(stderr);
	//Now, iterate through the z_list
	struct z_structure* cur;
	for (i=0; i<rpars.nz; i++) {
		fprintf(stderr, "Z-slice %d of %d.\n", i, rpars.nz);
		fflush(stderr);
		cur = z_list[i];
		if (cur==NULL) {
			continue;
		}
		//fast y, x, z
		//Seek to read in 3 z-slices, special cases if z==0 or z==rpars.nz-1
		if (cur->sgtindx->zsgt!=0) {
			long loc = lseek(media_fp, (long)3*sizeof(float)*(cur->sgtindx->zsgt-1)*rpars.nx*rpars.globny, SEEK_SET);
			//fprintf(stderr, "Seeking to location %ld.\n", loc);
			if (cur->sgtindx->zsgt==rpars.nz-1) { //last slice
				//Read slice before, this slice
				chunk_reed(media_fp, vel_data, (long)3*rpars.nx*rpars.globny*sizeof(float)*2);
			} else {
				//Read 3 slices
				//fprintf(stderr, "Reading %ld bytes.\n", (long)3*rpars.nx*rpars.globny*sizeof(float)*3);
				chunk_reed(media_fp, vel_data, (long)3*rpars.nx*rpars.globny*sizeof(float)*3);
			}
			//fprintf(stderr, "Done reading.\n");
		} else { //z==0
			lseek(media_fp, (long)3*sizeof(float)*(cur->sgtindx->zsgt)*rpars.nx*rpars.globny, SEEK_SET);
			//read in this slice, slice after
	                chunk_reed(media_fp, vel_data, (long)3*rpars.nx*rpars.globny*sizeof(float)*2);
		}

		while (cur!=NULL) {
			struct sgtindex* sip = cur->sgtindx;
			struct sgtheader* shp = cur->sgthdr;
                        //printf("Processing point (%d, %d %d)\n", sip->xsgt, sip->ysgt, sip->zsgt);
			fflush(stdout);
			sip->h = spacing;
	                //sgt header
	                shp->indx = sip->indx;
	                shp->geoproj = sgtmast.geoproj;
	                shp->modellon = sgtmast.modellon;
	                shp->modellat = sgtmast.modellat;
	                shp->modelrot = sgtmast.modelrot;
	                shp->xshift = sgtmast.xshift;
	                shp->yshift = sgtmast.yshift;
	                shp->xazim = 90.0 + sgtmast.modelrot;

	                shp->nt = sgtmast.nt;
	                shp->dt = dt*sgt_tinc;
			shp->h = spacing;

	                shp->xsrc = xsrc;
                	shp->ysrc = ysrc;
                	shp->zsrc = 1;

                	shp->xsgt = sip->xsgt;
	              	shp->ysgt = sip->ysgt;
            		shp->zsgt = sip->zsgt;

	                shp->xmom = 0.0;
	                shp->ymom = 0.0;
                	shp->zmom = 0.0;

	                if (strcmp(component,"x")==0) {
	                        shp->xmom = moment;
	                } else if (strcmp(component,"y")==0) {
	                        shp->ymom = moment;
	                } else if (strcmp(component,"z")==0) {
	                        shp->zmom = moment;
	                }

			//Note that flo is now the frequency the source was filtered at
	                shp->tst = -1.0/flo;
        
	                xx = shp->xsrc*spacing;
        	        yy = shp->ysrc*spacing;
        
	                gcproj(&xx,&yy,&(shp->src_lon),&(shp->src_lat),&rpars.erad,&rpars.g0,&rpars.b0,rpars.amat,rpars.ainv,0);
	                shp->src_dep = (shp->zsrc - 1)*spacing;
        
	                xx = shp->xsgt*spacing;
	                yy = shp->ysgt*spacing;

	                gcproj(&xx,&yy,&(shp->sgt_lon),&(shp->sgt_lat),&rpars.erad,&rpars.g0,&rpars.b0,rpars.amat,rpars.ainv,0);

	                xsgt = shp->xsgt - shp->xsrc;
        	        ysgt = shp->ysgt - shp->ysrc;
        	        zsgt = shp->zsgt - shp->zsrc;
	                shp->cdist = (spacing)*sqrt(xsgt*xsgt + ysgt*ysgt + zsgt*zsgt);

                	//To figure out lam, need to read in velocity model
			int vp_index, vs_index, rho_index, rho_x_plus_one, rho_y_plus_one, rho_z_plus_one;
			//if z!=0, want to read the first plane (since we do a z=z-1 substitution
			//if z==0, still want to read the first plane, since the first plane is z==0
			vp_index = 3*(shp->xsgt*rpars.globny + shp->ysgt);
			vs_index = vp_index + 1;
			if (shp->ysgt==0) {
				//take all the values from y=1 instead, because of copy in genmodel.c, line 3324.
                                rho_index = 3*(shp->xsgt*rpars.globny + shp->ysgt + 1) + 2;
                                rho_x_plus_one = 3*((shp->xsgt+1)*rpars.globny + shp->ysgt + 1) + 2;
                                rho_y_plus_one = 3*(shp->xsgt*rpars.globny + shp->ysgt + 1) + 2;
                                rho_z_plus_one = 3*(rpars.nx*rpars.globny + shp->xsgt*rpars.globny + shp->ysgt + 1) + 2;
			} else {
				rho_index = vp_index + 2;
                                rho_x_plus_one = 3*((shp->xsgt+1)*rpars.globny + shp->ysgt) + 2;
                                rho_y_plus_one = 3*(shp->xsgt*rpars.globny + shp->ysgt + 1) + 2;
                                rho_z_plus_one = 3*(rpars.nx*rpars.globny + shp->xsgt*rpars.globny + shp->ysgt) + 2;
			}

			/*if (shp->xsgt==0 && shp->ysgt==0 && shp->zsgt==5) {
				printf("offsets are vp=%d, vs=%d, rho=%d, rhox+1=%d, rhoy+1=%d, rhoz+1=%d\n", vp_index, vs_index, rho_index, rho_x_plus_one, rho_y_plus_one, rho_z_plus_one);
				printf("Vp=%f, Vs=%f, rho=%f\n", vel_data[vp_index]*scale_factor, vel_data[vs_index]*scale_factor, vel_data[rho_index]*scale_factor);
			}*/

			if (mu_corr_flag==1) {
				float vs_corr = vel_data[vs_index]*scale_factor*(1+(log(fh/fp)/(qs*3.14159)));
				mu = vel_data[rho_index]*scale_factor*vs_corr*vs_corr;
			} else {
		                mu = vel_data[vs_index]*scale_factor*vel_data[vs_index]*scale_factor*vel_data[rho_index]*scale_factor;
			}
			lam = vel_data[vp_index]*scale_factor*vel_data[vp_index]*scale_factor*vel_data[rho_index]*scale_factor - 2*mu;

	                shp->mu = convert_factor*mu;
                	shp->lam = convert_factor*lam;
                	float rho_zero, rho_x, rho_y, rho_z, bx, by, bz;
                	rho_zero = vel_data[rho_index]*scale_factor;
	                if (shp->xsgt==rpars.nx-1 || shp->zsgt==rpars.nz-1) {
	                        shp->rho = rho_zero;
	                } else {
        	                rho_x = vel_data[rho_x_plus_one]*scale_factor;
        	                rho_y = vel_data[rho_y_plus_one]*scale_factor;
        	                rho_z = vel_data[rho_z_plus_one]*scale_factor;
        	                bx = 2.0/(rho_zero+rho_x);
        	                by = 2.0/(rho_zero+rho_y);
        	                bz = 2.0/(rho_zero+rho_z);
        	                shp->rho = 3.0/(bx+by+bz);
                	}
                        /*if (shp->xsgt==0 && shp->ysgt==0 && shp->zsgt==5) {
				printf("rho_zero=%f, rho_x=%f, rho_y=%f, rho_z=%f\n", rho_zero, rho_x, rho_y, rho_z);
			}*/
			struct z_structure* next = cur->next;
			free(cur);
			cur = next;
		}
	}


/*
	for (i=0; i<sgtmast.globnp; i++) {
		fgets(buffer, 1024, coordfile_fp);
                sscanf(buffer, "%d %d %d %Ld %f %f %f",&sgtindx[i].xsgt,&sgtindx[i].ysgt,&sgtindx[i].zsgt,&sgtindx[i].indx,&lonsgt,&latsgt,&depsgt);
		if ((i+1)%100==0) {
			printf("%d of %d\n", i+1, sgtmast.globnp);
		}
		if (sgtindx[i].ysgt!=old_y) {
			//read in velocity chunk
			int y_diff = sgtindx[i].ysgt - old_y;
			size_t pos = lseek(vp_fp, (y_diff-1)*rpars.nx*rpars.nz*sizeof(float), SEEK_CUR);
			//printf("file position: %ld\n", pos);
                        lseek(vs_fp, (y_diff-1)*rpars.nx*rpars.nz*sizeof(float), SEEK_CUR);
			if (old_y==rpars.globny-1) {
				//only read 1 y-slice, so only skip back 1 y-slice
				lseek(rho_fp, (y_diff-1)*rpars.nx*rpars.nz*sizeof(float), SEEK_CUR);
			} else {
	                        lseek(rho_fp, (y_diff-2)*rpars.nx*rpars.nz*sizeof(float), SEEK_CUR);
			}
		
			reed(vp_fp, vp, rpars.nx*rpars.nz*sizeof(float));
                	reed(vs_fp, vs, rpars.nx*rpars.nz*sizeof(float));
			if (sgtindx[i].ysgt==rpars.globny-1) {
				//Last y point, so only read 1 in
	                        reed(rho_fp, rho, rpars.nx*rpars.nz*sizeof(float));
			} else {
	                	reed(rho_fp, rho, 2*rpars.nx*rpars.nz*sizeof(float));
			}
			old_y = sgtindx[i].ysgt;
		}
		sgtindx[i].h = spacing;
		//sgt header
  	        sgthead[i].indx = sgtindx[i].indx;
	       	sgthead[i].geoproj = sgtmast.geoproj;
        	sgthead[i].modellon = sgtmast.modellon;
         	sgthead[i].modellat = sgtmast.modellat;
         	sgthead[i].modelrot = sgtmast.modelrot;
         	sgthead[i].xshift = sgtmast.xshift;
         	sgthead[i].yshift = sgtmast.yshift;
         	sgthead[i].xazim = 90.0 + sgtmast.modelrot;

         	sgthead[i].nt = sgtmast.nt;
         	sgthead[i].dt = dt*sgt_tinc;
         	sgthead[i].h = spacing;

         	sgthead[i].xsrc = xsrc;
	        sgthead[i].ysrc = ysrc;
         	sgthead[i].zsrc = 1;

         	sgthead[i].xsgt = sgtindx[i].xsgt;
         	sgthead[i].ysgt = sgtindx[i].ysgt;
         	sgthead[i].zsgt = sgtindx[i].zsgt;

       	  	sgthead[i].xmom = 0.0;
         	sgthead[i].ymom = 0.0;
         	sgthead[i].zmom = 0.0;
		if (strcmp(component,"x")==0) {
			sgthead[i].xmom = moment;
		} else if (strcmp(component,"y")==0) {
			sgthead[i].ymom = moment;
		} else if (strcmp(component,"z")==0) {
			sgthead[i].zmom = moment;
		}

       	 	sgthead[i].tst = -1.0/flo;
	
	        xx = sgthead[i].xsrc*spacing;
	        yy = sgthead[i].ysrc*spacing;
	
	        gcproj(&xx,&yy,&(sgthead[i].src_lon),&(sgthead[i].src_lat),&rpars.erad,&rpars.g0,&rpars.b0,rpars.amat,rpars.ainv,0);
	        sgthead[i].src_dep = (sgthead[i].zsrc - 1)*spacing;
	
	        xx = sgthead[i].xsgt*spacing;
		yy = sgthead[i].ysgt*spacing;

         	gcproj(&xx,&yy,&(sgthead[i].sgt_lon),&(sgthead[i].sgt_lat),&rpars.erad,&rpars.g0,&rpars.b0,rpars.amat,rpars.ainv,0);
	        sgthead[i].sgt_dep = depsgt;

	        xsgt = sgthead[i].xsgt - sgthead[i].xsrc;
       	 	ysgt = sgthead[i].ysgt - sgthead[i].ysrc;
       	 	zsgt = sgthead[i].zsgt - sgthead[i].zsrc;
	        sgthead[i].cdist = (spacing)*sqrt(xsgt*xsgt + ysgt*ysgt + zsgt*zsgt);

		//To figure out lam, need to read in velocity model and do some pointer reassignment
		velocity_index = sgthead[i].zsgt*rpars.nx + sgthead[i].xsgt;
                int x_plus_one_index, y_plus_one_index, z_plus_one_index;
                x_plus_one_index = sgthead[i].zsgt*rpars.nx + sgthead[i].xsgt + 1;
                y_plus_one_index = rpars.nz*rpars.nx + sgthead[i].zsgt*rpars.nx + sgthead[i].xsgt;
                z_plus_one_index = (sgthead[i].zsgt+1)*rpars.nx + sgthead[i].xsgt;
		int v_index = velocity_index;
		if (sgthead[i].zsgt!=0) {
			velocity_index = (sgthead[i].zsgt-1)*rpars.nx + sgthead[i].xsgt;
		        x_plus_one_index = (sgthead[i].zsgt-1)*rpars.nx + sgthead[i].xsgt + 1;
	                y_plus_one_index = rpars.nz*rpars.nx + (sgthead[i].zsgt-1)*rpars.nx + sgthead[i].xsgt;
        	        z_plus_one_index = (sgthead[i].zsgt-1+1)*rpars.nx + sgthead[i].xsgt;

		}

		mu = vs[velocity_index]*vs[velocity_index]*rho[velocity_index];
		lam = vp[velocity_index]*vp[velocity_index]*rho[velocity_index] - 2*mu;

		sgthead[i].mu = convert_factor*mu;
		sgthead[i].lam = convert_factor*lam;
		//sgthead[i].rho = rho[velocity_index];
		float rho_zero, rho_x, rho_y, rho_z, bx, by, bz;
//		x_plus_one_index = sgthead[i].zsgt*rpars.nx + sgthead[i].xsgt + 1;
//		y_plus_one_index = rpars.nz*rpars.nx + sgthead[i].zsgt*rpars.nx + sgthead[i].xsgt;
//		z_plus_one_index = (sgthead[i].zsgt+1)*rpars.nx + sgthead[i].xsgt;
                rho_zero = rho[velocity_index];
		if (sgthead[i].xsgt==rpars.nx-1 || sgthead[i].zsgt==rpars.nz-1) {
			sgthead[i].rho = rho_zero;
		} else {
			rho_x = rho[x_plus_one_index];
			rho_y = rho[y_plus_one_index];
			rho_z = rho[z_plus_one_index];
			bx = 2.0/(rho_zero+rho_x);
			by = 2.0/(rho_zero+rho_y);
			bz = 2.0/(rho_zero+rho_z);
			sgthead[i].rho = 3.0/(bx+by+bz);
		}
	}
	close(vp_fp);
	close(vs_fp);
	close(rho_fp);
	*/
	close(media_fp);
	fclose(coordfile_fp);
	fd_out = cropfile_rw(outfile);
	//printf("Size of sgtmaster: %ld\n", sizeof(struct sgtmaster));
	//printf("Size of sgtindex: %ld\n", sizeof(struct sgtindex));
	//printf("Size of sgtheader: %ld\n", sizeof(struct sgtheader));
	rite(fd_out,&sgtmast, sizeof(struct sgtmaster));
	rite(fd_out,sgtindx, sgtmast.globnp*sizeof(struct sgtindex));
	rite(fd_out,sgthead, ((long long)sgtmast.globnp)*sizeof(struct sgtheader));
	close(fd_out);

	//free(vp);
	//free(vs);
	//free(rho);
	free(vel_data);
	free(z_list);
	free(sgthead);
	free(sgtindx);
}
