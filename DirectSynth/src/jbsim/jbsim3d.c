#include "include.h"
#include "defs.h"
#include "structure.h"
#include "rupgen_api.h"
#include "function.h"

const int X_COMP_FLAG = 1;
const int Y_COMP_FLAG = 2;
const int Z_COMP_FLAG = 4;

float** jbsim3d_synth(float*** seis_return, struct seisheader* header, char stat[], float slon, float slat,
		int ntout, char seis_file[], char rup_geom_file[], struct sgtfileparams* sfp,
		struct sgtparams* sgtparms, struct sgtmaster sgtmast, struct sgtindex* sgtindx,
		struct geoprojection geop, long long* indx_master, int num_sgts, long long** sgts_by_handler, int* num_sgts_by_handler,
		int num_sgt_handlers, int num_rup_vars, struct rupture_variation* rup_vars, int my_id)
{
FILE *fopfile(), *fpr;
struct sgtfileparams sgtextract;
struct sgtfileparams sgtfilepar = *sfp;
struct sgtindex eqindx, statindx;
float *sgtbuf;
struct sgtheader *sgthead;

float **gfmech;
float **stf, **seis, **subseis, *se, *sn, *sv;
float rt, scale;
float elon, elat, edep;
float vslip, **space;
float z0, strike, dip, rake;
int fdw, ip, maxmech, nstf, ntsum, maxnt, ig;
int i,j;
float mindt, maxdelta, fweight;

char string[256], outfile[128];
char outdir[256], sgtdir[256], sname[8], rupmodfile[256];

struct standrupformat* srf;
struct srf_planerectangle **prect_ptr;
struct srf_prectsegments **prseg_ptr;
struct srf_allpoints **apnts_ptr;
struct srf_apointvalues **apval_ptr;

srf = check_malloc(sizeof(struct standrupformat) * num_rup_vars);
prect_ptr = check_malloc(sizeof(struct srf_planerectangle) * num_rup_vars);
prseg_ptr = check_malloc(sizeof(struct srf_prectsegments) * num_rup_vars);
apnts_ptr = check_malloc(sizeof(struct srf_allpoints) * num_rup_vars);
apval_ptr = check_malloc(sizeof(struct srf_apointvalues) * num_rup_vars);

struct mechparam mechpar;

int nm, non_exact;

int extract_sgt = 0;

float* tmom = check_malloc(sizeof(float) * num_rup_vars);
for (i=0; i<num_rup_vars; i++) {
	tmom[i] = 0.0;
}
float dtout = -1.0;
float slip_conv = 1.0;  /* input slip in cm on each subfault */
float tstart = 0.0;

float sdep = 0.0;

int* apv_off = check_malloc(sizeof(int)*num_rup_vars);
for (i=0; i<num_rup_vars; i++) {
	apv_off[i] = 0;
}
int inbin = 0;

int intmem = 0;
long memlen;

int ptol;
int print_tol = 25;

int slip_id, hypo_id;

//Default 1 GB
int max_buf_mb = 1*1024;
long long MAX_BUFFER_SIZE;

int timing = 0;
struct timeval tv_start, tv_end;

sname[0] = '\0';

sgtextract.xfile[0] = '\0';
sgtextract.yfile[0] = '\0';
sgtextract.zfile[0] = '\0';

sgtextract.xfp = NULL;
//sgtextract.xfdr = -1;
sgtextract.yfp = NULL;
//sgtextract.yfdr = -1;
sgtextract.zfp = NULL;
//sgtextract.zfdr = -1;

sprintf(outdir,".");
sprintf(sgtdir,".");

getpar("max_buf_mb","d",&max_buf_mb);
MAX_BUFFER_SIZE = 1024*(long long)1024*max_buf_mb;

//We need this to be required, because we need the dt argument to pass to rupgen_genslip
mstpar("dtout","f",&dtout);
getpar("tstart","f",&tstart);
getpar("sname","s",sname);

rg_stats_t stats;
set_memcached_server("127.0.0.1");

char rupture_spacing_string[20];
int rupture_spacing;
sprintf(rupture_spacing_string, "random");
getpar("rupture_spacing", "s", rupture_spacing_string);
if (strcmp(rupture_spacing_string,"random")==0) {
        rupture_spacing = RUPGEN_RANDOM_HYPO;
        printf("Using random spacing.\n");
} else if (strcmp(rupture_spacing_string,"uniform")==0) {
        rupture_spacing = RUPGEN_UNIFORM_HYPO;
        printf("Using uniform spacing.\n");
} else {
        fprintf(stderr, "rupture_spacing argument %s must be one of 'random' or 'uniform', aborting.", rupture_spacing_string);
        exit(5);
}
//Generate all ruptures
for (i=0; i<num_rup_vars; i++) {
	//printf("Generating rupture from file %s, slip %d, hypo %d, spacing %d, dtout %f\n", rup_geom_file, rup_vars[i].slip_id, rup_vars[i].hypo_id, rupture_spacing, dtout);
	rupgen_genslip(rup_geom_file, rup_vars[i].slip_id, rup_vars[i].hypo_id, &stats, &(srf[i]), rupture_spacing, dtout);
}

//write_srf(&srf, "out.srf", 0);

for (i=0; i<num_rup_vars; i++) {
	prect_ptr[i] = &(srf[i]).srf_prect;
	prseg_ptr[i] = prect_ptr[i]->prectseg;
	apnts_ptr[i] = &(srf[i]).srf_apnts;
	apval_ptr[i] = apnts_ptr[i]->apntvals;

	apv_off[i] = 0;
}

int num_comps = 0;
if (sgtfilepar.xfile[0]!='\0') {
	num_comps++;
}
if (sgtfilepar.yfile[0]!='\0') {
        num_comps++;
}
if (sgtfilepar.zfile[0]!='\0') {
        num_comps++;
}


fprintf(stdout,"Constructing synthetic this rupture\n");

long current_offset = 0;
memlen = sizeof(struct sgtmaster) + (sgtmast.globnp)*(sizeof(struct sgtindex) + sizeof(struct sgtheader) + 18*(long long)(sgtmast.nt)*sizeof(float));

fprintf(stdout,"Total memory for SGTs= %.2f Mb\n",memlen*1.0e-06);

//sgthead = (struct sgtheader *) check_malloc ((sgtmast.globnp)*sizeof(struct sgtheader));
long long sgtbufsize = 18*(long long)sgtmast.globnp*sgtmast.nt*sizeof(float);

fprintf(stderr, "MAX_BUFFER_SIZE = %ld\n", MAX_BUFFER_SIZE);
fprintf(stderr, "sgtbufsize = %ld\n", sgtbufsize);

//Need to consider buffer size too
int max_points_per_request = MAX_BUFFER_SIZE/(18/3*4*sgtmast.nt*sizeof(float)+sizeof(struct sgtheader));

int request_from_handler_id = 0;
int starting_index = 0;
int ending_index = 0;
while (num_sgts_by_handler[request_from_handler_id]==0 && request_from_handler_id<num_sgt_handlers) {
	request_from_handler_id++;
}

//if (1) {

sgtbuf = check_malloc((long long)18*max_points_per_request*sgtmast.nt*sizeof(float));
sgthead = check_malloc(max_points_per_request*sizeof(struct sgtheader));

gfmech = NULL;
space = NULL;
seis = NULL;
subseis = NULL;
stf = NULL;
//Which ones were used for 1 particular SRF - will either use for setting up sgthead and sgtbuf pointers, or copying splits over for later
int sgt_used_offset[4];
//Which SRF indices we should process this cycle
//Try using srf[0] as the indicator
int* srf_indices_to_process = check_malloc(srf[0].srf_apnts.np*sizeof(int));
int num_srf_pts_found = 0;
//Track srf points which are found
char* srf_pts_used = check_malloc(srf[0].srf_apnts.np*sizeof(char));
//Keep data for points which are missed in first pass
int* srf_indices_split = 0;
struct sgtheader* sgt_heads_split = 0;
float* sgtbuf_split = 0;
int num_srfs_split = 0;
int num_data_pts_split = 0;
memset(srf_pts_used, 0, srf[0].srf_apnts.np);
while (1) {
	printf("%d) ending_index = %d, num_sgts_by_handler[%d] = %d\n", my_id, ending_index, request_from_handler_id, num_sgts_by_handler[request_from_handler_id]);
	if (ending_index==num_sgts_by_handler[request_from_handler_id]) {
		//Move to next handler
		do {
			request_from_handler_id++;
		} while (num_sgts_by_handler[request_from_handler_id]==0 && request_from_handler_id<num_sgt_handlers);
		if (request_from_handler_id>=num_sgt_handlers) {
			//No more requests, we are done
			break;
		}
		starting_index = 0;
		//We will increment ending_index outside the if statement
		ending_index = 0;
	} else {
		starting_index = ending_index;
	}

	ending_index += max_points_per_request;
	if (ending_index>num_sgts_by_handler[request_from_handler_id]) {
		ending_index = num_sgts_by_handler[request_from_handler_id];
	}

	memset(srf_indices_to_process, 0, srf[0].srf_apnts.np);

	request_sgt(sgthead, sgtbuf, num_comps, request_from_handler_id, sgts_by_handler, starting_index, ending_index, sgtmast.nt, my_id);

	if (debug) {
		struct mallinfo m_info;
		m_info = mallinfo();
		char buf[512];
		sprintf(buf, "non_mmap_alloc=%d, free_blocks=%d, fastbin_free_blocks=%d, mmap_alloc=%d, bytes_mmap_alloc=%d, highwater=%d, fastbin_alloc=%d, total_alloced=%d, total_free=%d, releasable=%d", m_info.arena, m_info.ordblks, m_info.smblks, m_info.hblks, m_info.hblkhd, m_info.usmblks, m_info.fsmblks, m_info.uordblks, m_info.fordblks, m_info.keepcost);
		write_log(buf);
	}

	//read_part_sgt(&sgtfilepar,&sgtmast,sgtindx,sgthead,sgtbuf,starting_pt,ending_pt);

	//Generate list of SRF points whose SGTs we have
	num_srf_pts_found = 0;
	if (debug) write_log("Processing SGTs from the request.");
	for (ip=0; ip<srf[0].srf_apnts.np; ip++) {
	//for(ip=0; ip<1; ip++) {
		//Only check the points we haven't found yet
		if (srf_pts_used[ip]==0) {
			for (i=0; i<4; i++) {
				sgt_used_offset[i] = -1;
			}
			//For each sgt that sgtparms says it needs
			int num_found = 0;
			for (i=0; i<sgtparms[ip].nsgt; i++) {
				//Check that it falls w/i in the indices, otherwise it's not even worth looking for
				if (sgtparms[ip].indx[i]>=sgts_by_handler[request_from_handler_id][starting_index] && sgtparms[ip].indx[i]<=sgts_by_handler[request_from_handler_id][ending_index-1]) {
					//if (sgtparms[ip].indx[i]>=sgtindx[starting_pt].indx && sgtparms[ip].indx[i]<=sgtindx[ending_pt-1].indx) {
					num_found++;
					//can't use master_ip to calculate offset
					//master_ip list has only the points used;
					//might be points in SGT file we don't use
					//int start_index = starting_pt;
					int start_index = starting_index;
					//int end_index = ending_pt;
					int end_index = ending_index;
					int mid_index = (start_index + end_index) / 2;
					while(end_index>=start_index) {
						mid_index = (start_index + end_index) / 2;
						//if (sgtparms[ip].indx[i]>sgtindx[mid_index].indx) {
						if (sgtparms[ip].indx[i]>sgts_by_handler[request_from_handler_id][mid_index]) {
							start_index = mid_index+1;
							//							} else if (sgtparms[ip].indx[i]<sgtindx[mid_index].indx) {
						} else if (sgtparms[ip].indx[i]<sgts_by_handler[request_from_handler_id][mid_index]) {
							end_index = mid_index-1;
						} else {
							break;
						}
					}

					sgt_used_offset[i] = mid_index - starting_index;
					//sgt_used_offset[i] = sgtparms[ip].master_ip[i] - starting_pt;
				}
			}
			if (num_found==sgtparms[ip].nsgt) {
				srf_pts_used[ip] = 1;
				srf_indices_to_process[num_srf_pts_found] = ip;
				num_srf_pts_found++;
				//All the SGTs for this point are here
				//set up sgtparams pointers
				for (i=0; i<sgtparms[ip].nsgt; i++) {
					//fprintf(stderr, "Setting up pointers for SGT with indx %ld\n", sgtparms[ip].indx[i]);
					int offset = sgt_used_offset[i];
					//set up sgthead
					//fprintf(stderr, "sgtparms head offset: %d\n", offset);
					sgtparms[ip].sgt_head_ptrs[i] = sgthead+offset;
					//fprintf(stderr, "xmom from sgtparms[%d] = %f, at location %ld\n", ip, sgtparms[ip].sgt_head_ptrs[i]->xmom, sgtparms[ip].sgt_head_ptrs[i]);
					//set up data
					//fprintf(stderr, "sgtparms buf offset: %ld\n", (long long)18*sgtmast.nt*offset);
					sgtparms[ip].sgtbuf_ptrs[i] = sgtbuf+(long long)18*sgtmast.nt*offset;
					//fprintf(stderr, "sgtbuf location: %ld\n", sgtparms[ip].sgtbuf_ptrs[i]);

				}
			} else if (num_found>0) {
				//We're going to come back to this point later
				//Copy the point and header info into a buffer
				fprintf(stderr, "SRF point %d is split, saving for later.\n", ip);
				int found = 0;
				for (i=0; i<num_srfs_split; i++) {
					if (srf_indices_split[i]==ip) {
						found = 1;
						break;
					}
				}
				if (!found) {
					srf_indices_split = check_realloc(srf_indices_split, sizeof(int)*(num_srfs_split+1));
					srf_indices_split[num_srfs_split] = ip;
					num_srfs_split++;
				}

				sgt_heads_split = check_realloc(sgt_heads_split, sizeof(struct sgtheader)*(long long)(num_data_pts_split+num_found));
				sgtbuf_split = check_realloc(sgtbuf_split, sizeof(float)*(long long)(num_data_pts_split+num_found)*18*sgtmast.nt);

				for (i=0; i<sgtparms[ip].nsgt; i++) {
					if (sgt_used_offset[i]!=-1) {
						//Copy this
						memcpy(sgt_heads_split+num_data_pts_split, sgthead+sgt_used_offset[i], sizeof(struct sgtheader));
						memcpy(sgtbuf_split+num_data_pts_split*18*sgtmast.nt, sgtbuf+(long long)sgt_used_offset[i]*18*sgtmast.nt, 18*sgtmast.nt*sizeof(float));
						num_data_pts_split++;
					}
				}
			}

		}
	}

	if (gfmech==NULL) {
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

		header->dt = dtout;
		header->nt = ntout;

		maxmech = 3;
		mechpar.nmech = 1;
		mechpar.flag[0] = U1FLAG;
		mechpar.flag[1] = 0;
		mechpar.flag[2] = 0;

		fprintf(stderr, "gfmech\n");
		gfmech = check_malloc(sizeof(float*) * num_rup_vars);
		space = check_malloc(sizeof(float*) * num_rup_vars);
		seis = check_malloc(sizeof(float*) * num_rup_vars);
		subseis = check_malloc(sizeof(float*) * num_rup_vars);
		stf = check_malloc(sizeof(float*) * num_rup_vars);

		for (i=0; i<num_rup_vars; i++) {
			gfmech[i] = (float *) check_malloc (maxmech*12*ntsum*sizeof(float));
			space[i] = (float *) check_malloc (2*ntsum*sizeof(float));

			seis[i] = (float *) check_malloc (3*ntout*sizeof(float));
			subseis[i] = (float *) check_malloc (maxmech*3*ntout*sizeof(float));
			stf[i] = (float *) check_malloc (ntout*sizeof(float));

			zapit(seis[i],3*ntout);
			ptol = print_tol;
			tmom[i] = 0.0;
		}
	}

	//fprintf(stderr, "Found %d srf points, processing them.\n", num_srf_pts_found);

	for(ip=0; ip<num_srf_pts_found; ip++) {
		int srf_index = srf_indices_to_process[ip];
		/*if (sgtparms[srf_index].indx[0]!=52705670033) {
			continue;
		}*/
		int j;
		for (i=0; i<num_rup_vars; i++) {
			zapit(subseis[i],maxmech*3*ntout);

			get_srfpars(&(srf[i]),apv_off[i],srf_index,&rt,&vslip,&mechpar.stk,&mechpar.dip,&mechpar.rak,&mechpar);
			scale = slip_conv*apval_ptr[i][srf_index].area;
	                //fprintf(stderr, "Scale = %e, vslip = %e\n", scale, vslip);
			mech_sgt(gfmech[i],&sgtparms[srf_index],ntsum,mechpar,&scale);
			tmom[i] = tmom[i] + vslip*scale;
                   	//fprintf(stderr, "tmom = %e, scale = %e\n", tmom[i], scale);
                   	/*for (j=0; j<ntout; j++) {
                   	     printf("gfmech[%d] = %e\n", j, gfmech[i][j]);
                   	}*/
			sum_sgt(subseis[i],ntout,gfmech[i],&sgtparms[srf_index],ntsum,&rt,&tstart,mechpar);
                        /*for (j=0; j<ntout; j++) {
                                printf("subseis[%d] = %e\n", j, subseis[i][j]);
                        }*/

			srf_stf(&(srf[i]),apv_off[i],srf_index,seis[i],subseis[i],stf[i],ntout,&dtout,mechpar,space[i]);
			/*int j;
			for (j=0; j<ntout; j++) {
				printf("subseis[%d] = %e\n", j, subseis[i][j]);
			}*/
		}
	}

}
//Once we've finished the while loop, iterate over all the points we rolled over
num_srf_pts_found = 0;
//fprintf(stderr, "%d split points to process.\n", num_srfs_split);
for (ip=0; ip<num_srfs_split; ip++) {
	int srf_index = srf_indices_split[ip];
	srf_pts_used[srf_index] = 1;
	for (i=0; i<4; i++) {
		sgt_used_offset[i] = -1;
	}
	//Find these points - can't assume any order in the split buffers, so linear search
	//The split buffers should be very small anyway
	for (i=0; i<sgtparms[srf_index].nsgt; i++) {
		//populate sgt_used_offset
		for (j=0; j<num_data_pts_split; j++) {
			if (sgt_heads_split[j].indx == sgtparms[srf_index].indx[i]) {
				sgt_used_offset[i] = j;
				break;
			}
		}
		if (sgt_used_offset[i]==-1) {
			fprintf(stderr, "Error:  Couldn't find SGT indx %ld in the split list, aborting.\n", sgtparms[srf_index].indx[i]);
			exit(-2);
		}
		//set up pointers in sgtparams
		sgtparms[srf_index].sgt_head_ptrs[i] = sgt_heads_split + sgt_used_offset[i];
		sgtparms[srf_index].sgtbuf_ptrs[i] = sgtbuf_split + sgt_used_offset[i]*18*sgtmast.nt;
	}
	//process this point
	for (i=0; i<num_rup_vars; i++) {
		zapit(subseis[i],maxmech*3*ntout);

		get_srfpars(&(srf[i]),apv_off[i],srf_index,&rt,&vslip,&mechpar.stk,&mechpar.dip,&mechpar.rak,&mechpar);
		scale = slip_conv*apval_ptr[i][srf_index].area;

		mech_sgt(gfmech[i],&sgtparms[srf_index],ntsum,mechpar,&scale);
		tmom[i] = tmom[i] + vslip*scale;

		sum_sgt(subseis[i],ntout,gfmech[i],&sgtparms[srf_index],ntsum,&rt,&tstart,mechpar);
		srf_stf(&(srf[i]),apv_off[i],srf_index,seis[i],subseis[i],stf[i],ntout,&dtout,mechpar,space[i]);
	}
}
int all_points_processed = 1;
/*for (ip=0; ip<srf[0].srf_apnts.np; ip++) {
	if (srf_pts_used[ip]==0) {
		fprintf(stderr, "SRF point %d was never processed.\n", ip);
		all_points_processed = 0;
	}
}
if (!all_points_processed) {
	if (debug) close_log();
	MPI_Finalize();
	exit(-1);
}*/
free(srf_indices_to_process);
free(srf_pts_used);

free(srf_indices_split);
free(sgt_heads_split);
free(sgtbuf_split);

//free(sgtparms);
//free(sgtindx);
free(sgtbuf);
free(sgthead);
//free(indx_master);
for (i=0; i<num_rup_vars; i++) {
	free(gfmech[i]);
	free(space[i]);
	free(subseis[i]);
	free(stf[i]);
}
free(gfmech);
free(space);
free(subseis);
free(stf);

*seis_return = seis;

//Accumulate data for writing in single buffer
char* writing_buffer = check_malloc((sizeof(struct seisheader) + sizeof(float)*3*ntout)*num_rup_vars);
int entry_size = 0;
int offset = 0;

if (sgtfilepar.xfile[0]!='\0') {
	if (sgtfilepar.yfile[0]!='\0') {
		if (sgtfilepar.zfile[0]!='\0') { //x,y,z
			header->comps = (X_COMP_FLAG | Y_COMP_FLAG | Z_COMP_FLAG);
			entry_size = 3*sizeof(float)*ntout;
		} else { //x,y
			header->comps = X_COMP_FLAG | Y_COMP_FLAG;
			entry_size = 2*sizeof(float)*ntout;
			offset = ntout;
			//(*seis_return)[i] = sn;
		}
	}
	else if (sgtfilepar.zfile[0]!='\0') { //x,z
                header->comps = X_COMP_FLAG | Z_COMP_FLAG;
    			entry_size = 2*sizeof(float)*ntout;
	} else { //x
        		header->comps = X_COMP_FLAG;
        		offset = ntout;
        		//(*seis_return)[i] = sn;
	}
} else if (sgtfilepar.yfile[0]!='\0') {
	if (sgtfilepar.zfile[0]!='\0') { //y,z
		printf("Can't output merged Y and Z components.\n");
		exit(2);
	} else { //y
		header->comps = Y_COMP_FLAG;
		entry_size = sizeof(float)*ntout;
		offset = 2*ntout;
		//(*seis_return)[i] = se;
	}
} else if (sgtfilepar.zfile[0]!='\0') { //z
	header->comps = Z_COMP_FLAG;
	entry_size = sizeof(float)*ntout;
}

for (i=0; i<num_rup_vars; i++) {
	header->rup_var_id = rup_vars[i].rup_var_id;
	memcpy(writing_buffer+(sizeof(struct seisheader)+entry_size)*i, header, sizeof(struct seisheader));
	memcpy(writing_buffer+(sizeof(struct seisheader)+entry_size)*i + sizeof(struct seisheader), seis[i]+offset, entry_size);
	(*seis_return)[i] = seis[i]+offset;
	if (debug) {
		char buf[256];
		sprintf(buf, "Copying %d bytes from offset %d to offset %d", entry_size, offset, (sizeof(struct seisheader)+entry_size)*i);
		write_log(buf);
	}
}

/*for (i=0; i<num_rup_vars; i++) {

	sv = seis[i];
	sn = seis[i] + ntout;
	se = seis[i] + 2*ntout;

	if(sname[0] == '\0') {
		strncpy(sname,stat,7);
		sname[7] = '\0';
	}

	header->comps = 0;
	header->rup_var_id = rup_vars[i].rup_var_id;*/

//	char* last_slash = strrchr(seis_file, '/');
//	if (last_slash!=NULL) { //means there was a slash in the filename, might need to create a path
//		char* path = check_malloc(sizeof(char)*strlen(seis_file)+1);
//		memset(path, '\0', strlen(seis_file)+1);
//		strncpy(path, seis_file, last_slash-seis_file);
//		makedir(path);
//	}
//
//	if (merge_output==0) {
//		if (sgtfilepar.xfile[0]!='\0') {
//			header->comps |= X_COMP_FLAG;
//			write_seis(header, seis_file,sname,"000",sn,&dtout,ntout,&tstart,output_binary);
//			(*seis_return)[i] = sn;
//		}
//		if (sgtfilepar.yfile[0]!='\0') {
//			header->comps |= Y_COMP_FLAG;
//			write_seis(header, seis_file,sname,"090",se,&dtout,ntout,&tstart,output_binary);
//			(*seis_return)[i] = se;
//		}
//		if (sgtfilepar.zfile[0]!='\0') {
//			header->comps |= Z_COMP_FLAG;
//			write_seis(header, seis_file,sname,"ver",sv,&dtout,ntout,&tstart,output_binary);
//			(*seis_return)[i] = sv;
//		}
//	} else { //merging output
/*		if (sgtfilepar.xfile[0]!='\0') {
			if (sgtfilepar.yfile[0]!='\0') {
				if (sgtfilepar.zfile[0]!='\0') { //x,y,z
					header->comps = (X_COMP_FLAG | Y_COMP_FLAG | Z_COMP_FLAG);
					memcpy(writing_buffer+)
//					send_data_file(header, seis_file, header->source_id, header->rupture_id, sv, sizeof(float)*3*ntout, my_id);
//					write_seis(header, seis_file,sname,"grm",sv,&dtout,3*ntout,&tstart,output_binary);
				} else { //x,y
					header->comps = X_COMP_FLAG | Y_COMP_FLAG;
					send_data_file(header, seis_file, header->source_id, header->rupture_id, sn, sizeof(float)*2*ntout, my_id);
//					write_seis(header, seis_file,sname,"grm",sn,&dtout,2*ntout,&tstart,output_binary);
					(*seis_return)[i] = sn;
				}
			}
			else if (sgtfilepar.zfile[0]!='\0') { //x,z
		                header->comps = X_COMP_FLAG | Z_COMP_FLAG;
						send_data_file(header, seis_file, header->source_id, header->rupture_id, sv, sizeof(float)*2*ntout, my_id);
//        		        write_seis(header, seis_file,sname,"grm",sv,&dtout,2*ntout,&tstart,output_binary);
			} else { //x
                		header->comps = X_COMP_FLAG;
    					send_data_file(header, seis_file, header->source_id, header->rupture_id, sn, sizeof(float)*ntout, my_id);
//                		write_seis(header, seis_file,sname,"grm",sn,&dtout,ntout,&tstart,output_binary);
                		(*seis_return)[i] = sn;
			}
		} else if (sgtfilepar.yfile[0]!='\0') {
			if (sgtfilepar.zfile[0]!='\0') { //y,z
				printf("Can't output merged Y and Z components.\n");
				exit(2);
			} else { //y
				header->comps = Y_COMP_FLAG;
				send_data_file(header, seis_file, header->source_id, header->rupture_id, se, sizeof(float)*ntout, my_id);
//				write_seis(header, seis_file,sname,"grm",se,&dtout,ntout,&tstart,output_binary);
				(*seis_return)[i] = se;
			}
		} else if (sgtfilepar.zfile[0]!='\0') { //z
			header->comps = Z_COMP_FLAG;
			send_data_file(header, seis_file, header->source_id, header->rupture_id, sv, sizeof(float)*ntout, my_id);
//			write_seis(header, seis_file,sname,"grm",sv,&dtout,ntout,&tstart,output_binary);
		}
	}*/

//Add 1 to ending rup_var_id because it's exclusive
send_data_cluster(seis_file, header->source_id, header->rupture_id, rup_vars[0].rup_var_id, rup_vars[num_rup_vars-1].rup_var_id+1, writing_buffer, (sizeof(struct seisheader)+entry_size)*num_rup_vars, my_id);

free(writing_buffer);

fprintf(stdout,"Total moment= %13.5e\n",tmom);

for (i=0; i<num_rup_vars; i++) {
	free_srf_ptrs(&(srf[i]));
}
free(srf);
free(prect_ptr);
free(prseg_ptr);
free(apnts_ptr);
free(tmom);
free(apv_off);

return seis;
}
