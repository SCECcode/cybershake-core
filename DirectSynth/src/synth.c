/*
 * synth.c
 *
 *  Created on: Jan 16, 2015
 *      Author: scott
 */

#include "include.h"
#include "defs.h"
#include "structure.h"
#include "duration.h"
#include "functions.h"

/*
 * 		Get list of SGTs needed
		Construct requests of (process, SGTs) (if request is too big, split it into pieces)
		for each request:
			Request data from process (<N)
			Perform convolution on data and SRF
		Calculate PSA, RotD
		Send output data to process N to aggregate and write
 *
 */
int get_handler_for_sgt_index(long long indx, long long* sgt_cutoffs, int num_sgt_handlers, int my_id);
void send_data_file(struct seisheader* header, char data_filename[256], int src_id, int rup_id, void* buf, int data_size_bytes, int my_id);

//This duplicates some of seis_psa.c
int run_synth(task_info* t_info, int* proc_points, int num_sgt_handlers, char stat[64], float slat, float slon, int run_id, float det_max_freq, float stoch_max_freq, int run_PSA, int run_rotd, int run_duration, int my_id) {
    struct geoprojection geop;
    struct sgtindex statindx;
    int ip, i;
    int print_tol = 25;
    float elon, elat, edep;

    //Do some initializing
    set_geoproj(t_info->sgtmast,&geop);
    struct sgtindex eqindx;
    eqindx.h = t_info->sgtindx[0].h;
    statindx.h = t_info->sgtindx[0].h;
    float sdep = 0.0;
    get_indx(&slon,&slat,&sdep,&statindx,&geop);

    struct rup_geom_point* rg_points;
    int num_srf_pts = parse_rup_geom(t_info->task->rupture_filename, &rg_points);

    //This is freed inside jbsim3d
    struct sgtparams* sgtparms = (struct sgtparams *) check_malloc ((num_srf_pts)*sizeof(struct sgtparams));

    int ptol = print_tol;
    float maxdelta = 0.0;
    float fweight = 1.0;
    int non_exact = 0;
	int nm, memlen;

	if (debug) write_log("Converting SGT lat/lon/dep coordinates to X/Y/Z.");
	long long* indx_master;
    float find_sgt_time = 0.0;
    struct timeval start_tv, end_tv;
    for(ip=0;ip<num_srf_pts;ip++) {
    	elon = rg_points[ip].lon;
        elat = rg_points[ip].lat;
        edep = rg_points[ip].dep;
        get_indx(&elon,&elat,&edep,&eqindx,&geop);
	
	gettimeofday(&start_tv, NULL);
        find_sgt(&sgtparms[ip],t_info->sgtmast,t_info->sgtindx,&eqindx,&statindx,&maxdelta,&fweight);
	gettimeofday(&end_tv, NULL);
	find_sgt_time += (end_tv.tv_sec - start_tv.tv_sec) + (end_tv.tv_usec - start_tv.tv_usec)/1000000.0;

        if (sgtparms[ip].nsgt != 1)
        	non_exact++;

           if((float)(100.0*(float)(ip+1)/(float)(num_srf_pts)) >= ptol) {
              fprintf(stdout," %3d percent done (%d of %d)\n",ptol,ip,num_srf_pts);
              ptol = ptol + print_tol;
           }
    }
    if (debug) {
	char buf[256];
	sprintf(buf, "Spent %f seconds in find_sgt function.", find_sgt_time);
	write_log(buf);
    }

	indx_master = (long long *) check_malloc (4*(num_srf_pts)*sizeof(long long));
	get_master_list_opt(sgtparms,num_srf_pts,indx_master,&nm);
	indx_master = (long long *) check_realloc (indx_master,nm*sizeof(long long));

	fprintf(stdout,"nm= %d non_exact= %d\n",nm,non_exact);

	//synthesis
	struct seisheader header;
	strcpy(header.version, "12.10");
	strcpy(header.site_name, stat);
	memset(header.padding, 0, 8);
	header.source_id = t_info->task->source_id;
	header.rupture_id = t_info->task->rupture_id;
	header.det_max_freq = det_max_freq;
	header.stoch_max_freq = stoch_max_freq;
	//get dt, nt, num_comps from jbsim3d function

	if (debug) write_log("Creating rupture variations from task list.");
	//Create array of rupture_variation structures from task list
	int num_rup_vars = t_info->task->ending_rv_id - t_info->task->starting_rv_id;
	rup_var_struct* rup_vars = check_malloc(sizeof(rup_var_struct) * num_rup_vars);
	for (i=0; i<num_rup_vars; i++) {
		rup_vars[i].rup_var_id = t_info->task->starting_rv_id + i;
		rup_vars[i].slip_id = rup_vars[i].rup_var_id/t_info->task->num_hypos;
		rup_vars[i].hypo_id = rup_vars[i].rup_var_id % t_info->task->num_hypos;
	}

	//Create array of SGT requests
	//one entry in array for each SGT handler
	if (debug) write_log("Creating array of SGT requests.");
	int step_size = nm/num_sgt_handlers;
	if (step_size<10) {
		step_size = 10;
	}
	long long** sgts_by_handler = check_malloc(sizeof(long long*)*num_sgt_handlers);
	int* num_sgts_by_handler = check_malloc(sizeof(int)*num_sgt_handlers);
	for (i=0; i<num_sgt_handlers; i++) {
		sgts_by_handler[i] = check_malloc(sizeof(long long)*step_size);
		num_sgts_by_handler[i] = 0;
	}
	//Convert proc_point values, which are the index into sgtindx, into the long long SGT index values, so we can search on them and determine the handler
	long long* sgt_cutoffs = check_malloc(sizeof(long long)*(num_sgt_handlers+1));
	//Ignore processor 0, since it doesn't have any points
	sgt_cutoffs[0] = 0;
	for (i=1; i<num_sgt_handlers; i++) {
		sgt_cutoffs[i] = (t_info->sgtindx)[proc_points[i]].indx;
	}
	//for the sgt_cutoffs, make it one unit past the last index
	sgt_cutoffs[num_sgt_handlers] = (t_info->sgtindx)[proc_points[num_sgt_handlers]-1].indx + 1;
	for (i=0; i<nm; i++) {
		int handler = get_handler_for_sgt_index(indx_master[i], sgt_cutoffs, num_sgt_handlers, my_id);
		sgts_by_handler[handler][num_sgts_by_handler[handler]] = indx_master[i];
		num_sgts_by_handler[handler]++;
		if (num_sgts_by_handler[handler] % step_size == 0) {
			sgts_by_handler[handler] = check_realloc(sgts_by_handler[handler], sizeof(long long)*(num_sgts_by_handler[handler]+step_size));
		}
	}
	free(sgt_cutoffs);
	//Now shrink to fit
	for (i=0; i<num_sgt_handlers; i++) {
		int shrink_to = sizeof(long long)*num_sgts_by_handler[i];
		//realloc gets upset if we use size 0
		if (shrink_to<1) {
			shrink_to = 1;
		}
		sgts_by_handler[i] = check_realloc(sgts_by_handler[i], shrink_to);
	}

	if (debug) {
		char buf[256];
		for (i=0; i<num_sgt_handlers; i++) {
			sprintf(buf, "Will request %d SGTs from handler %d.", num_sgts_by_handler[i], i);
			write_log(buf);
		}
	}	

	//Create array of seis ptrs
	float** seis = check_malloc(sizeof(float*) * num_rup_vars);

	//Arguments to jbsim3d_synth
	int ntout = -1;
	mstpar("ntout","d",&ntout);
	char seis_filename[256];
	sprintf(seis_filename, "Seismogram_%s_%d_%d_%d.grm", stat, run_id, t_info->task->source_id, t_info->task->rupture_id);

	if (debug) write_log("Starting jbsim3d_synth.");
	float** original_seis = jbsim3d_synth(&seis, &header, stat, slon, slat, ntout, seis_filename,
			t_info->task->rupture_filename, t_info->sgtfilepar, sgtparms, *(t_info->sgtmast), t_info->sgtindx, geop, indx_master, nm, sgts_by_handler,
			num_sgts_by_handler, num_sgt_handlers, num_rup_vars, rup_vars, my_id);

	for(i=0; i<num_sgt_handlers; i++) {
		free(sgts_by_handler[i]);
	}
	free(num_sgts_by_handler);
	free(sgts_by_handler);

	if (run_PSA==0) {
		return 0;
	} else {
		printf("Calculating PSA.\n");
	}

	//psa
	int nx, ny, npts;
	float dt;
	char seis_units[120];
	char output_units[120];
	char output_type[120];
	char period[120];
	float filter_high_hz;
	char byteswap[120];
	char input_file[256];
	char output_file[256];
	int seis_units_len, output_units_len, output_type_len, period_len, byteswap_len, input_file_len, output_file_len;

	memset(seis_units, ' ', 120);
	memset(output_units, ' ', 120);
	memset(output_type, ' ', 120);
	memset(period, ' ', 120);
	memset(byteswap, ' ', 120);
	memset(input_file, ' ', 256);
	memset(output_file, ' ', 256);

	//Put this back in b/c we called setpar and endpar in the rupture variation generator
	setpar(t_info->argc, t_info->argv);
	mstpar("simulation_out_pointsX","d",&nx);
	mstpar("simulation_out_pointsY","d",&ny);
	mstpar("simulation_out_timesamples","d",&npts);
	mstpar("simulation_out_timeskip","f",&dt);
	mstpar("surfseis_rspectra_seismogram_units","s",seis_units);
	mstpar("surfseis_rspectra_output_units","s",output_units);
	mstpar("surfseis_rspectra_output_type","s",output_type);
	mstpar("surfseis_rspectra_period","s",period);
	mstpar("surfseis_rspectra_apply_filter_highHZ","f",&filter_high_hz);
	mstpar("surfseis_rspectra_apply_byteswap","s",byteswap);

	//Make sure dt and dtout (from the header) are the same
	if (fabs(dt-header.dt)>0.00001) {
		fprintf(stderr, "%d) dt from 'simulation_out_timeskip' argument is %f, but dt from 'dtout' argument and used in seismogram header is %f.  These must match, aborting.\n", my_id, dt, header.dt);
		MPI_Abort(MPI_COMM_WORLD, 7);
		exit(7);
	}

//	mstpar("out","s",&output_file);
	///Construct output filename
	output_file_len = sprintf(output_file, "PeakVals_%s_%d_%d_%d.bsa", stat, run_id, t_info->task->source_id, t_info->task->rupture_id);
	printf("%d) output_file = %s, len = %d\n", my_id, output_file, output_file_len);
	char output_file_ext[256];
	sprintf(output_file_ext, output_file);
	printf("%d) output_file_ext = %s\n", my_id, output_file_ext);
	
	//replace null terminator with space
	//fortran expects all the strings to be space-terminated
	seis_units_len = strlen(seis_units);
	seis_units[seis_units_len] = ' ';
	output_units_len = strlen(output_units);
	output_units[output_units_len] = ' ';
	output_type_len = strlen(output_type);
	output_type[output_type_len] = ' ';
	period_len = strlen(period);
	period[period_len] = ' ';
	byteswap_len = strlen(byteswap);
	byteswap[byteswap_len] = ' ';
	input_file_len = strlen(input_file);
	input_file[input_file_len] = ' ';
	//output_file_len = strlen(output_file);
	output_file[output_file_len] = ' ';
	printf("%d) output_file_len = %d\n", my_id, output_file_len);

	//This tells spectrad to just return the data, not to do the output itself
	int output_option = 2;
	//Supply an array to get the data1 back from spectrad
	float* psa_data = check_malloc(sizeof(float)*MAXPERIODS);
	//Accumulation array for writing to file
	int psa_rv_file_size = sizeof(struct seisheader) + sizeof(float)*nx*ny*NUM_SCEC_PERIODS;
	char* psa_file_buffer = check_malloc(psa_rv_file_size*num_rup_vars);
	int rotd_rv_file_size = sizeof(struct seisheader) + sizeof(int) + sizeof(struct rotD_record)*NUM_ROTD_PERIODS;
	char* rotd_file_buffer = check_malloc(rotd_rv_file_size*num_rup_vars);
	//Use nx*ny here to capture the number of components
	int duration_rv_file_size = sizeof(struct seisheader) + sizeof(int) + sizeof(struct duration_record)*nx*ny*NUM_DURATION_MEASURES;
	char* duration_file_buffer = check_malloc(duration_rv_file_size*num_rup_vars);

	if (MAXPERIODS<NUM_SCEC_PERIODS*nx*ny) {
		fprintf(stderr, "%d) Error!  Need to increase MAXPERIODS in defs.h and maxperiods in surfseis_rspectra.f.  Aborting.", my_id);
		MPI_Finalize();
		exit(7);
	}

	int rv;
	struct rotD_record* rotD_data = NULL;
	struct duration_record* duration_data = NULL;
	char rotd_filename[256];
	char duration_filename[256];
	for (rv=0; rv<num_rup_vars; rv++) {
		if (debug) {
			char buf[256];
			sprintf(buf, "Performing PSA for rv %d.", rv);
			write_log(buf);
		}
		header.rup_var_id = rup_vars[rv].rup_var_id;
		//printf("Entering spectrad.\n");
		       //printf("locations: %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld\n", &nx, &ny, &npts, &dt, seis_units, output_units, output_type, period, &filter_high_hz, byteswap, input_file, output_file, seis);
//		fflush(stdout);
		//spectrad_(&header, &nx, &ny, &npts, &dt, seis_units, seis_units_len, output_units, output_units_len, output_type, output_type_len, period, period_len, &filter_high_hz, byteswap, byteswap_len, input_file, input_file_len, output_file, output_file_len, seis, &using_pipe_forwarding);
		spectrad_(&header, &nx, &ny, &npts, &dt, seis_units, output_units, output_type, period, &filter_high_hz, byteswap, input_file, output_file, seis[rv], &output_option, psa_data,
					seis_units_len, output_units_len, output_type_len, period_len, byteswap_len, input_file_len, output_file_len);
		//Replace null terminator so send_data_file doesn't have an issue
		printf("%d) output_file: %s\n", my_id, output_file_ext);
		memcpy(psa_file_buffer+rv*psa_rv_file_size, &header, sizeof(struct seisheader));
		memcpy(psa_file_buffer+rv*psa_rv_file_size+sizeof(struct seisheader), psa_data, nx*ny*NUM_SCEC_PERIODS*sizeof(float));
		//send_data_file(&header, output_file_ext, t_info->task->source_id, t_info->task->rupture_id, psa_data, nx*ny*NUM_SCEC_PERIODS*sizeof(float), my_id);

		if (run_rotd) {
			if (debug) {
				char buf[256];
				sprintf(buf, "Performing RotD for rv %d.", rv);
				write_log(buf);
			}
			if (rotD_data==NULL) {
				rotD_data = check_malloc(sizeof(int) + sizeof(struct rotD_record)*NUM_ROTD_PERIODS);
			}
			//We offset by sizeof(int) because we'll put the # of periods into the first spot
			char* start_ptr = ((char*)(rotD_data))+sizeof(int);
			//printf("%d) rotD_data @ %ld, start_ptr @ %ld\n", my_id, rotD_data, start_ptr);
			//if (rotD_data==NULL) {
			//	rotD_data = check_malloc(sizeof(struct rotD_record)*NUM_ROTD_PERIODS);
			//}
			int rc = rotd(&header, seis[rv], (struct rotD_record*)start_ptr);
			//int rc = rotd(&header, seis[rv], rotD_data);
			if (rc!=0) {
				fprintf(stderr, "%d) Error in rotd code, aborting.\n", my_id);
				if (debug) close_log();
				MPI_Finalize();
				exit(rc);
			}
			//Also need to send # of periods as part of the file
			//Need to put it in an int for memcpy
			int num_rotd_periods = NUM_ROTD_PERIODS;
			//char* file_buffer = check_malloc(sizeof(int) + sizeof(struct rotD_record)*NUM_ROTD_PERIODS);
			//memcpy(file_buffer, &num_rotd_periods, sizeof(int));
			//memcpy(file_buffer+sizeof(int), rotD_data, sizeof(struct rotD_record)*NUM_ROTD_PERIODS);
			memcpy(rotD_data, &num_rotd_periods, sizeof(int));
			memcpy(rotd_file_buffer+rv*rotd_rv_file_size, &header, sizeof(struct seisheader));
			memcpy(rotd_file_buffer+rv*rotd_rv_file_size+sizeof(struct seisheader), rotD_data, sizeof(int)+sizeof(struct rotD_record)*NUM_ROTD_PERIODS);
			//send_data_file(&header, rotd_filename, t_info->task->source_id, t_info->task->rupture_id, rotD_data, sizeof(int)+NUM_ROTD_PERIODS*sizeof(struct rotD_record), my_id);
			//send_data_file(&header, rotd_filename, file_buffer, sizeof(int)+NUM_ROTD_PERIODS*sizeof(struct rotD_record), my_id);
			//free(file_buffer);
		}

		if (run_duration) {
			if (debug) {
				char buf[256];
				sprintf(buf, "Performing duration calculations for rv %d.", rv);
				write_log(buf);
			}
			if (duration_data==NULL) {
				duration_data = check_malloc(sizeof(int)+sizeof(struct duration_record)*nx*ny*NUM_DURATION_MEASURES);
			}
			//Put the # of measures into the first spot
			int num_dur_measures = NUM_DURATION_MEASURES;
			memcpy(duration_data, &num_dur_measures, sizeof(int));
			//Offset by int
			char* start_ptr = ((char*)(duration_data))+sizeof(int);
			int rc = duration(&header, seis[rv], (struct duration_record*)start_ptr);
			if (rc!=0) {
				fprintf(stderr, "%d) Error in duration code, aborting.\n", my_id);
				if (debug) close_log();
				MPI_Finalize();
				exit(rc);
			}
			struct duration_record* ad_rec = (struct duration_record*)start_ptr + 4;
			//printf("AD value after function: %d, %d, %d, %f\n", ad_rec->type, ad_rec->type_value, ad_rec->component, ad_rec->value);
			//Add header info to the file buffer
			memcpy(duration_file_buffer+rv*duration_rv_file_size, &header, sizeof(struct seisheader));
			memcpy(duration_file_buffer+rv*duration_rv_file_size+sizeof(struct seisheader), duration_data, sizeof(int)+sizeof(struct duration_record)*nx*ny*NUM_DURATION_MEASURES);
			//printf("Arias duration value in duration_file_buffer: %f\n", duration_file_buffer+rv*duration_rv_file_size+sizeof(struct seisheader)+sizeof(int)+76);
		}			
	}

	if (run_PSA) {
		send_data_cluster(output_file_ext, header.source_id, header.rupture_id, rup_vars[0].rup_var_id, rup_vars[num_rup_vars-1].rup_var_id+1, psa_file_buffer, num_rup_vars*psa_rv_file_size, my_id);
	}
	if (run_rotd) {
		sprintf(rotd_filename, "RotD_%s_%d_%d_%d.rotd", stat, run_id, t_info->task->source_id, t_info->task->rupture_id);
		send_data_cluster(rotd_filename, header.source_id, header.rupture_id, rup_vars[0].rup_var_id, rup_vars[num_rup_vars-1].rup_var_id+1, rotd_file_buffer, num_rup_vars*rotd_rv_file_size, my_id);
	}
	if (run_duration) {
		sprintf(duration_filename, "Duration_%s_%d_%d_%d.dur", stat, run_id, t_info->task->source_id, t_info->task->rupture_id);
		send_data_cluster(duration_filename, header.source_id, header.rupture_id, rup_vars[0].rup_var_id, rup_vars[num_rup_vars-1].rup_var_id+1, duration_file_buffer, num_rup_vars*duration_rv_file_size, my_id);
	}

	free(psa_file_buffer);
	free(rotd_file_buffer);
	free(duration_file_buffer);

	free(psa_data);
	free(rotD_data);
	free(duration_data);
	free(sgtparms);
	free(indx_master);
	free(rup_vars);
	for (rv=0; rv<num_rup_vars; rv++) {
		free(original_seis[rv]);
	}
	free(original_seis);

	free(rg_points);

	return 0;
}

//Determine which handler is responsible for this SGT point
int get_handler_for_sgt_index(long long indx, long long* sgt_cutoffs, int num_sgt_handlers, int my_id) {
	//SGTs were sorted by lexicographical order
	//Do binary search for index just smaller than
	int start_index = 1;
	int end_index = num_sgt_handlers;
	int mid_index = (start_index + end_index) / 2;
	while(end_index>start_index+1) {
		mid_index = (start_index + end_index) / 2;
		if (indx>sgt_cutoffs[mid_index]) {
			start_index = mid_index;
		} else if (indx<sgt_cutoffs[mid_index]) {
			end_index = mid_index;
		} else {
			start_index = mid_index;
			end_index = start_index+1;
			break;
		}
	}
	//Verify that this is right
	if (!(indx>=sgt_cutoffs[start_index]) || !(indx<sgt_cutoffs[end_index])) {
		fprintf(stderr, "%d) Error in determining SGT handler for index %ld: we thought it was on handler %d between indices %ld and %ld, but it's not, so aborting.\n", my_id, indx, start_index, sgt_cutoffs[start_index], sgt_cutoffs[end_index]);
		if (debug) close_log();
		MPI_Finalize();
		exit(5);
	}
	return start_index;
}


void request_sgt(struct sgtheader* sgthead, float* sgtbuf, int num_comps, int request_from_handler_id, long long** sgts_by_handler, int starting_index, int ending_index, int nt, int my_id) {
	int i, j, k;
	MPI_Datatype sgtheader_type, sgtdata_type;
	construct_sgtheader_datatype(&sgtheader_type);
	construct_sgtdata_datatype(&sgtdata_type, nt);

	handler_msg msg;
	msg.msg_src = my_id;
	msg.msg_type = SGT_REQUEST;
	msg.num_sgts_requested = ending_index - starting_index;
	if (debug) {
		char buf[256];
		sprintf(buf, "Requesting %d SGTs from handler %d.", msg.num_sgts_requested, request_from_handler_id);
		write_log(buf);
	}
	check_send(&msg, 3, MPI_INT, request_from_handler_id, SGT_REQUEST_TAG, MPI_COMM_WORLD, "Error sending SGT request, aborting.", my_id);
	if (debug) write_log("Sending list of SGTs.");
	check_send(sgts_by_handler[request_from_handler_id]+starting_index, msg.num_sgts_requested, MPI_LONG_LONG, request_from_handler_id, SGT_REQUEST_TAG, MPI_COMM_WORLD, "Error sending SGT request, aborting.", my_id);
	if (debug) write_log("Receiving SGT header data.");
	check_recv(sgthead, msg.num_sgts_requested, sgtheader_type, request_from_handler_id, SGT_HEADER_DATA_TAG, MPI_COMM_WORLD, "Error receiving SGT header information, aborting.", my_id);
	float* buffer = check_malloc(msg.num_sgts_requested*sizeof(float)*6*nt);
	for(i=0; i<num_comps; i++) {
		if (debug) {
			char buf[256];
			sprintf(buf, "Receiving raw SGT data for component %d.", i);
			write_log(buf);
		}
		check_recv(buffer, msg.num_sgts_requested, sgtdata_type, request_from_handler_id, SGT_RAW_DATA_TAG, MPI_COMM_WORLD, "Error receiving raw SGT data, aborting.", my_id);
		//copy it into sgtbuf in the correct location, with component offsets
		for(j=0; j<msg.num_sgts_requested; j++) {
			memcpy(sgtbuf+(long long)(j*18+6*i)*nt, buffer+(long long)j*6*nt, 6*nt*sizeof(float));
		}
	}
	float xmom = sgthead[0].xmom;
	//We are assuming here that the xmoment and ymoment are the same; if this is a bad assumption, we've already aborted in the SGT handlers
	float ymom = xmom;
	float zmom = ymom;
	if (num_comps>=2) {
		for(j=0; j<msg.num_sgts_requested; j++) {
			sgthead[j].ymom = ymom;
		}
	}
	if (num_comps==3) {
		for(j=0; j<msg.num_sgts_requested; j++) {
                        sgthead[j].zmom = zmom;
                }
	}

	free(buffer);
	MPI_Type_free(&sgtheader_type);
        MPI_Type_free(&sgtdata_type);
}

void send_data_cluster(char data_filename[256], int src_id, int rup_id, int start_rv, int end_rv, void* data, int data_size_bytes, int my_id) {
	master_msg msg;
	msg.msg_src = my_id;
	msg.msg_type = DATA_FILE;
	msg.msg_data = data_size_bytes;
	if (debug) {
		char buf[256];
		sprintf(buf, "Notifying master of request to send %d bytes for data file %s.", msg.msg_data, data_filename);
		write_log(buf);
	}
	check_send(&msg, 3, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD, "Error sending data incoming message to master, aborting.", my_id);
	data_file_metadata df;
	df.src_id = src_id;
	df.rup_id = rup_id;
	df.starting_rv = start_rv;
	df.ending_rv = end_rv;
	strcpy(df.filename, data_filename);
	if (debug) write_log("Sending data file metadata to master.");
	check_send(&df, 4*sizeof(int)+256, MPI_BYTE, 0, DATA_FILENAME_TAG, MPI_COMM_WORLD, "Error sending data file metadata to master, aborting.", my_id);
	if (debug) write_log("Sending data contents to master.");
	check_send(data, msg.msg_data, MPI_BYTE, 0, DATA_TAG, MPI_COMM_WORLD, "Error sending data contents to master, aborting.", my_id);
}

/*void send_data_file(struct seisheader* header, char data_filename[256], int src_id, int rup_id, void* buf, int data_size_bytes, int my_id) {
	master_msg msg;
	msg.msg_src = my_id;
	msg.msg_type = DATA_FILE;
	msg.msg_data = sizeof(struct seisheader) + data_size_bytes;
	if (debug) {
		char buf[256];
		sprintf(buf, "Notifying master of request to send %d bytes for data file %s.", msg.msg_data, data_filename);
		write_log(buf);
	}
	check_send(&msg, 3, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD, "Error sending data incoming message to master, aborting.", my_id);
	data_file_metadata df;
	df.src_id = src_id;
	df.rup_id = rup_id;
	strcpy(df.filename, data_filename);
	if (debug) write_log("Sending data file metadata to master.");
	check_send(&df, sizeof(int)+sizeof(int)+256, MPI_BYTE, 0, DATA_FILENAME_TAG, MPI_COMM_WORLD, "Error sending data file metadata to master, aborting.", my_id);
	char* data = check_malloc(sizeof(char)*msg.msg_data);
	memcpy(data, header, sizeof(struct seisheader));
	memcpy(data + sizeof(struct seisheader), buf, data_size_bytes);
	if (debug) write_log("Sending data contents to master.");
	check_send(data, msg.msg_data, MPI_BYTE, 0, DATA_TAG, MPI_COMM_WORLD, "Error sending data contents to master, aborting.", my_id);
	free(data);
}*/
