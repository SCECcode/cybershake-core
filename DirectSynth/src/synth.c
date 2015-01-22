/*
 * synth.c
 *
 *  Created on: Jan 16, 2015
 *      Author: scott
 */

#include "include.h"
#include "structure.h"
#include "functions.h"
#include "defs.h"

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
void send_data_file(struct seisheader* header, char data_filename[256], void* buf, int data_size_bytes, int my_id);

//This duplicates some of seis_psa.c
int run_synth(task_info* t_info, int* proc_points, int num_sgt_handlers, char stat[64], float slat, float slon, int run_id, float det_max_freq, float stoch_max_freq, int run_PSA, int run_rotd, int my_id) {
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
    for(ip=0;ip<num_srf_pts;ip++) {
    	elon = rg_points[ip].lon;
        elat = rg_points[ip].lat;
        edep = rg_points[ip].dep;
        get_indx(&elon,&elat,&edep,&eqindx,&geop);

        find_sgt(&sgtparms[ip],t_info->sgtmast,t_info->sgtindx,&eqindx,&statindx,&maxdelta,&fweight);

        if (sgtparms[ip].nsgt != 1)
        	non_exact++;

           if((float)(100.0*(float)(ip+1)/(float)(num_srf_pts)) >= ptol) {
              fprintf(stdout," %3d percent done (%d of %d)\n",ptol,ip,num_srf_pts);
              ptol = ptol + print_tol;
           }
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
	long long** sgts_by_handler = check_malloc(sizeof(long long*)*num_sgt_handlers);
	int* num_sgts_by_handler = check_malloc(sizeof(int)*num_sgt_handlers);
	for (i=0; i<num_sgt_handlers; i++) {
		sgts_by_handler[i] = check_malloc(sizeof(long long)*step_size);
		num_sgts_by_handler[i] = 0;
	}
	//Convert proc_point values, which are the index into sgtindx, into the long long SGT index values, so we can search on them and determine the handler
	long long* sgt_cutoffs = check_malloc(sizeof(long long)*num_sgt_handlers);
	for (i=0; i<num_sgt_handlers; i++) {
		sgt_cutoffs[i] = (t_info->sgtindx)[proc_points[i]].indx;
	}

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
		sgts_by_handler[i] = check_realloc(sgts_by_handler[i], sizeof(long long)*num_sgts_by_handler[i]);
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
	//setpar(argc, argv);
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

//	mstpar("out","s",&output_file);
	///Construct output filename
	sprintf(output_file, "PeakVals_%s_%d_%d_%d.bsa", stat, run_id, t_info->task->source_id, t_info->task->rupture_id);

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
	output_file_len = strlen(output_file);
	output_file[output_file_len] = ' ';

	//This tells spectrad to just return the data, not to do the output itself
	int output_option = 2;
	//Supply an array to get the data back from spectrad
	float* psa_data = check_malloc(sizeof(float)*MAXPERIODS);

	if (MAXPERIODS<NUM_SCEC_PERIODS*nx*ny) {
		fprintf(stderr, "%d) Error!  Need to increase MAXPERIODS in defs.h and maxperiods in surfseis_rspectra.f.  Aborting.", my_id);
		MPI_Finalize();
		exit(7);
	}

	int rv;
	struct rotD_record* rotD_records = NULL;
	char rotd_filename[256];
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
		send_data_file(&header, output_file, psa_data, header.comps*NUM_SCEC_PERIODS*sizeof(float), my_id);

		if (run_rotd) {
			if (debug) {
				char buf[256];
				sprintf(buf, "Performing RotD for rv %d.", rv);
				write_log(buf);
			}
			if (rotD_records==NULL) {
				rotD_records = check_malloc(sizeof(struct rotD_record)*NUM_ROTD_PERIODS);
			}
			int rc = rotd(&header, seis[rv], rotD_records);
			if (rc!=0) {
				fprintf(stderr, "%d) Error in rotd code, aborting.\n", my_id);
				if (debug) close_log();
				MPI_Finalize();
				exit(rc);
			}
			sprintf(rotd_filename, "RotD_%s_%d_%d_%d.rotd", stat, run_id, t_info->task->source_id, t_info->task->rupture_id);
			send_data_file(&header, rotd_filename, rotD_records, NUM_ROTD_PERIODS*sizeof(struct rotD_record), my_id);
		}
	}

	free(psa_data);
	free(rotD_records);

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
			start_index = mid_index+1;
		} else if (indx<sgt_cutoffs[mid_index]) {
			end_index = mid_index-1;
		} else {
			start_index = mid_index;
			end_index = start_index+1;
			break;
		}
	}
	//Verify that this is right
	if (!(indx>sgt_cutoffs[start_index]) || !(indx<sgt_cutoffs[end_index])) {
		fprintf(stderr, "%d) Error in determining SGT handler for index %ld: we thought it was on handler %d between indices %ld and %ld, but it's not, so aborting.\n", my_id, indx, start_index, sgt_cutoffs[start_index], sgt_cutoffs[end_index]);
		if (debug) close_log();
		MPI_Finalize();
		exit(5);
	}
	return start_index;
}


void request_sgt(struct sgtheader* sgthead, float* sgtbuf, int request_from_handler_id, long long** sgts_by_handler, int starting_index, int ending_index, int nt, int my_id) {
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
	check_send(sgts_by_handler[request_from_handler_id]+starting_index, ending_index, MPI_LONG_LONG, request_from_handler_id, SGT_REQUEST_TAG, MPI_COMM_WORLD, "Error sending SGT request, aborting.", my_id);
	if (debug) write_log("Receiving SGT header data.");
	check_recv(sgthead, msg.num_sgts_requested, sgtheader_type, request_from_handler_id, SGT_HEADER_DATA_TAG, MPI_COMM_WORLD, "Error receiving SGT header information, aborting.", my_id);
	if (debug) write_log("Receiving raw SGT data.");
	check_recv(sgtbuf, msg.num_sgts_requested, sgtdata_type, request_from_handler_id, SGT_RAW_DATA_TAG, MPI_COMM_WORLD, "Error receiving raw SGT data, aborting.", my_id);
}


void send_data_file(struct seisheader* header, char data_filename[256], void* buf, int data_size_bytes, int my_id) {
	master_msg msg;
	msg.msg_src = my_id;
	msg.msg_type = DATA_FILE;
	msg.msg_data = sizeof(struct seisheader) + data_size_bytes;
	if (debug) {
		char buf[256];
		sprintf(buf, "Notifying master of request to send %d bytes for data file %s.", data_size_bytes, data_filename);
		write_log(buf);
	}
	check_send(&msg, 3, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD, "Error sending data incoming message to master, aborting.", my_id);
	data_file df;
	strcpy(df.filename, data_filename);
	df.data = check_malloc(sizeof(char)*msg.msg_data);
	memcpy(df.data, header, sizeof(struct seisheader));
	memcpy(df.data + sizeof(struct seisheader), buf, data_size_bytes);
	if (debug) write_log("Sending data to master.");
	check_send(&msg, 256+msg.msg_data, MPI_BYTE, 0, DATA_TAG, MPI_COMM_WORLD, "Error sending data contents to master, aborting.", my_id);
	free(df.data);
}
