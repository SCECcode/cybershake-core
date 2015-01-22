/*
 * log.c
 *
 *  Created on: Jan 22, 2015
 *      Author: scott
 */

#include "include.h"

FILE* log_fp = NULL;

void open_log(int my_id) {
	if (log_fp==NULL) {
		char log_filename[256];
		sprintf(log_filename, "log.%d", my_id);
		log_fp = fopen(log_filename, "w");
	}
}

void write_log(char* string) {
	struct timeval tv;
	time_t cur_time;
	struct tm* time_struct;
	char time_buf[60];
	gettimeofday(&tv, NULL);
	cur_time = tv.tv_sec;
	time_struct = localtime(&cur_time);
	strftime(time_buf, 59, "%b %d %T", time_struct);
	fprintf(log_fp, "%s.%06d> %s\n", time_buf, tv.tv_usec, string);
}

void close_log() {
	fflush(log_fp);
	fclose(log_fp);
}
