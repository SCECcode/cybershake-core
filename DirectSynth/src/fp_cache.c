/*
 * fp_cache.c
 *
 *  Created on: Jan 14, 2015
 *      Author: scott
 */

#include "include.h"
#include "structure.h"
#include "defs.h"
#include "functions.h"

const int CACHE_SIZE = 100;
fp_cache_entry* fp_cache = NULL;

void create_fp_cache() {
	if (debug) write_log("Creating fp cache.");
	fp_cache = check_malloc(sizeof(fp_cache_entry)*CACHE_SIZE);
	int i;
	for (i=0; i<CACHE_SIZE; i++) {
		fp_cache[i].filename[0] = '\0';
		fp_cache[i].fp = NULL;
	}
}

FILE* find_and_use_fp(char* filename) {
	if (fp_cache==NULL) {
		create_fp_cache();
	}
	//Locate filename in the cache
	//If it's not present, create a new entry
	//Move filename entry to index 0
	//Adjust everyone else accordingly
	int cur_index = find_fp(filename);
	printf("Found %s at location %d\n", filename, cur_index);
	int i;
	if (cur_index==-1) {
		//Didn't find it
		//Close index 99, if it's valid
		//move everyone down
		if (fp_cache[CACHE_SIZE-1].fp!=NULL) {
			fflush(fp_cache[CACHE_SIZE-1].fp);
			fclose(fp_cache[CACHE_SIZE-1].fp);
		}
		for (i=CACHE_SIZE-2; i>=0; i--) {
			fp_cache[i+1] = fp_cache[i];
		}
	        //Put our new entry in 0
	        strcpy(fp_cache[0].filename, filename);
	        printf("Opening output file %s.\n", fp_cache[0].filename);
	        fp_cache[0].fp = fopfile(fp_cache[0].filename, "a");
		//Move the pointer to the end
		fseek(fp_cache[0].fp, 0, SEEK_END);
	} else if (cur_index>0) {
		//We found it, but not at the top
		fp_cache_entry tmp_cache;
		tmp_cache.fp = fp_cache[cur_index].fp;
		strcpy(tmp_cache.filename, fp_cache[cur_index].filename);
		//Move some people down, but nothing's being evicted
		for (i=cur_index-1; i>=0; i--) {
			printf("Moving %s from %d to %d.\n", fp_cache[i].filename, i, i+1);
			fflush(stdout);
			fp_cache[i+1] = fp_cache[i];
		}
		//Put the one we want in position 0
		fp_cache[0].fp = tmp_cache.fp;
		strcpy(fp_cache[0].filename, tmp_cache.filename);
	}
	//If cur_index==0, that's fine; it's already at the top of the list
	return fp_cache[0].fp;
}

void fsync_and_close(char* filename) {
	int index = find_fp(filename);
	if (index==-1) {
		//This filename has already been purged, don't worry about it
		return;
	}
	int fd = fileno(fp_cache[index].fp);
	fsync(fd);
	//Close the file pointer, move everyone up
	fclose(fp_cache[index].fp);
	//Clear entry
	fp_cache[index].fp = NULL;
	fp_cache[index].filename[0] = '\0';
	int i = index+1;
	while (i<CACHE_SIZE && fp_cache[i].fp!=NULL) {
		fp_cache[i-1] = fp_cache[i];
		i++;
	}
	//Set the last one to empty
	fp_cache[i-1].fp = NULL;
	fp_cache[i-1].filename[0] = '\0';
}


int find_fp(char* filename) {
	int i;
	for (i=0; i<CACHE_SIZE; i++) {
		if (fp_cache[i].fp==NULL) {
			//All the ones after this won't be active, either
			break;
		} else if (strcmp(fp_cache[i].filename,filename)==0) {
			return i;
		}
	}
	return -1;
}

void remove_fp_cache() {
	int i;
	for (i=0; i<CACHE_SIZE; i++) {
		if (fp_cache[i].fp!=NULL) {
			fflush(fp_cache[i].fp);
			fclose(fp_cache[i].fp);
		} else {
			//The rest of them will be empty
			break;
		}
	}
	free(fp_cache);
}
