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

void create_fp_cache(fp_cache_entry** fp_cache) {
	*fp_cache = check_malloc(sizeof(fp_cache_entry)*CACHE_SIZE);
	int i;
	for (i=0; i<CACHE_SIZE; i++) {
		(*fp_cache)[i].filename[0] = '\0';
		(*fp_cache)[i].fp = NULL;
	}
}

void find_and_use_fp(fp_cache_entry* fp_cache, char* filename) {
	//Locate filename in the cache
	//If it's not present, create a new entry
	//Move filename entry to index 0
	//Adjust everyone else accordingly
	int cur_index = find_fp(fp_cache, filename);
	int i;
	if (cur_index==0) {
		return;
	} else if (cur_index==-1) {
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
	} else {
		//Move some people down, but nothing's being evicted
		for (i=cur_index-1; i>=0; i++) {
			fp_cache[i+1] = fp_cache[i];
		}
	}
	//Put our new entry in 0
	strcpy(fp_cache[0].filename, filename);
	fp_cache[0].fp = fopen(fp_cache[0].filename, "a");
}


int find_fp(fp_cache_entry* fp_cache, char* filename) {
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

void remove_fp_cache(fp_cache_entry** fp_cache) {
	int i;
	for (i=0; i<CACHE_SIZE; i++) {
		if ((*fp_cache)[i].fp!=NULL) {
			fflush((*fp_cache)[i].fp);
			fclose((*fp_cache)[i].fp);
		}
	}
	free(*fp_cache);
}
