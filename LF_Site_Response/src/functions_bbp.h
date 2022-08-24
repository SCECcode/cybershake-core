#ifndef FUNCTIONS_BBP_H
#define FUNCTIONS_BBP_H

#include "include.h"
#include "structure.h"
#include "WccFormat/Progs/function.h"

void *check_malloc(size_t len);
void *check_realloc(void *ptr,size_t len);

float ucvm_vsD(float lon, float lat, char* model, int depth);

#endif
