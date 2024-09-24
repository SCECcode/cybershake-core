#!/bin/bash

icc -o test test.c -I../include/ -I/../../../../software/Getpar/getpar/include -L../lib -lrupgen -L../../../../software/Getpar/getpar/lib -lget -I../../../../utils/libmemcached_1.0.18/include -lm -L${TACC_FFTW3_LIB} -lfftw3f -lm -L../../../../utils/libmemcached_1.0.18/lib -lmemcached

icc -o get_num_rvs get_num_rvs.c -I../include/ -I/../../../../software/Getpar/getpar/include -L../lib -lrupgen -L../../../../software/Getpar/getpar/lib -lget -I../../../../utils/libmemcached_1.0.18/include -lm -L${TACC_FFTW3_LIB} -lfftw3f -lm -L../../../../utils/libmemcached_1.0.18/lib -lmemcached

