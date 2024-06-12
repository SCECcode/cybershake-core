#!/bin/bash

#gcc -o test test.c -I../include/ -I/../../../../software/Getpar/getpar/include -L../lib -lrupgen -L../../../../software/Getpar/getpar/lib -lget -I../../../../utils/libmemcached_1.0.18/include -lm -L${OLCF_FFTW_ROOT}/lib -lfftw3f -lm -L../../../../utils/libmemcached_1.0.18/lib -lmemcached

gcc -o test_rvfrac test_rvfrac.c -I../include/ -I/../../../../software/Getpar/getpar/include -L../lib -lrupgen -L../../../../software/Getpar/getpar/lib -lget -I../../../../utils/libmemcached_1.0.18/include -lm -L${TACC_FFTW3_LIB} -lfftw3f -lm -L../../../../utils/libmemcached_1.0.18/lib -lmemcached

#gcc -o test_opt test.c -I../src_opt/ -I/../../../../software/Getpar/getpar/include -L../src_opt/ -lrupgen -L../../../../software/Getpar/getpar/lib -lget -I../../../../utils/libmemcached_1.0.18/include -lm -L${OLCF_FFTW_ROOT}/lib -lfftw3f -lm -L../../../../utils/libmemcached_1.0.18/lib -lmemcached

#gcc -o test_opt_omp test.c -I../src_opt/ -I/../../../../software/Getpar/getpar/include -L../src_opt/ -lrupgen -L../../../../software/Getpar/getpar/lib -lget -I../../../../utils/libmemcached_1.0.18/include -lm -L${OLCF_FFTW_ROOT}/lib -lfftw3f_omp -lfftw3f -lm -L../../../../utils/libmemcached_1.0.18/lib -lmemcached -fopenmp

#gcc -o test_with_seed test_with_seed.c -I../include/ -I/../../../../software/Getpar/getpar/include -L../lib -lrupgen -L../../../../software/Getpar/getpar/lib -lget -I../../../../utils/libmemcached_1.0.18/include -lm -L${OLCF_FFTW_ROOT}/lib -lfftw3f -lm -L../../../../utils/libmemcached_1.0.18/lib -lmemcached

#gcc -o test_with_seed_gsf test_with_seed_gsf.c -I../include/ -I/../../../../software/Getpar/getpar/include -L../lib -lrupgen -L../../../../software/Getpar/getpar/lib -lget -I../../../../utils/libmemcached_1.0.18/include -lm -L${OLCF_FFTW_ROOT}/lib -lfftw3f -lm -L../../../../utils/libmemcached_1.0.18/lib -lmemcached

