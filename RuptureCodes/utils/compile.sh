#!/bin/bash

gcc -O3 -o extract_rup_geom extract_rup_geom.c -lrupgen -L/work2/00349/scottcal/frontera/CyberShake/software/RuptureCodes/RupGen-api-5.5.2/lib/ -lm -lfftw3f -L$TACC_FFTW3_LIB -lmemcached -L/work2/00349/scottcal/frontera/CyberShake/utils/libmemcached_1.0.18/lib

gcc -O3 -o extract_rup_geom_from_gsf extract_rup_geom_from_gsf.c -lrupgen -L/work2/00349/scottcal/frontera/CyberShake/software/RuptureCodes/RupGen-api-5.5.2/lib/ -lm -lfftw3f -L$TACC_FFTW3_LIB -lmemcached -L/work2/00349/scottcal/frontera/CyberShake/utils/libmemcached_1.0.18/lib

