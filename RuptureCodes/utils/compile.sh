#!/bin/bash

gcc -O3 -o extract_rup_geom extract_rup_geom.c -lrupgen -L/gpfs/alpine/proj-shared/geo112/CyberShake/software/RuptureCodes/RupGen-api-5.4.2/lib/ -lm -lfftw3f -lmemcached -L/gpfs/alpine/proj-shared/geo112/CyberShake/utils/libmemcached_1.0.18/lib

gcc -O3 -o extract_rup_geom_from_gsf extract_rup_geom_from_gsf.c -lrupgen -L/gpfs/alpine/proj-shared/geo112/CyberShake/software/RuptureCodes/RupGen-api-5.4.2/lib/ -lm -lfftw3f -lmemcached -L/gpfs/alpine/proj-shared/geo112/CyberShake/utils/libmemcached_1.0.18/lib

