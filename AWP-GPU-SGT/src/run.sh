#!/bin/bash
mpiexec -np 4 ./pmcl3d --NX 448 --NY 448 -Z 1024 -x 2 -y 2 --NSKPX 2 --NSKPY 2 --NBGZ 1 --NEDZ 1 --INSRC FAULTPOW --MEDIASTART 0 --NSRC 1 --NST 91 --READ_STEP 91
