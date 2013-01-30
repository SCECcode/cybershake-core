#!/bin/csh

if ($# != 11) then
  echo "Usage:  genmod.csh <site> <model_cords file> <nx> <ny> <nz> <CS_PATH> <SCRATCH_PATH> <LOG_ROOT> <MODELS> <MPI_CMD> <JOB_ID>"
  exit 1
endif

set SITE = $1
set MODELCORDS = $2
set NX = $3
set NY = $4
set NZ = $5
set CS_PATH = $6
set SCRATCH_PATH = $7
set LOG_ROOT = $8
set MODELS = $9
set MPI_CMD = $10
set JOB_ID = $11

#set BPATH = `pwd`
set BPATH = ${CS_PATH}/UCVM
set BPROG = bin/ucvm-mpi
set MODELDIR = ${BPATH}/DataFiles

#set OUTDIR = ../../data/LayerBin/${SITE}
set OUTDIR = `pwd`
set SCRATCH_OUTDIR = ${SCRATCH_PATH}/job-${JOB_ID}/data/LayerBin
set MAIN_RUNDIR = ./

umask 002

mkdir -p ${OUTDIR}
\mkdir -p ${SCRATCH_OUTDIR}
ls -lt ${SCRATCH_OUTDIR}

set FILEROOT = v4_sgt
#set LOGDIR = ../../logs/GenLog/${SITE}
set LOGDIR = ${LOG_ROOT}/GenLog/${SITE}

#set MODELCORDS = ../../data/ModelParams/${SITE}/model_coords_GC_${SITE}
#set DEPFILE = ../../data/ModelParams/${SITE}/zdepths.meter
set DEPFILE = zdepths.meter

if ($MPI_CMD == "mpirun") then
	set NP = `wc -l < ${PBS_NODEFILE}`
	set MPI_CMD = "${MPI_CMD} -np ${NP} -machinefile ${PBS_NODEFILE}"
else if ($MPI_CMD == "aprun") then
	@ NP = $PBS_NUM_NODES * $PBS_NUM_PPN
	if ${MODELS} == "cvmh" then
	        set MPI_CMD = "${MPI_CMD} -n ${NP}"
	else
	        set MPI_CMD = "${MPI_CMD} -n ${NP}"
	endif
endif

echo "${MPI_CMD} ${BPATH}/${BPROG} nx=${NX} ny=${NY} nz=${NZ} cordfile=${MODELCORDS} depfile=${DEPFILE} modeldir=${MODELDIR} outdir=${SCRATCH_OUTDIR} fileroot=${FILEROOT} logdir=${LOGDIR} models=${MODELS}"

${MPI_CMD} ${BPATH}/${BPROG} nx=${NX} ny=${NY} nz=${NZ} cordfile=${MODELCORDS} depfile=${DEPFILE} modeldir=${MODELDIR} outdir=${SCRATCH_OUTDIR} fileroot=${FILEROOT} logdir=${LOGDIR} models=${MODELS}

set RC = $?

if (${RC} == 0) then
  echo "Station= ${SITE}: Done   " `date`
else
  echo "Station= ${SITE}: FAILED WITH CODE ${RC}!   " `date`
  exit 1
endif

\cp -r ${SCRATCH_OUTDIR}/* ${OUTDIR}

