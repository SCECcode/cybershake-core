#!/bin/csh

if ($# != 12) then
  echo "Usage:  single_exe.csh <site> <model_cords file> <nx> <ny> <nz> <CS_PATH> <SCRATCH_PATH> <LOG_ROOT> <MODELS> <MPI_CMD> <JOB_ID> <format>"
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
set FORMAT = $12

#set BPATH = `pwd`
set BPATH = ${CS_PATH}/UCVM
set BPROG = bin/ucvm-single-mpi
set MODELDIR = ${BPATH}/DataFiles

#set OUTDIR = ../../data/LayerBin/${SITE}
set OUTDIR = `pwd`
set SCRATCH_OUTDIR = ${SCRATCH_PATH}/job-${JOB_ID}/data/LayerBin
set MAIN_RUNDIR = ./

umask 002

mkdir -p ${OUTDIR}
#\mkdir -p ${SCRATCH_OUTDIR}
#ls -lt ${SCRATCH_OUTDIR}

#set FILEROOT = v4_sgt
set FILEROOT = v_sgt-${SITE}
if ($FORMAT == "awp") then
	set FILEROOT = awp.${SITE}.media
endif
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
	#if ($FORMAT == "rwg") then
	#	while (`expr $NY \* $NZ % $NP` != 0)
	#		@ NP = $NP - 1
	#	end
	#else if ($FORMAT == "awp") then
	#	while (`expr $NX \* $NZ % $NP` != 0)
	#		@ NP = $NP - 1
	#	end
	#endif
	if ${MODELS} == "cvmh" then
	        set MPI_CMD = "${MPI_CMD} -n ${NP}"
	else
	        set MPI_CMD = "${MPI_CMD} -n ${NP}"
	endif
endif

set VPMIN = 1700.0
set VSMIN = 500.0
set DENMIN = 1700.0

#Add to LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /work/00940/tera3d/CyberShake/software/UCVM/ucvm_12.2.0/lib:/work/00940/tera3d/CyberShake/software/UCVM/cvmh_11.9.1/lib:/work/00940/tera3d/CyberShake/software/UCVM/cvms/lib:/work/00940/tera3d/CyberShake/software/UCVM/euclid3-1.3/libsrc:/work/00940/tera3d/CyberShake/software/UCVM/proj_4.7.0/lib

echo "${MPI_CMD} ${BPATH}/${BPROG} nx=${NX} ny=${NY} nz=${NZ} cordfile=${MODELCORDS} depfile=${DEPFILE} modeldir=${MODELDIR} outfile=${OUTDIR}/${FILEROOT} models=${MODELS} min_vp=${VPMIN} min_vs=${VSMIN} min_rho=${DENMIN} format=${FORMAT} logdir=${LOGDIR}"

${MPI_CMD} ${BPATH}/${BPROG} nx=${NX} ny=${NY} nz=${NZ} cordfile=${MODELCORDS} depfile=${DEPFILE} modeldir=${MODELDIR} outfile=${OUTDIR}/${FILEROOT} models=${MODELS} min_vp=${VPMIN} min_vs=${VSMIN} min_rho=${DENMIN} format=${FORMAT} logdir=${LOGDIR}
endif

set RC = $?

if (${RC} == 0) then
  echo "Station= ${SITE}: Done   " `date`
else
  echo "Station= ${SITE}: FAILED WITH CODE ${RC}!   " `date`
  exit 1
endif

if ($FORMAT == "awp") then
	#Move to awp.<site>.media
	mv ${OUTDIR}/${FILEROOT} ${OUTDIR}/awp.${SITE}.media
endif

#echo "\\cp -r ${SCRATCH_OUTDIR}/* ${OUTDIR}"
#\cp -r ${SCRATCH_OUTDIR}/* ${OUTDIR}

