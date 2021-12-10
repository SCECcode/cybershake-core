#!/bin/csh

if ($# != 13 && $# != 14) then
  echo "Usage:  single_exe.csh <site> <model_cords file> <nx> <ny> <nz> <CS_PATH> <SCRATCH_PATH> <LOG_ROOT> <MODELS> <MPI_CMD> <JOB_ID> <format> <surface_cvm_depth> [min vs]"
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
set SURF_CVM_DEPTH = $13

if ($# == 14) then
	set VSMIN = $14
	if ($VSMIN == 900.0) then
		set VPMIN = 1800.0
		set DENMIN = 2000.0
	else
		set VPMIN = `echo "3.4 * $VSMIN" | bc`
		set DENMIN = $VPMIN
	endif
else
	set VSMIN = 500.0
	set VPMIN = 1700.0
	set DENMIN = 1700.0
endif

echo "VSMIN=$VSMIN, VPMIN=$VPMIN, DENMIN=$DENMIN"

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
	@ NP = $PBS_NUM_NODES * 8
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
	set MPI_CMD = "${MPI_CMD} -S 4"
else if ($MPI_CMD == "jsrun") then
	#We're running on Summit
	#LSB_DJOB_NUMPROC is 42*num_nodes + 1
	@ NP = $LSB_DJOB_NUMPROC - 1
	@ NNODES = $NP / 42
	@ NUM_RESOURCE_SETS = $NNODES * 42
	#-a 1: 1 MPI tasks per core
	#-c 1: 1 core per resource set
	#-r 42: 42 resource sets per node
	set MPI_CMD = "${MPI_CMD} -n $NUM_RESOURCE_SETS -a 1 -c 1 -r 42 -l CPU-CPU"
	echo "LSB_DJOB_NUMPROC=$LSB_DJOB_NUMPROC, NP=$NP, NNODES=$NNODES, NUM_RESOURCE_SETS=$NUM_RESOURCE_SETS"
endif

#set VPMIN = 1700.0
#set VSMIN = 500.0
#set DENMIN = 1700.0

#Add to LD_LIBRARY_PATH
#setenv LD_LIBRARY_PATH /work/00940/tera3d/CyberShake/software/UCVM/ucvm_12.2.0/lib:/work/00940/tera3d/CyberShake/software/UCVM/cvmh_11.9.1/lib:/work/00940/tera3d/CyberShake/software/UCVM/cvms/lib:/work/00940/tera3d/CyberShake/software/UCVM/euclid3-1.3/libsrc:/work/00940/tera3d/CyberShake/software/UCVM/proj_4.7.0/lib

set UCVM_HOME = /gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0
setenv LD_LIBRARY_PATH ${UCVM_HOME}/lib/euclid3/lib:${UCVM_HOME}/lib/proj-4/lib:${UCVM_HOME}/model/cvms426/lib:${UCVM_HOME}/model/cencal/lib:${UCVM_HOME}/model/cvms5/lib:${UCVM_HOME}/model/cca/lib:$LD_LIBRARY_PATH
setenv PROJ_LIB /gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/lib/proj-4/share/proj

#source /lustre/atlas/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/setup.csh
echo $LD_LIBRARY_PATH
#ldd ${BPATH}/${BPROG}

#On Summit, must unload darshan-runtime module or will hang on MPI_Finalize
module unload darshan-runtime

echo "${MPI_CMD} ${BPATH}/${BPROG} nx=${NX} ny=${NY} nz=${NZ} cordfile=${MODELCORDS} depfile=${DEPFILE} modeldir=${MODELDIR} outfile=${OUTDIR}/${FILEROOT} models=${MODELS} min_vp=${VPMIN} min_vs=${VSMIN} min_rho=${DENMIN} format=${FORMAT} logdir=${LOGDIR} surface_cvm_depth=${SURF_CVM_DEPTH}"

${MPI_CMD} ${BPATH}/${BPROG} nx=${NX} ny=${NY} nz=${NZ} cordfile=${MODELCORDS} depfile=${DEPFILE} modeldir=${MODELDIR} outfile=${OUTDIR}/${FILEROOT} models=${MODELS} min_vp=${VPMIN} min_vs=${VSMIN} min_rho=${DENMIN} format=${FORMAT} logdir=${LOGDIR} surface_cvm_depth=${SURF_CVM_DEPTH}

set RC = $?

if (${RC} == 0) then
  echo "Station= ${SITE}: Done   " `date`
else
  echo "Station= ${SITE}: FAILED WITH CODE ${RC}!   " `date`
  exit 1
endif

if ($FORMAT == "awp") then
	#Move to awp.<site>.media if it's different
	mv ${OUTDIR}/${FILEROOT} ${OUTDIR}/awp.${SITE}.media
endif

#echo "\\cp -r ${SCRATCH_OUTDIR}/* ${OUTDIR}"
#\cp -r ${SCRATCH_OUTDIR}/* ${OUTDIR}

