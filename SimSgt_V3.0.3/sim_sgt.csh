#!/bin/csh

# Make output files group writeable
umask 002

# Check parameters
if ($# != 15) then
  echo "Usage: $0 SITE NX NY NZ MODEL_LAT MODEL_LON MODEL_ROT XSRC YSRC RUN CS_PATH SCRATCH_PATH TMP_PATH MPI_CMD JOB_ID"
  echo
  echo "SITE: The site ID (e.g. LADT)"
  echo "NX:"
  echo "NY:"
  echo "NZ:"
  echo "MODEL_LAT:"
  echo "MODEL_LON:"
  echo "MODEL_ROT:"
  echo "XSRC: the x grid coordinate of the site"
  echo "YSRC: the y grid coordinate of the site"
  echo "RUN: The run to do ('x', 'y', or 'z')"
  exit 1
endif

# Read parameters
set SITE = $1
set NX = $2
set NY = $3
set NZ = $4
set MODEL_LAT = $5
set MODEL_LON = $6
set MODEL_ROT = $7
set XSRC = $8
set YSRC = $9
set RUN = $10
set CS_PATH = $11
set SCRATCH_PATH = $12
set TMP_PATH = $13
set MPI_CMD = $14
set JOB_ID = $15

set PWD = `pwd`

# Determine which run to do
switch ($RUN)
  case x:
    set R = 1
    breaksw
  case y:
    set R = 2
    breaksw
  case z:
    set R = 3
    breaksw
  default:
    echo "Invalid run: $RUN"
    exit 1
    breaksw
endsw

#set RUNS = ( fx fy fz )
set RUNS = ( fx fy )
set XMOM = ( 1.0e+20 0.0 0.0 )
set YMOM = ( 0.0 1.0e+20 0.0 )
set ZMOM = ( 0.0 0.0 1.0e+20 )

set H = 0.2
set NT = 20000
#set NT = 1000
set DT = 0.01

set BIN_DIR = ${CS_PATH}/SimSgt_V3.0.3/bin

#set SRCLOCS_DIR = ../../data/SgtInfo/FdLocs
#set SGTLOCS_DIR = ${CS_PATH}/data/SgtInfo/SgtCords
set SGTLOCS_DIR = $PWD
#set LOG_DIR = ../../logs/SgtLog/${SITE}
set LOG_DIR = ${CS_PATH}/logs/SgtLog_V3/$SITE-${JOB_ID}
#set MAIN_OUTPDIR = ${CS_PATH}/data/OutBin/${SITE}
set MAIN_OUTPDIR = $PWD
#set VMODDIR = ${CS_PATH}/data/Models/${SITE}
set VMODDIR = $PWD

#set TMP_ROOT = ${TMP_PATH}/job-${JOB_ID}
set TMP_ROOT = $PWD
#set SCRATCH_ROOT = ${SCRATCH_PATH}/job-${JOB_ID}
set SCRATCH_ROOT = $PWD

#set TMP_SEISDIR = ${TMP_ROOT}/OutBin/${SITE}
set TMP_SEISDIR = ${TMP_ROOT}
#set SCRATCH_VMODDIR = ${SCRATCH_ROOT}/Models/${SITE}
set SCRATCH_VMODDIR = ${SCRATCH_ROOT}

\mkdir -p ${LOG_DIR} ${MAIN_OUTPDIR} ${SCRATCH_ROOT} ${TMP_ROOT} \
	${SCRATCH_VMODDIR} ${TMP_SEISDIR}

# Copy VM files to scratch
#echo "copy command: cp ${VMODDIR}/v_sgt-* ${SCRATCH_VMODDIR}"
#ls ${VMODDIR}
#ls ${VMODDIR}/v_sgt-*
#\cp "${VMODDIR}/v_sgt-*" ${SCRATCH_VMODDIR}
#echo "Model files: ${SCRATCH_VMODDIR}"
#ls -lt $SCRATCH_VMODDIR
#set DIR = $VMODDIR
#echo "${DIR}: "
#ls -lt $DIR
#set DDIR = `dirname $DIR`
#while ("$DDIR" != "$DIR")
#	echo "${DDIR}: "
#	ls -lt $DDIR
#	set DIR = $DDIR
#	set DDIR = `dirname $DIR`
#end

# Read fdsrcloc params
#set XSRC = `gawk '{print $1;}' ${SRCLOCS_DIR}/${SITE}.fdloc `
#set YSRC = `gawk '{print $2;}' ${SRCLOCS_DIR}/${SITE}.fdloc `

set NAME = ${SITE}-${RUNS[$R]}

if ($MPI_CMD == "mpirun") then
        set NP = `wc -l < $PBS_NODEFILE`
        set MPI_CMD = "${MPI_CMD} -np $NP -machinefile $PBS_NODEFILE"
else if ($MPI_CMD == "aprun") then
	@ NP = $PBS_NUM_NODES * $PBS_NUM_PPN
	set MPI_CMD = "${MPI_CMD} -n $NP"
endif

set NX_PROC = 25
set NY_PROC = 40
set NZ_PROC = 4

limit coredumpsize unlimited

#$MPI_CMD ${BIN_DIR}/emod3d-mpi par=${CS_PATH}/SimSgt/e3d-D.par name=${NAME} logdir=${LOG_DIR} pmodfile=v_sgt-${SITE}.p nx=${NX} ny=${NY} nz=${NZ} h=${H} nt=${NT} dt=${DT} modellat=${MODEL_LAT} modellon=${MODEL_LON} modelrot=${MODEL_ROT} smodfile=v_sgt-${SITE}.s dmodfile=v_sgt-${SITE}.d xmom=$XMOM[$R] ymom=$YMOM[$R] zmom=$ZMOM[$R] xsrc=$XSRC ysrc=$YSRC zsrc=1 sgtout=1 sgt_tinc=10 sgtcords=${SGTLOCS_DIR}/${SITE}.cordfile seisdir=${TMP_SEISDIR} main_dump_dir=${MAIN_OUTPDIR} vmoddir=${SCRATCH_VMODDIR}

echo "$MPI_CMD ${BIN_DIR}/emod3d-mpi par=${CS_PATH}/SimSgt_V3.0.3/e3d-D.par name=${NAME} logdir=${LOG_DIR} pmodfile=v_sgt-${SITE}.p nx=${NX} ny=${NY} nz=${NZ} h=${H} nt=${NT} dt=${DT} modellat=${MODEL_LAT} modellon=${MODEL_LON} modelrot=${MODEL_ROT} smodfile=v_sgt-${SITE}.s dmodfile=v_sgt-${SITE}.d xmom=$XMOM[$R] ymom=$YMOM[$R] zmom=$ZMOM[$R] xsrc=$XSRC ysrc=$YSRC zsrc=1 sgtout=1 sgt_tinc=10 sgtcords=${SGTLOCS_DIR}/${SITE}.cordfile seisdir=${TMP_SEISDIR} main_dump_dir=${MAIN_OUTPDIR} vmoddir=${SCRATCH_VMODDIR} nproc_x=${NX_PROC} nproc_y=${NY_PROC} nproc_z=${NZ_PROC}"

$MPI_CMD ${BIN_DIR}/emod3d-mpi par=${CS_PATH}/SimSgt_V3.0.3/e3d-D.par name=${NAME} logdir=${LOG_DIR} pmodfile=v_sgt-${SITE}.p nx=${NX} ny=${NY} nz=${NZ} h=${H} nt=${NT} dt=${DT} modellat=${MODEL_LAT} modellon=${MODEL_LON} modelrot=${MODEL_ROT} smodfile=v_sgt-${SITE}.s dmodfile=v_sgt-${SITE}.d xmom=$XMOM[$R] ymom=$YMOM[$R] zmom=$ZMOM[$R] xsrc=$XSRC ysrc=$YSRC zsrc=1 sgtout=1 sgt_tinc=10 sgtcords=${SGTLOCS_DIR}/${SITE}.cordfile seisdir=${TMP_SEISDIR} main_dump_dir=${MAIN_OUTPDIR} vmoddir=${SCRATCH_VMODDIR} nproc_x=${NX_PROC} nproc_y=${NY_PROC} nproc_z=${NZ_PROC}

# Check return code
set RC = $?
if ($RC == 0) then
  echo "Station= ${SITE}, Run= ${RUNS[$R]}: Done   " `date`
else
  echo "Station= ${SITE}, Run= ${RUNS[$R]}: FAILED WITH CODE ${RC}!   " `date`
  exit -1
endif

# Add group write permissions
chmod -R g+w ${MAIN_OUTPDIR}

# Remove temporary files
#rm -rf ${TMP_ROOT} ${SCRATCH_ROOT}

exit 0
