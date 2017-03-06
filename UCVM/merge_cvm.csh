#!/bin/csh

umask 002

if ($# != 10) then
  echo "Usage: $0 <site> <gridfile> <velocity p-file> <velocity s-file> <velocity d-file> <CS_PATH> <SCRATCH_PATH> <LOG_ROOT> <MPI_CMD> <JOB_ID>"
  exit 1
endif

set SITE = $1
set GRIDFILE = $2
set PMOD = $3
set SMOD = $4
set DMOD = $5
set CS_PATH = $6
set SCRATCH_PATH = $7
set LOG_ROOT = $8
set MPI_CMD = $9
set JOB_ID = $10

set BPATH = ${CS_PATH}/UCVM/bin
set BPROG = merge_cvm

#set LAYDIR = ${CS_PATH}/data/LayerBin/${SITE}
set LAYDIR = `pwd`
#set OUTDIR = ${CS_PATH}/data/Models/${SITE}
set OUTDIR = `pwd`
set LOGDIR = ${LOG_ROOT}/GenLog/${SITE}
#set LOGDIR = /home/rcf-104/CyberShake2007/logs/GenLog/${SITE}

#set FILELIST = ../../data/Models/${SITE}/mod.filelist
set FILELIST = mod.filelist

set SCRATCH_ROOT = ${SCRATCH_PATH}/job-$JOB_ID
set SCRATCH_LAYDIR = ${SCRATCH_ROOT}/data/LayerBin/${SITE}
set SCRATCH_OUTDIR = ${SCRATCH_ROOT}/data/Models/${SITE}

# Make all the output and working directories
mkdir -p ${OUTDIR} ${LOGDIR} ${SCRATCH_OUTDIR} ${SCRATCH_LAYDIR}

# Copy the input files to scratch
# Let's try symlinking instead
#cp -r ${LAYDIR}/* ${SCRATCH_LAYDIR}
foreach i ( ${LAYDIR}/* )
ln -s $i ${SCRATCH_LAYDIR}/`basename $i`
end

#set PMOD = v4_sgt-${SITE}.p
#set SMOD = v4_sgt-${SITE}.s
#set DMOD = v4_sgt-${SITE}.d

set VPMAX = 16.7
set VSMAX = 13.9
set DENMAX = 12.9

set VPMIN = 1.7
set VSMIN = 0.5
set DENMIN = 1.7

# Construct list of input files
ls ${SCRATCH_LAYDIR}/v4_sgt-*.cvm > $FILELIST

if ($MPI_CMD == "mpirun") then
        set NP = `wc -l < $PBS_NODEFILE`
        set MPI_CMD = "${MPI_CMD} -np $NP -machinefile $PBS_NODEFILE"
else if ($MPI_CMD == "aprun") then
	@ NP = $PBS_NUM_NODES * $PBS_NUM_PPN
        set MPI_CMD = "${MPI_CMD} -n $NP"
endif

# Merge input files
echo "$MPI_CMD $BPATH/$BPROG gridfile=$GRIDFILE filelist=$FILELIST outdir=$SCRATCH_OUTDIR pmodfile=$PMOD smodfile=$SMOD dmodfile=$DMOD pmin=$VPMIN smin=$VSMIN dmin=$DENMIN pmax=$VPMAX smax=$VSMAX dmax=$DENMAX fdversion=1.14 logdir=$LOGDIR check_poisson_ratio=1"


$MPI_CMD $BPATH/$BPROG gridfile=$GRIDFILE filelist=$FILELIST outdir=$SCRATCH_OUTDIR pmodfile=$PMOD smodfile=$SMOD dmodfile=$DMOD pmin=$VPMIN smin=$VSMIN dmin=$DENMIN pmax=$VPMAX smax=$VSMAX dmax=$DENMAX fdversion=1.14 logdir=$LOGDIR check_poisson_ratio=1


set RC = $?
if ($RC == 0) then
  echo "Merge for ${SITE}: Done   " `date`
else
  echo "Merge for ${SITE}: FAILED WITH CODE ${RC}!   " `date`
  exit 1
endif

# Copy the output from temp storage
cp -r $SCRATCH_OUTDIR/* $OUTDIR

# Give group read and write to the output
chmod -R g+rw $OUTDIR

# Clean up temp files
rm -rf ${SCRATCH_ROOT}

echo "Done   " `date`

exit 0
