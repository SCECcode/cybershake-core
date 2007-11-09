#!/usr/bin/python

import sys

if len(sys.argv) < 2:
    print "Usage: gen_merge_cvm.py <site>"
    print "Example: gen_merge_cvm.py USC"
    sys.exit(1)

site = sys.argv[1]

output = open("../V4-WrapC/PBS/merge_cvm." + site + ".pbs", "w")
output.write('#!/bin/csh\n')

output.write('\n#*** The "#PBS" lines must come before any non-blank non-comment lines ***\n')

output.write('\n#PBS -e merge_cvm.' + site + '.e -o merge_cvm.' + site + '.o\n')
output.write('#PBS -l walltime=01:00:00,nodes=2:myri:ppn=2\n')
output.write('#PBS -A lc_pjm\n')

output.write('\nif ($?PBS_JOBID) then           # if this is a PBS job\n')
output.write('  echo "Starting" $PBS_JOBID `date`\n')
output.write('  echo "Initiated on `hostname`"\n')
output.write('  echo ""\n')
output.write('  cd "$PBS_O_WORKDIR"           # connect to working directory of qsub\n')
output.write('else\n')
output.write('  echo "This script must be run as a PBS job"\n')
output.write('  exit 1\n')
output.write('endif\n')

output.write('\nif ($?PBS_NODEFILE) then\n')
output.write('  #count the number of processors assigned by PBS\n')
output.write('  set NP = `wc -l < $PBS_NODEFILE `\n')
output.write('  echo "Running on $NP processors: "`cat $PBS_NODEFILE`\n')
output.write('else\n')
output.write('  echo "This job did not have the number of nodes specified with"\n')
output.write('  echo "the node= resource"\n')
output.write('  exit 1\n')
output.write('endif\n')

output.write('\n#The following should contain your program and any arguments\n')

output.write('\nset BPATH = `pwd`\n')
output.write('set BPROG = bin/merge_cvm\n')

output.write('\nset SITE = %s\n' % site)

output.write('\nset LAYDIR = ../../data/LayerBin/${SITE}\n')
output.write('set OUTDIR = ../../data/Models/${SITE}\n')

output.write('\nset SCRATCH_ROOT = /scratch/pbsjob-${PBS_JOBID}\n')
output.write('set SCRATCH_LAYDIR = ${SCRATCH_ROOT}/data/LayerBin/${SITE}\n')
output.write('set SCRATCH_OUTDIR = ${SCRATCH_ROOT}/data/Models/${SITE}\n')

output.write('\nmkdir -p ${OUTDIR}\n')
output.write('\mkdir -p ${SCRATCH_OUTDIR}\n')
output.write('mkdir -p ${SCRATCH_LAYDIR}\n')
output.write('ls -lt $SCRATCH_OUTDIR\n')

output.write('\ncp -r ${LAYDIR}/* ${SCRATCH_LAYDIR}\n')
output.write('ls -lt $SCRATCH_LAYDIR\n')

output.write('\nset FILEROOT = v4_sgt\n')
output.write('set LOGDIR = ../../logs/GenLog/${SITE}\n')

output.write('\nset FILELIST = ../../data/Models/${SITE}/mod.filelist\n')
output.write('set GRIDFILE = ../../data/ModelParams/${SITE}/gridfile-${SITE}\n')

output.write('\nset PMOD = v4_sgt-${SITE}.p\n')
output.write('set SMOD = v4_sgt-${SITE}.s\n')
output.write('set DMOD = v4_sgt-${SITE}.d\n')

output.write('\nset VPMAX = 16.7\n')
output.write('set VSMAX = 13.9\n')
output.write('set DENMAX = 12.9\n')

output.write('\nset VPMIN = 1.7\n')
output.write('set VSMIN = 0.5\n')
output.write('set DENMIN = 1.7\n')
output.write('\nset inroot = ${SCRATCH_LAYDIR}/${FILEROOT}-\n')

#This is not python-ed away because it has to be done at runtime, since $PBS_JOBID is not known until runtime.
output.write('''\nset PROCS = `ls ${inroot}*.cvm | gawk -F"/" '{n=split($NF,a,".");m=split(a[n-1],b,"-");print b[m];}' `\n''')

output.write('\n\\rm $FILELIST\n')
output.write('foreach proc ( $PROCS )\n')

output.write('\necho "${inroot}${proc}.cvm" >> $FILELIST\n')

output.write('\nend\n')

output.write('\nmpiexec $BPATH/$BPROG \ \n')
output.write('         gridfile=$GRIDFILE filelist=$FILELIST \ \n')
output.write('         outdir=$SCRATCH_OUTDIR pmodfile=$PMOD smodfile=$SMOD dmodfile=$DMOD \ \n')
output.write('         pmin=$VPMIN smin=$VSMIN dmin=$DENMIN \ \n')
output.write('         pmax=$VPMAX smax=$VSMAX dmax=$DENMAX \ \n')
output.write('         fdversion=1.14 logdir=$LOGDIR \n')

output.write('\n\cp -r $SCRATCH_OUTDIR/* $OUTDIR\n')

output.write('\necho "Done   " `date`\n')
output.write('\nexit\n')
