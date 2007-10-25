#!/usr/bin/python

if len(sys.argv) < 2
    print "Usage: gen_merge_cvm.py <site>"
    print "Example: gen_merge_cvm.py USC"
    exit(1)

site = sys.argv[1]
output = open("../V4-WrapC/merge_cvm." + site + ".pbs", "w")
output.writeLine('#!/bin/csh')
output.writeLine('#*** The "#PBS" lines must come before any non-blank non-comment lines ***')
output.writeLine('#PBS -e merge_cvm.' + site + '.e -o merge_cvm.' + site + '.o')
output.writeLine('#PBS -l walltime=01:00:00,nodes=2:myri:ppn=2')
output.writeLine('\nif ($?PBS_JOBID) then           # if this is a PBS job')
output.writeLine('  echo "Starting" $PBS_JOBID `date`')
output.writeLine('  echo "Initiated on `hostname`"')
output.writeLine('  echo ""')
output.writeLine('  cd "$PBS_O_WORKDIR"           # connect to working directory of qsub')
output.writeLine('else')
output.writeLine('  echo "This script must be run as a PBS job"')
output.writeLine('  exit 1')
output.writeLine('endif')
output.writeLine('\nif ($?PBS_NODEFILE) then')
output.writeLine('  #count the number of processors assigned by PBS')
output.writeLine('  set NP = `wc -l < $PBS_NODEFILE `')
output.writeLine('  echo "Running on $NP processors: "`cat $PBS_NODEFILE`')
output.writeLine('else')
output.writeLine('  echo "This job did not have the number of nodes specified with"')
output.writeLine('  echo "the node= resource"')
output.writeLine('  exit 1')
output.writeLine('endif')
output.writeLine('\n#The following should contain your program and any arguments')
set BPATH = `pwd`
set BPROG = bin/merge_cvm

set LAYDIR = LayerBin
set OUTDIR = Models

set SCRATCH_ROOT = /scratch/pbsjob-${PBS_JOBID}
set SCRATCH_LAYDIR = ${SCRATCH_ROOT}/${LAYDIR}
set SCRATCH_OUTDIR = ${SCRATCH_ROOT}/${OUTDIR}

set MAIN_LAYDIR = ./$LAYDIR
set MAIN_MODDIR = ../

set SCRATCH_LAYDIR = ${MAIN_LAYDIR}
set SCRATCH_OUTDIR = ${MAIN_MODDIR}/${OUTDIR}

\mkdir -p ${SCRATCH_OUTDIR}
ls -lt $SCRATCH_OUTDIR

#\cp -r ${MAIN_LAYDIR} ${SCRATCH_ROOT}
ls -lt $SCRATCH_LAYDIR

set FILEROOT = v4_sgt
set LOGDIR = GenLog

set FILELIST = mod.filelist
set GRIDFILE = ../ModelParams/gridfile-h0.2

set PMOD = v4_sgt-h0.2.p
set SMOD = v4_sgt-h0.2.s
set DMOD = v4_sgt-h0.2.d

output.writeLine('set VPMAX = 16.7')
output.writeLine('set VSMAX = 13.9')
output.writeLine('set DENMAX = 12.9')
output.writeLine('\nset VPMIN = 1.7')
output.writeLine('set VSMIN = 0.5')
output.writeLine('set DENMIN = 1.7')
output.writeLine('\nset inroot = ${SCRATCH_LAYDIR}/${FILEROOT}-')

set PROCS = `ls ${inroot}*.cvm | gawk -F"/" '{n=split($NF,a,".");m=split(a[n-1],b,"-");print b[m];}' `

gawk -F"/" {
    n=split($NF,a,".");
    m=split(a[n-1],b,"-");
    print b[m];
} `


output.writeLine('\rm $FILELIST')
foreach proc ( $PROCS )

echo "${inroot}${proc}.cvm" >> $FILELIST

end

output.writeLine('mpiexec $BPATH/$BPROG \ ')
output.writeLine('         gridfile=$GRIDFILE filelist=$FILELIST \ ')
output.writeLine('         outdir=$SCRATCH_OUTDIR pmodfile=$PMOD smodfile=$SMOD dmodfile=$DMOD \ ')
output.writeLine('         pmin=$VPMIN smin=$VSMIN dmin=$DENMIN \ ')
output.writeLine('         pmax=$VPMAX smax=$VSMAX dmax=$DENMAX \ ')
output.writeLine('         fdversion=1.14 logdir=$LOGDIR ')

#\cp -r $SCRATCH_OUTDIR $MAIN_MODDIR

output.writeLine('\necho "Done   " `date`')
output.writeLine('\nexit')
