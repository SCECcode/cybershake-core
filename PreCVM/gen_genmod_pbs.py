#!/usr/bin/python

def genDepFile("zdepths.meter", gridout):
    KM2METER = 1000.0
    input = open("../../data/ModelParams/gridout_" + site)
    iz = 0
    go = 0
    for line in input:
        split...
        
        
        
        
        
        
        
        
        
    BEGIN{iz=0;go=0;} {
    split($1,a,"=");
    if(a[1]=="nz") {
        nz=a[2];
        go=1;
    } else if(go==1 && iz<nz) {
        iz++;
        dep[iz]=$2*cnv;
        if(dep[iz]<0.0)
            dep[iz]=0.0;
    }
} END {
    for(iz=1;iz<=nz;iz++) {
        printf "%14.2f\n",dep[iz];
    }
}


if len(sys.argv) < 2:
    print "Usage:  gen_genmod_pbs.py <site>"
    print "Example: gen_genmod_pbs.py USC"
    exit(1)

site = sys.argv[1]
output = open("../V4-WrapC/genmod." + site + ".pbs", "w")
output.writeLine('#!/bin/csh')
output.writeLine('#*** The "#PBS" lines must come before any non-blank non-comment lines ***')
output.writeLine('#PBS -e genmod.' + site + '.e -o genmod.' + site + '.o')
output.writeLine('#PBS -l walltime=05:00:00,nodes=50:myri:ppn=2')
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
output.writeLilne('  set NP = `wc -l < $PBS_NODEFILE `')
output.writeLine('  echo "Running on $NP processors: "`cat $PBS_NODEFILE`')
output.writeLine('else')
output.writeLine('  echo "This job did not have the number of nodes specified with"')
output.writeLine('  echo "the node= resource"')
output.writeLine('  exit 1')
output.writeLine('endif')
output.writeLine('\n#The following should contain your program and any arguments')
output.writeLine('set BPATH = `pwd`')
output.writeLine('set BPROG = bin/ver4C-mpi')
output.writeLine('set MODELDIR = ${BPATH}/DataFiles')
output.writeLine('\nset NX = 
output.writeLine('set NY = 
output.writeLine('set NZ = 
output.writeLine('\nset OUTDIR = ../../data/LayerBin')
output.writeLine('set SCRATCH_OUTDIR = /scratch/pbsjob-${PBS_JOBID}/data/LayerBin')
output.writeLine('set MAIN_RUNDIR = ./')
output.writeLine('\mkdir -p ${SCRATCH_OUTDIR}')
output.writeLine('ls -lt $SCRATCH_OUTDIR')
output.writeLine('\nset FILEROOT = v4_sgt')
output.writeLine('set LOGDIR = ../../logs/GenLog')
output.writeLine('\nset MODELCORDS = ../../data/ModelParams/' + site + '/model_coords_GC_' + site)
output.writeLine('set DEPFILE = zdepths.meter')

gridout = '../../data/ModelParams/' + site + '/gridout_' + site
genDepFile("zdepths.meter", gridout)

gawk -v cnv=$KM2METER 'BEGIN{iz=0;go=0;}{ \
split($1,a,"="); \
if(a[1]=="nz"){nz=a[2];go=1;} \
else if(go==1 && iz<nz){iz++;dep[iz]=$2*cnv; if(dep[iz]<0.0)dep[iz]=0.0;}} \
END { \
for(iz=1;iz<=nz;iz++){printf "%14.2f\n",dep[iz];}}' \
$GRIDOUT > $DEPFILE

output.writeLine('mpiexec ${BPATH}/${BPROG} \ ')
output.writeLine('         nx=$NX ny=$NY nz=$NZ \ ')
output.writeLine('         cordfile=$MODELCORDS \ ')
output.writeLine('         depfile=$DEPFILE \ ')
output.writeLine('         modeldir=$MODELDIR \ ')
output.writeLine('         outdir=$SCRATCH_OUTDIR \ ')
output.writeLine('         fileroot=$FILEROOT \ ')
output.writeLine('         logdir=$LOGDIR')
output.writeLine('\n\cp -r $SCRATCH_OUTDIR $OUTDIR')
output.writeLine('\necho "Done   " `date`')
output.writeLine('\nexit')




BEGIN{iz=0;go=0;} {
    split($1,a,"=");
    if(a[1]=="nz") {
        nz=a[2];
        go=1;
    } else if(go==1 && iz<nz) {
        iz++;
        dep[iz]=$2*cnv;
        if(dep[iz]<0.0)
            dep[iz]=0.0;
    }
} END {
    for(iz=1;iz<=nz;iz++) {
        printf "%14.2f\n",dep[iz];
    }
}