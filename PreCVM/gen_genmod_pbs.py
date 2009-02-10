#!/usr/bin/env python

import sys

def genDepFile(zdepth, gridout):
    '''this replicates the gawk code from RWG in genmod.XXX.pbs -- it writes a list of all depths from the gridout file to the zdepth file.  It returns the number of steps in the X, Y, and Z directions as a tuple.'''
    KM2METER = 1000.0
    input = open(gridout)
    output = open(zdepth, "w")
    inputContents = input.readlines()
    intNX = int((inputContents[1].split("="))[1])
    intNY = int((inputContents[1+intNX+2].split("="))[1])
    intNZ = int((inputContents[1+intNX+2+intNY+2].split("="))[1])
    for i in range(len(inputContents)-intNZ, len(inputContents)):
        #print float(inputContents[i].split()[1])
        output.write("%14.2f\n" % (float(inputContents[i].split()[1])*KM2METER))
    output.flush()
    output.close()
    return [intNX, intNY, intNZ]


if len(sys.argv) < 2:
    print "Usage:  gen_genmod_pbs.py <site>"
    print "Example: gen_genmod_pbs.py USC"
    sys.exit(1)

site = sys.argv[1]

zdepthFile = "../../data/ModelParams/" + site + "/zdepths.meter"
gridout = '../../data/ModelParams/' + site + '/gridout_' + site
numSteps = genDepFile(zdepthFile, gridout)

output = open("../V4-WrapC/PBS/genmod." + site + ".pbs", "w")
output.write('#!/bin/csh\n')

output.write('\n#*** The "#PBS" lines must come before any non-blank non-comment lines ***\n')

output.write('\n#PBS -e genmod.' + site + '.e -o genmod.' + site + '.o\n')
output.write('#PBS -l walltime=05:00:00,nodes=50:myri:ppn=2\n')
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
output.write('set BPATH = `pwd`\n')
output.write('set BPROG = bin/ver4C-mpi\n')
output.write('set MODELDIR = ${BPATH}/DataFiles\n')

output.write('set SITE = %s\n' % site)
output.write('set NX = %d\n' % numSteps[0])
output.write('set NY = %d\n' % numSteps[1])
output.write('set NZ = %d\n' % numSteps[2])

output.write('\nset OUTDIR = ../../data/LayerBin/${SITE}\n')
output.write('set SCRATCH_OUTDIR = /scratch/pbsjob-${PBS_JOBID}/data/LayerBin\n')
output.write('set MAIN_RUNDIR = ./\n')

output.write('\nmkdir -p ${OUTDIR}\n')
output.write('\mkdir -p ${SCRATCH_OUTDIR}\n')
output.write('ls -lt $SCRATCH_OUTDIR\n')

output.write('\nset FILEROOT = v4_sgt\n')
output.write('set LOGDIR = ../../logs/GenLog/${SITE}\n')

output.write('\nset MODELCORDS = ../../data/ModelParams/${SITE}/model_coords_GC_${SITE}\n')
output.write('set DEPFILE = ' + zdepthFile + '\n')
output.write('mpiexec ${BPATH}/${BPROG} \ \n')
output.write('         nx=$NX ny=$NY nz=$NZ \ \n')
output.write('         cordfile=$MODELCORDS \ \n')
output.write('         depfile=$DEPFILE \ \n')
output.write('         modeldir=$MODELDIR \ \n')
output.write('         outdir=$SCRATCH_OUTDIR \ \n')
output.write('         fileroot=$FILEROOT \ \n')
output.write('         logdir=$LOGDIR\n')

output.write('\n\cp -r $SCRATCH_OUTDIR/* $OUTDIR\n')

output.write('\necho "Done   " `date`\n')

output.write('\nexit\n')

output.flush()
output.close()
