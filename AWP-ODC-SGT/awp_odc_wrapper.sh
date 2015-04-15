#!/bin/bash

#Blue Waters wrapper for AWP-ODC, so that we can use jobtype=single and return codes will be picked up correctly

if [ $# -lt 2 ]; then
	echo "Usage: $0 <num cores> <IN3D file>"
	exit 1
fi

NUM_CORES=$1
IN3D_IN=$2

echo "aprun -n $NUM_CORES /projects/sciteam/jmz/CyberShake/software/AWP-ODC-SGT/bin/pmcl3d $IN3D_IN"
aprun -n $NUM_CORES /projects/sciteam/jmz/CyberShake/software/AWP-ODC-SGT/bin/pmcl3d $IN3D_IN

RC=$?

if [ $RC -eq 0 ]; then
	echo "AWP-ODC done " `date`
else
	echo "AWP-ODC FAILED WITH CODE $RC! " `date`
	exit $RC
fi

