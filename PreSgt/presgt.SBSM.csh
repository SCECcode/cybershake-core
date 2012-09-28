#!/bin/csh

# TWEAK THESE VALUES
set SITE = SBSM
set LAT = 34.064986
set LON = -117.29201
set MODEL_LAT = 34.47073
set MODEL_LON = -117.87082
set MODEL_ROT = -55.0
set NX = 2100
set NY = 2950
set NZ = 200 
# END TWEAKING

set BIN = `dirname $0`

$BIN/presgt.csh $SITE $LON $LAT $NX $NY $NZ $MODEL_LAT $MODEL_LON $MODEL_ROT
exit $?
