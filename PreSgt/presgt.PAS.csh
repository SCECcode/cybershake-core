#!/bin/csh

# TWEAK THESE VALUES
set SITE = PAS
set LAT = 34.1484266
set LON = -118.17119 
set MODEL_LAT = 34.61206
set MODEL_LON = -117.89577
set MODEL_ROT = -55.0
set NX = 2200
set NY = 2850
set NZ = 200 
# END TWEAKING

set BIN = `dirname $0`

$BIN/presgt.csh $SITE $LON $LAT $NX $NY $NZ $MODEL_LAT $MODEL_LON $MODEL_ROT
exit $?
