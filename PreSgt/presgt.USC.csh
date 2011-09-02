#!/bin/csh

# TWEAK THESE VALUES
set SITE = USC
set LAT = 34.01919
set LON = -118.28631
set MODEL_LAT = 34.01919
set MODEL_LON = -118.28631
set MODEL_ROT = -90.0
set NX = 2500
set NY = 2500
set NZ = 200
# END TWEAKING

set BIN = `dirname $0`

$BIN/presgt.csh $SITE $LON $LAT $NX $NY $NZ $MODEL_LAT $MODEL_LON $MODEL_ROT
exit $?
