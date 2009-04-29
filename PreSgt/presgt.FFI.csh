#!/bin/csh

# TWEAK THESE VALUES
set SITE = FFI
set LAT = 34.33603
set LON = -118.50862
set MODEL_LAT = 34.64731
set MODEL_LON = -117.95753
set MODEL_ROT = -55.0
set NX = 2200
set NY = 2900
set NZ = 200
# END TWEAKING

set BIN = `dirname $0`

$BIN/presgt.csh $SITE $LON $LAT $NX $NY $NZ $MODEL_LAT $MODEL_LON $MODEL_ROT
exit $?
