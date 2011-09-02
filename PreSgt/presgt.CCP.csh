#!/bin/csh

# TWEAK THESE VALUES
set SITE = CCP
set LAT = 34.054886
set LON = -118.41302
set MODEL_LAT = 34.52762
set MODEL_LON = -117.97144
set MODEL_ROT = -55.0
set NX = 2100
set NY = 2850
set NZ = 200
# END TWEAKING

set BIN = `dirname $0`

$BIN/presgt.csh $SITE $LON $LAT $NX $NY $NZ $MODEL_LAT $MODEL_LON $MODEL_ROT
exit $?
