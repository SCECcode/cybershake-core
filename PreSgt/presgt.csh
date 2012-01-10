#!/bin/csh

if ( $# != 9 ) then
	echo "Usage: `basename $0` SITE LON LAT NX NY NZ MODEL_LAT MODEL_LON MODEL_ROT"
	exit 1
endif

set SITE = $1
set LON = $2
set LAT = $3
set NX = $4
set NY = $5
set NZ = $6
set MODEL_LAT = $7
set MODEL_LON = $8
set MODEL_ROT = $9
set BIN = `dirname $0`

echo "Calculating source location..."
$BIN/gen_fdsrcloc.csh $SITE $LON $LAT
if ( $? != 0 ) then
	echo "ERROR: Failed to calculate site coordinates!"
	exit 1
endif

echo "Generating fault list..."
$BIN/gen_fault_list.csh $SITE
if ( $? != 0 ) then
	echo "ERROR: Failed to generate fault list!"
	exit 1
endif

echo "Generating sgt grid..."
$BIN/gen_sgtgrid.csh $SITE $NX $NY $NZ $MODEL_LAT $MODEL_LON $MODEL_ROT
if ( $? != 0 ) then
	echo "ERROR: Failed to generate SGT grid!"
	exit 1
endif

echo "Completed $SITE successfully."
