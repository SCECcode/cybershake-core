#!/bin/csh

if ( $# != 3 ) then
	echo "Usage: $0 SITE LON LAT"
	echo ""
	echo "SITE: the site identifier (e.g. LADT)"
	echo "LON: the longitude of the site"
	echo "LAT: the latitude of the site"
	echo ""
	echo "example: $0 GID -117.832423 34.34234234"
	exit 1
endif

set SITE = $1
set LON = $2
set LAT = $3

set OUTDIR = ../../data/SgtInfo/FdLocs
set CORDFILE = ../../data/ModelParams/${SITE}/model_coords_GC_${SITE}
set SRC_CORDS = $OUTDIR/${SITE}.fdloc

# Make sure the coordinate file exists
if ( ! -r $CORDFILE ) then
	echo "Cannot read $CORDFILE, have model coordinates been generated for $SITE?"
	exit 1
endif

# Create output directory if it doesn't already exist
mkdir -p $OUTDIR

# Remove output file if it already exists
rm -f $SRC_CORDS

# Find the closest x,y coordinate to the site's lat,lon
echo "Calculating nearest coordinate..."
set RESULT = `gawk -v lon=$LON -v lat=$LAT \
	'BEGIN { minr=1.0e+15; } \
	{ \
		dx = lon - $1; \
		dy = lat - $2; \
		dr = dx*dx + dy*dy; \
		if(dr < minr) { ix=$3; iy=$4; minr=dr; }\
	} \
	END { printf "%d %d\n",ix,iy; }' $CORDFILE`

set RC = $?
if ( $RC != 0 ) exit $RC
echo $RESULT

echo "Writing output..."
echo $RESULT > $SRC_CORDS

chgrp lc_pjm $SRC_CORDS
chmod g+w $SRC_CORDS

echo "Done."

exit 0
