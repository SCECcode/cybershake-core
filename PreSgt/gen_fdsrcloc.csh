#!/bin/csh

set SITE = LADT

set OUTDIR = ../../data/FdLocs/${SITE}

mkdir -p $OUTDIR

#set SUFX = _sml-h0.2
#set SUFX = -h0.2

set CORDFILE = ../../data/ModelParams/${SITE}/model_coords_GC_${SITE}
set SITEFILE = ../../data/StatInfo/cybersites.ll
set STATS = `gawk '{print $3;}' $SITEFILE`
set SLON  = `gawk '{print $1;}' $SITEFILE`
set SLAT  = `gawk '{print $2;}' $SITEFILE`

@ s = 0
foreach stat ( $STATS )
@ s ++

set SRC_CORDS = $OUTDIR/${stat}.fdloc

gawk -v lon=$SLON[$s] -v lat=$SLAT[$s] 'BEGIN{minr=1.0e+15;}{ \
dx = lon - $1; dy = lat - $2; dr = dx*dx + dy*dy; \
if(dr < minr) { ix=$3; iy=$4; minr=dr;}} \
END{ printf "%d %d\n",ix,iy;}' $CORDFILE > $SRC_CORDS

end
