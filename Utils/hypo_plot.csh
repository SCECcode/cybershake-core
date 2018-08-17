#!/bin/csh

mkdir PlotFiles

set SRFDIR = .
set ERF = 244_5

set PSFILE = PlotFiles/${ERF}_hypos.ps

# Get SRF filenames to loop over from listing of all realizations
#
set FILES = `ls $SRFDIR/${ERF}_v3.3.1-r*.srf `

# Get fault length and width from header of 1st SRF file
# #
# # Later, in for loop over SRFs, each shypo and dhypo will be extracted
# #
#

set FLEN = `gawk '{if(NR==3)print $5;}' $SRFDIR/${ERF}_v3.3.1-r000000.srf `
set FWID = `gawk '{if(NR==3)print $6;}' $SRFDIR/${ERF}_v3.3.1-r000000.srf `

echo $ERF $FLEN $FWID

set KMINCH = 10.0

set XZERO = 1.5
set YZERO = 3.2

set XLAB = "along strike (km)"
set YLAB = "W (km)"

set XINCH = `echo $FLEN ${KMINCH} | gawk '{printf "%f\n",$1/$2;}'`
set YINCH = `echo $FWID ${KMINCH} | gawk '{printf "%f\n",$1/$2;}'`

set SCALE = "$XINCH/-$YINCH"

set REGION = "0/$FLEN/0/$FWID"
set ATTRIB = "-JX$SCALE -R$REGION"

gmtset FRAME_WIDTH 0.05 LABEL_FONT_SIZE 11 ANOT_FONT_SIZE 11 PLOT_DEGREE_FORMAT D PAGE_ORIENTATION PORTRAIT TICK_LENGTH 0.03 D_FORMAT %lg

pstext -JX8.5/11.0 -R0/8.5/0/11 -N -G0/0/0 -K -X0.0 -Y0.0 << END > $PSFILE
#1.2 1.0 20 0 1 1 {source}
#1.2 0.8 10 0 0 1 {SLIPS[1]}
END

pstext -JX8.5/11.0 -R0/8.5/0/11 -N -G0/0/0 -O -K -X$XZERO -Y$YZERO << END >> $PSFILE
END

foreach file ( $FILES )
echo $file

# Get shypo and dhypo from header of this SRF, adjust shypo to start from 0.0
# #

gawk -v fl=$FLEN '{if(NR==4)printf "%f %f\n",$4+0.5*fl,$5;}' $file | \
psxy $ATTRIB -Sx0.10 -W3/255/0/0 -G255/0/0 -K -O >> $PSFILE

end

psxy $ATTRIB -Sx0.10 -W1/255/0/0 -G255/0/0 -B5:"$XLAB":/5:"$YLAB":WSen -K -O << END >> $PSFILE
END

psxy $ATTRIB -W5 -G0/0/0 -O << END >> $PSFILE
END

