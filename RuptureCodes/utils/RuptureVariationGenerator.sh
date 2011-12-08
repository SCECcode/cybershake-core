#!/bin/sh

# TODO:  add generation of hypocenters into this script
# TODO:  add insertion of hypocenters into the RV insertion code

TOP_LEVEL_OF_RUPTURES=/home/rcf-104/CyberShake2007/ruptures/RuptureVariations_35_V3.2
RUPTURE_VARIATION_GENERATOR=/home/rcf-104/CyberShake2007/ruptures/ruptureCodes/GenRandV3.2/genslip-v3.2

export LD_LIBRARY_PATH=/home/rcf-104/CyberShake2007/ruptures/ruptureCodes/GenRandV3.2:$LD_LIBRARY_PATH

echo " "
echo "   RuptureVariationsGenerator.sh started. " `date`
echo " "

#
# genslip-mreal outfile=$SRFDIR/${NAME[$m]} < $ERFDIR/$file.erf
#

cd $TOP_LEVEL_OF_RUPTURES
directory_count=0
error_count=0
echo "Finding directories - `date`"
for directory in `find ${TOP_LEVEL_OF_RUPTURES} -type d`; do
  cd ${directory}
  if [ `ls *.txt* | wc -l` -eq 1 ]; then
    echo "Generating variations for ${directory}"
    RUPTURE_FILE=`ls`
    #${RUPTURE_VARIATION_GENERATOR} outfile=${RUPTURE_FILE}.variation < $RUPTURE_FILE > ${RUPTURE_FILE}.variation.output 2>&1
    ${RUPTURE_VARIATION_GENERATOR} infile=$RUPTURE_FILE outfile=${RUPTURE_FILE}.variation &> ${RUPTURE_FILE}.variation.output
    RC=$?
    if [ $RC -ne 0 ]; then
	echo "Error generating ruptures for rupture file $RUPTURE_FILE, exiting."
	exit 2
    fi
  else 
    echo "Multiple files in ${directory}. Skipping"
    error_count=`expr ${error_count} + 1`
  fi   

  # Print out a message for every 100 directories processed 
  directory_count=`expr ${directory_count} + 1`
  if [ `expr ${directory_count} % 100` -eq 0 ]; then
    echo "${directory_count} files copied - `date`"
  fi

done  

echo "${directory_count} directories processed. ${error_count} errors encountered."

echo 
echo "   All done. Goodbye.   " `date` 
echo 
