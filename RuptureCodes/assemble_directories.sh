#!/bin/sh

# This script takes a directory full of rupture variations, creates a bunch 
# of subdirectories based on the rupture source number and moves those 
# rupture files into the subdirectories.
#
# Rupture files need to have the following naming convention.
#  <S#>_<R#>.txt
#  For example,
#  1_5.txt
#  3_1.txt
#  3_2.txt
#  100_1.txt
#
#  Created directories will be of the form <S#>/<R#>


if [ $# -ne 2 ]; then
  echo "Usage: $0 <source directory> <destination directory>"
  exit
fi

SOURCE_DIR=$1
DEST_DIR=$2

count=0
echo "Getting file list - `date`"
for file in `ls ${SOURCE_DIR}`; do

  if [ ${count} -eq 0 ]; then
    echo " Started copying files - `date`"
  fi

  S=`echo ${file} | cut -d _ -f 1`
  R=`echo ${file} | cut -d _ -f 2 | cut -d . -f 1`
 
  mkdir -p ${DEST_DIR}/${S}/${R} 
  cp ${SOURCE_DIR}/${file} ${DEST_DIR}/${S}/${R}

  count=`expr ${count} + 1`
  if [ `expr ${count} % 1000` -eq 0 ]; then
    echo "${count} files copied - `date`"
  fi

done
