#!/bin/bash

# Process options
FLAGS=""

# Pass along any arguments to UCVM
while getopts 'b:e:' OPTION
do
  if [ $OPTARG != "" ]; then
      FLAGS="${FLAGS} -$OPTION $OPTARG"
  else
      FLAGS="${FLAGS} -$OPTION"
  fi
done
shift $(($OPTIND - 1))


if [ $# -lt 2 ]; then
	printf "Usage: %s: [options] <infile> <outfile>\n" $(basename $0) >&2    
    	exit 1
fi


IN_FILE=$1
OUT_FILE=$2

echo "./dem_query ${FLAGS} < ${IN_FILE} > ${OUT_FILE}" >> run.log
./dem_query ${FLAGS} < ${IN_FILE} > ${OUT_FILE}
if [ $? -ne 0 ]; then
    exit 1
fi

exit 0
