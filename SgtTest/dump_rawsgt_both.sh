#!/bin/bash

if [ $# -ne 3 ]; then
	echo "Usage: $0 <xprefix> <yprefix> <nproc>"
	exit 1
fi

X_PREF=$1
Y_PREF=$2
NPROC=$3

SCRIPT_DIR=`dirname $0`
$SCRIPT_DIR/dump_rawsgt.py $X_PREF $NPROC
$SCRIPT_DIR/dump_rawsgt.py $Y_PREF $NPROC

