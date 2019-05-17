#!/bin/bash
# This scripts counts the number of occurrences of function calls.
# Call it by passing in the function names you want to count.
# For example,
#
# numcalls cudaMalloc cudaFree 
#

count_calls () {
        let num=`grep $1\s*\( *.c *.cu -r | wc -l`
        echo $1: $num
}

for call in "$@";
do
        count_calls $call
done
