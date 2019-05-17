#!/usr/bin/bash
# This script adds a call to cudaMemset after each cudaMalloc.
#
# For example, 
#
#       cudaMalloc((void**)&d_sgt_sta, num_bytes);
#
# will be changed to
#
#       cudaMalloc((void**)&d_sgt_sta, num_bytes);
#       cudaMemset(d_sgt_sta, 0, num_bytes);
#

add_memset () {
        sed -i "s/^\(\s*\)\(${1}(\s*(void\*\*)\s*\&\s*\(\w*\)\s*,\s*\(\w*\));\)/\1\2\n\1cudaMemset(\3, 0, \4);/" ${2}
}

for file in *.c *.cu;
do
        sed -i '/^\s*cudaMemset(.*);/d' $file
        add_memset cudaMalloc $file
        add_memset cudaMallocHost $file
done
