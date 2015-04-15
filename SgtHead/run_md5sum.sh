#!/bin/bash

if [ $# -lt 2 ]; then
	echo "Usage: $0 <file>"
fi

aprun -n 1 md5sum $1 > $1.md5
