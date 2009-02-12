#!/bin/csh

if ($# < 3) then
	echo "Usage: $0 <prefix> <nproc> <cs_path>"
	exit(-1)
endif

set FILE_PREFX = $1
set NPROC = $2
set CS_PATH = $3
mkdir DumpLogs

\rm nanfile
touch nanfile
set np = 0
while ( $np < $NPROC )

set infile = `echo $FILE_PREFX $np | gawk '{printf "%s%.4d.e3d\n",$1,$2;}'`
set dumpfile = `echo $np | gawk '{printf "DumpLogs/dump.%.4d\n",$1;}'`

${CS_PATH}/software/SgtTest/dump_rawsgt infile=$infile >& $dumpfile
grep "nan" $dumpfile >> nanfile

@ np ++
end

if (-z nanfile) then
	echo "No nans detected, cleaning up."
	rm -rf DumpLogs
	exit(0)
else
	echo "Nans detected in SGTs."
	exit(2)
endif
