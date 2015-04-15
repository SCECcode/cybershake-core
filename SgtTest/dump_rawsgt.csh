#!/bin/csh

if ($# < 3) then
	echo "Usage: $0 <prefix> <nproc> <cs_path>"
	exit(-1)
endif

set FILE_PREFX = $1
set NPROC = $2
set CS_PATH = $3

# Each component gets its own nan file
set NANFILE = nanfile_$FILE_PREFX

mkdir DumpLogs

\rm $NANFILE
touch $NANFILE
set np = 0
while ( $np < $NPROC )

set infile = `echo $FILE_PREFX $np | gawk '{printf "%s%.5d.e3d\n",$1,$2;}'`
set dumpfile = `echo $np | gawk '{printf "DumpLogs/dump.%.5d\n",$1;}'`

${CS_PATH}/SgtTest/bin/dump_rawsgt infile=$infile >& $dumpfile
grep "nan" $dumpfile >> $NANFILE

@ np ++
end

if (-z $NANFILE) then
	echo "No nans detected, cleaning up."
	#rm -rf DumpLogs
	exit(0)
else
	echo "Nans detected in SGTs."
	exit(2)
endif
