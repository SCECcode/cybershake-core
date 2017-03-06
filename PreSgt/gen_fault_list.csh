#!/bin/csh

if ( $# != 1 ) then
	echo "Usage: `basename $0` SITE"
	echo ""
	echo "SITE: the site id (e.g. LADT)"
	exit 1
endif

set SITE = $1
set ERF_ID = 29
set RUP_PATH = ../../ruptures/RuptureVariations
set OUTPUT = ../../data/SgtInfo/FaultList/$SITE.faultlist

# Remove output if it already exists
rm -f $OUTPUT

java -classpath .:faultlist/mysql-connector-java-3.1.6-bin.jar faultlist/CreateFaultList $SITE $ERF_ID $RUP_PATH $OUTPUT

set RC = $?

chgrp lc_pjm $OUTPUT
chmod g+rw $OUTPUT

exit $RC
