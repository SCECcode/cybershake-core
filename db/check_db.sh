#!/bin/bash

set -o errexit

cd /home/scec-00/cybershk/db
java -classpath .:mysql-connector-java-5.0.5-bin.jar:commons-cli-1.0.jar CheckDBDataForSite $@
