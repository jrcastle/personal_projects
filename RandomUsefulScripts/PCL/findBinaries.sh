#!/bin/sh
RUNNUM="$1"
NUMFILES="$2"
COUNTER=0
while [ $COUNTER -lt $NUMFILES ]; do
	find /afs/cern.ch/work/j/jcastle/CMSSW_7_5_5_patch1/src/PbPb3.8T_2015_Alignment/Results$RUNNUM/MinBias_2015_$COUNTER/ -name 'milleBinary*' ! -empty -type f -exec mv -vv {} ./ \;
	let COUNTER=COUNTER+1
done

ls -1 milleBinary_* > BinaryFiles.txt
