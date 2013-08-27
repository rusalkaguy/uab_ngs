#!/bin/bash

#
# Check integrety of Complete Genomics Hard Drive
#

MANIFEST_LIST=`/usr/bin/find . -name manifest.all`
STATUS=manifest.status.txt
echo "sample	status	ok	failed" > $STATUS

for m in ${MANIFEST_LIST}; do
	D=`dirname $m`
	echo "--------------------------- $m ----------------------------"
	pushd $D
	sha256sum -c manifest.all | tee manifest.all.out
	RC=$?	
	COUNT_FAIL=`grep -c " FAILED" manifest.all.out`
	COUNT_OK=`grep -c " OK" manifest.all.out`
	if [ $COUNT_FAIL != 0 ]; then 
		SAMPLE_STATUS="FAILED"
	else
		SAMPLE_STATUS="OK"
	fi
	echo "$D	$SAMPLE_STATUS	$COUNT_OK	$COUNT_FAIL" >> ../../$STATUS
	popd 
done

	
echo "=========================== DONE ==============================="

