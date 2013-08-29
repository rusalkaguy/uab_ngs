#!/bin/bash

#
# Check integrety of Complete Genomics Hard Drive
#

if [ -z "$1" ]; then
	# find and fork manifests
	MANIFEST_LIST=`/usr/bin/find . -name manifest.all`
	for m in ${MANIFEST_LIST}; do
		D=`dirname $m`
		ASM=`basename $D`
		DD=`dirname $D`
		DID=`basename $DD`
		STATUS=../../$DID.$ASM.status.txt
		echo "--------------------------- $m ----------------------------"
		pushd $D
		/home/curtish/uab_ngs/complete_genomics/verify_cg_manifests.sh $PWD/manifest.all $PWD/$STATUS &
		popd 
	done
	wait
else
	# process manifest
	sha256sum -c $1 > manifest.all.out
	RC=$?	
	COUNT_FAIL=`grep -c " FAILED" manifest.all.out`
	COUNT_OK=`grep -c " OK" manifest.all.out`
	if [ $COUNT_FAIL != 0 ]; then 
		SAMPLE_STATUS="FAILED"
	else
		SAMPLE_STATUS="OK"
	fi
	echo "sample	status	ok	failed" > $2
	echo "$D	$SAMPLE_STATUS	$COUNT_OK	$COUNT_FAIL" >> $2
fi
	
echo "=========================== DONE ==============================="

