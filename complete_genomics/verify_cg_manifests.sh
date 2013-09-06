#!/bin/bash

ls -lst  | email -V -r vera.dpo.uab.edu -f ${USER}@${HOSTNAME}.ad.uab.edu -n “$0” -s "$PWD starting manifest verification" ${USER}@uab.edu

#
# Check integrety of Complete Genomics Hard Drive
#

 if [ -z "$1" ]; then
	# find and fork manifests
	STARTDATE=`date`
	echo $STARTDATE
	MANIFEST_LIST=`/usr/bin/find . -name manifest.all`
	for m in ${MANIFEST_LIST}; do
		D=`dirname $m`
		ASM=`basename $D`
		DD=`dirname $D`
		DID=`basename $DD`
		STATUS=../../$DID.$ASM.status.txt
		echo "--------------------------- $m ----------------------------"
		pushd $D
		/home/curtish/uab_ngs/complete_genomics/verify_cg_manifests.sh $PWD/manifest.all $PWD/$STATUS 
		popd 
	done
	wait
	echo "START: $STARTDATE"
	echo "END  : "`date`
else
	if [ ! -e "$2" ]; then 
		SSTARTDATE=`date`
		echo $SSTARTDATE
		# if summary file does not exist, then process
		
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
		echo "SSTART: $SSTARTDATE"
		echo "SEND  : "`date`
	fi
fi
	
echo "=========================== DONE ==============================="

(echo "START: $STARTDATE\nEND  : "`date`; grep -v sample *.status.txt) | \
	email -V -r vera.dpo.uab.edu -f ${USER}@${HOSTNAME}.ad.uab.edu -n “$0” \
		-s "$PWD completed manifest verification" ${USER}@uab.edu


