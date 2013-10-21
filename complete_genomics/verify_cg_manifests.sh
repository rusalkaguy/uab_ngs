#!/bin/bash
#
# Check integrety of Complete Genomics Hard Drive
#
#
# USAGES
# 1. verify_cg_manifests.sh # find MANIFESTs in PWD
# 2. verify_cg_manifests.sh drive_path # cd to $drive_path, find MANIFESTs
# 3. verify_cg_manifests.sh manifest.all_path status_out_path # verify MANIFEST and summarize to out

SCRIPT_NAME=`basename $0`

send_mail () { # args: subject
which email
if [ $? == 0 ]; then
	email -V -r vera.dpo.uab.edu \
		-f ${USER}@${HOSTNAME}.ad.uab.edu \
		-n $SCRIPT_NAME \
		-s "$1" \
		${USER}@uab.edu 
else
	which mailsend
	if [ $? == 0 ]; then
		mailsend -smtp vera.dpo.uab.edu \
			-f ${USER}@${HOSTNAME}.ad.uab.edu \
			-name "$SCRIPT_NAME" \
			-sub "$1" \
			-t ${USER}@uab.edu
	else
		echo "ERROR: can't find either email or mailsend!!!"
	fi
fi	
}


CHECKSUM_OUT="manifest.all.$HOSTNAME.out"
if [ -z "$2" ]; then
	# cd to target directory, if one is given
	if [ -n "$1" ]; then 
		cd $1
	fi
	# email start
	ls -lst | send_mail "$PWD starting drive manifest verification" 

	# find and process manifests
	STARTDATE=`date`
	echo $STARTDATE
	DRIVE=`basename $PWD`
	STATUS_DIR=sha256sum.$HOSTNAME
	mkdir -p $STATUS_DIR
	echo "mkdir -p $STATUS_DIR RC=$?"
	MANIFEST_LIST=`/usr/bin/find . -name manifest.all`
	for m in ${MANIFEST_LIST}; do
		export D=`dirname $m`; echo "D=$D"
		ASM=`basename $D`; echo "ASM=$ASM"
		DD=`dirname $D`; echo "DD=$DD"
		DID=`basename $DD`; echo "DID=$DID"
		STATUS=../../$STATUS_DIR/$DID.txt
		echo "--------------------------- $m ----------------------------"
		pushd $D
		~/uab_ngs/complete_genomics/verify_cg_manifests.sh $PWD/manifest.all $PWD/$STATUS 
		popd 
	done
	wait
	echo "START: $STARTDATE"
	echo "END  : "`date`
	echo "=========================== DONE ==============================="
	(echo "START: $STARTDATE"; echo "END  : "`date`; ls -lst $STATUS_DIR; sort -r $STATUS_DIR/*.txt | uniq | perl -pe '$_="> $_";') | \
		send_mail "$PWD completed drive manifest verification"

else
	if [ ! -e "$2" ]; then 
		# email start
		ls -lst | send_mail "$PWD starting sample manifest verification" 
		SSTARTDATE=`date`
		echo $SSTARTDATE
		# if summary file does not exist, then process
		
		# process manifest
		echo "[$PWD] sha256sum -c $1 > $CHECKSUM_OUT"
		#sha256sum -c $1 > $CHECKSUM_OUT
		RC=$?	
		COUNT_FAIL=`grep -c " FAILED" $CHECKSUM_OUT`
		COUNT_OK=`grep -c " OK" $CHECKSUM_OUT`
		if [ $COUNT_FAIL != 0 ]; then 
			SAMPLE_STATUS="FAILED"
		else
			SAMPLE_STATUS="OK"
		fi
		echo "sample	status	ok	failed" > $2
		echo "$D	$SAMPLE_STATUS	$COUNT_OK	$COUNT_FAIL" >> $2
		echo "SSTART: $SSTARTDATE"
		echo "SEND  : "`date`
		(echo "START: $SSTARTDATE"; echo "END  : "`date`; sort -r $2 | perl -pe '$_="> $_";' ) | \
			send_mail "$PWD completed sample manifest verification" 
	else
		(echo "START: $SSTARTDATE"; echo "END  : "`date`; sort -r $2| perl -pe '$_="> $_";') | \
			send_mail "$PWD SKIPPED sample manifest verification" 
	fi
	echo "=========================== DONE ==============================="

fi
	


