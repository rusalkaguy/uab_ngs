#!/bin/bash

# Ingenuity account info
IVA_KEY=$1
IVA_USER=$2
IVA_PASSWORD=$3
if [ -z "$IVA_PASSWORD" ]; then
	echo "USAGE: $0 IVA_KEY IVA_USER IVA_PASSWORD"
	exit 1
fi

#Get list of samples. Define variables. Make a directory called status in the current directory. Pass UAB_Bridges_SampleManifest into variable samplefile, then output that
SAMPLE_BASE_DIR=/base_directory
STATUS_DIR=./status
mkdir -p $STATUS_DIR
SAMPLEFILE=sample_list.txt
SAMPLE_LIST=`cat $SAMPLEFILE`
MANIFEST_FILE=./manifest.xls  #Pass the path and handle for each junction file into  JUNCTION_FILE

#Email User to confirm samples sent
SCRIPT_NAME=`basename $0`   #returns the name of the process running, without the file path directories
send_mail () { # args: subject
which email #returns the full path of shell command email
if [ $? == 0 ]; then #check if the return code of the previous is zero (i.e. everything went OK)
        email -V -r vera.dpo.uab.edu \
                -f ${USER}@${HOSTNAME}.ad.uab.edu \
                -n $SCRIPT_NAME \
                -s "$1" \
                ${USER}@uab.edu
# the previous five lines list the arguments for the email command variables USER and HOSTNAME very helpful in this kind of context
else
        which mailsend #returns path of shell command mailsend
        if [ $? == 0 ]; then #check if the return code of mailsend 
                mailsend -smtp vera.dpo.uab.edu \
                        -f ${USER}@${HOSTNAME}.ad.uab.edu \
                        -name "$SCRIPT_NAME" \
                        -sub "$1" \
                        -t ${USER}@uab.edu
# the previous 5 lines list the arguments for the sendmail command. Here, USER and HOSTNAME are once again used but note that the inputs are very different than for email
        else
                echo "ERROR: can't find either email or mailsend!!!"
        fi
fi
}
cat $0 |  send_mail "program backup of $0"
# This line pipes a printed version of the output of upload_bridges.sh (via cat) through send_mail - thus sending a mail with whats going on in our program to us.



#Iterate over samples
for SAMPLE in $SAMPLE_LIST; do #define a variable SAMPLE as each component of SAMPLE_LIST
	echo "===============================" #create a delimiter for human-readability
	echo -n $SAMPLE":" #prints each SAMPLE as a string
	SAMPLE_DIR=`ls -1d $SAMPLE_BASE_DIR/rpk-slb-?/*/*/$SAMPLE` #Pass each file path into the variable SAMPLE_DIR
	if [[ -z "$SAMPLE_DIR" || ! -e $SAMPLE_DIR ]]; then #If either {the length of sample DIR==0 |  OR |  SAMPLE_DIR does NOT exist}, then echo MISSING
		echo "MISSING"
	else
		echo "found $SAMPLE_DIR" #Otherwise, echo that SAMPLE_DIR was found.
		PARENT_DIR=`dirname $SAMPLE_DIR` #Pass just the path name of SAMPLE_DIR into PARENT_DIR
		ASM_NAME=`basename $PARENT_DIR` #Return just PARENT_DIR name and pass it into ASM_NAME
		ASM_DIR=$SAMPLE_DIR/ASM 
# Now lets create the names of the 3 files that we will be uploading
		VAR_FILE=$ASM_DIR/var-$ASM_NAME.tsv.bz2 #Pass the handle for each Var file into VAR_FILE
		JUNCTION_FILE=$ASM_DIR/SV/highConfidenceJunctionsBeta-$ASM_NAME.tsv #Pass the path and handle for each junction file into  JUNCTION_FILE
                COVERAGE_FILE=$ASM_DIR/CNV/depthOfCoverage_100000-$ASM_NAME.tsv #And for the coverage files
#The next 7 lines simply allow the user to check to see whether the file names are being appropriately constructed
		echo "SAMPLE=$SAMPLE" 
		echo -n "PARENT_DIR=$PARENT_DIR "; if  [ -e  "$PARENT_DIR" ]; then echo "(found)"; else echo "MISSING"; fi
		echo -n "ASM_NAME=$ASM_NAME "; if  [ -e  "$ASM_NAME" ]; then echo "(found)"; else echo "MISSING"; fi
		echo -n "ASM_DIR=$ASM_DIR "; if  [ -e  "$ASM_DIR" ]; then echo "(found)"; else echo "MISSING"; fi
		echo -n "VAR_FILE=$VAR_FILE "; if  [ -e  "$VAR_FILE" ]; then echo "(found)"; else echo "MISSING"; fi 
		echo -n "JUNCTION_FILE=$JUNCTION_FILE "; if  [ -e  "$JUNCTION_FILE" ]; then echo "(found)"; else echo "MISSING"; fi
		echo -n "COVERAGE_FILE=$COVERAGE_FILE "; if  [ -e  "$COVERAGE_FILE" ]; then echo "(found)"; else echo "MISSING"; fi
		for FILE_NAME in $VAR_FILE $JUNCTION_FILE $COVERAGE_FILE; do #Create a variable FILE_NAME for each file to be transferred
			#check if files have already been uploaded
			STATUS_FILE=$STATUS_DIR/`basename $FILE_NAME` #Pass the path and name of FILE_NAME into a new file in STATUS_DIR for verification purposes later
			if [ -e $STATUS_FILE ]; then #If status file exits, then do nothing and echo SKIP
				echo "SKIP: Already Uploaded: $FILENAME"
			else 
				#Otherwise, upload the file. The next 5 lines are the arguments for curl. 
				curl --ftp-ssl --pubkey ing-pub-key.pem \
					 -T $FILE_NAME \
					-u ${IVA_USER}:$PASSWORD \
					--insecure \
					ftp://ftps2.ingenuity.com/for_${IVA_USER}_FTP_upload_$IVA_KEY/
				RC=$? #Ask, What is the status of the previous Curl command? What happened?
				if [ $RC==0 ]; then #If this query returns 0, the command went OK, so echo File Transferred
					echo "File Transferred"
					md5sum $FILE_NAME > $STATUS_FILE #Generate a Checksum and pass it into STATUS_FILE
					echo $FILE_NAME | send_mail "Uploaded: "`basename $FILE_NAME` #Send and email notification the file was uploaded and the file name of the transferred file.
				else
					echo "ERROR: Transfer failed" #Otherwise, notify user of failure and send failure message.
					echo $FILE_NAME | send_mail "Failed: "`basename $FILE_NAME`
				fi
			fi
		done
	fi
done
