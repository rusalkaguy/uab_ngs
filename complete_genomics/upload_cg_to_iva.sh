#!/bin/bash
######################################################################
# 
# Upload 3 CG files for a list of samples
#
# 3 files = {var,highConfidenceJunctionsBeta,depthOfCoverage_100000}
#
# Requires:
# * Ingenuity account info {KEY, ACCOUNT, PASSWORD}
# * SMTP server
# * notification email address
# * location of Complete Genomic disk images
#
# SECURITY WARNING: password will appear on curl commandline
#    do not use this script on a non-secure or multi-user server
# 
# Created: 2013-11-08 Vincent Laufer & Curtis Hendrickson
# 
######################################################################
SCRIPT_NAME=`basename $0`   #returns the name of the process running, without the file path directories

# PARAMETERS
IVA_KEY=$1       # hash value
IVA_USER=$2      # usually email address
IVA_PASSWORD=$3  
SMTP_SERVER=$4   # remote/local mail server
CG_DISK_DIR=$5   # path to disk images of CG data.
EMAIL_CC=$6      # carbon copy status emails
# USAGE MESSAGE
if [ -z "$IVA_PASSWORD" ]; then
	echo "USAGE: $SCRIPT_NAME IVA_KEY IVA_USER IVA_PASSWORD SMTP_SERVER CG_DISK_DIR EMAIL_CC"
	exit 1
fi

#Get list of samples. Define variables. Make a directory called status in the current directory. Pass UAB_Bridges_SampleManifest into variable samplefile, then output that
SAMPLE_BASE_DIR=$CG_DISK_DIR
STATUS_DIR=./status
mkdir -p $STATUS_DIR
if [ $? != 0 ]; then echo "ERROR: couldn't create $STATUS_DIR"; exit 1; fi
SAMPLEFILE=sample_list.txt  # on id per line
if [ ! -e $SAMPLE_FILE ]; then echo "ERROR: couldn't find sample_list.txt to read sample IDs"; exit 1; fi
SAMPLE_LIST=`cat $SAMPLEFILE`
MANIFEST_FILE=./manifest.xls  # manifest file to upload giving sample metadata

#Email User to confirm samples sent
send_mail () { # args: subject
    which email #returns the full path of shell command email
    if [ $? == 0 ]; then #check if the return code of the previous is zero (i.e. everything went OK)
        email -V -r $SMTP_SERVER \
	    -f ${USER}@`domainname -f` \
	    -n $SCRIPT_NAME \
	    -s "$1" \
	    $IVA_USER $EMAIL_CC
    else
        which mailsend #returns path of shell command mailsend
        if [ $? == 0 ]; then #check if the return code of mailsend 
                # use mailsend to send email
                mailsend -smtp $SMTP_SERVER \
		    -f ${USER}@$`domainname -f` \
		    -name "$SCRIPT_NAME" \
		    -sub "$1" \
		    -t $IVA_USER \
		    -cc $EMAIL_CC
          
        else
                echo "ERROR: can't find either email or mailsend!!!"
        fi
    fi
}


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
			if [ ! -e $FILE_NAME ]; then # if file doesn't exist, complain
			    echo "ERROR: file not found: $FILE_NAME"
			else
			    if [ -e $STATUS_FILE ]; then #If status file exits, then do nothing and echo SKIP
				echo "SKIP: Already Uploaded: $FILE_NAME"
			    else 
				#Otherwise, upload the file. The next 5 lines are the arguments for curl. 
				curl --ftp-ssl --pubkey ing-pub-key.pem \
				    -T $FILE_NAME \
				    -u ${IVA_USER}:$PASSWORD \
				    --insecure \
				    ftp://ftps2.ingenuity.com/for_${IVA_USER}_FTP_upload_$IVA_KEY/ \
				    2>&1 | tee ${STATUS_FILE}.out
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
			fi
		done
	fi
done
