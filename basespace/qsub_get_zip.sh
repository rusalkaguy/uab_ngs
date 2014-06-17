#!/bin/bash
#
# command-line download from BaseSpace
# Source : http://joshquick.github.io/2013/11/14/downloading-from-basespace/
# as summarized by rkumar@uab.edu
# 2014-04-24 script by curtish@uab.edu
#
# get USER/PASS
#
BS_FILE=~/.basespace
if [ ! -e "$BS_FILE" ]; then
    echo "ERROR: $BS_FILE doesn't exist"
    echo "Please create and chmod 700."
    echo "Two required keys: user and password, one per line"
    echo "values, seprated by a tab, on the same line."
    echo "Please create and chmod 700."
    exit 1
fi

BS_USER=`grep ^user ~/.basespace | cut -f 2`
if [ -z "$BS_USER" ]; then echo "ERROR: 'user' key not found in $BS_FILE"; exit 1; else echo "USER=$BS_USER"; fi
BS_PASS=`grep ^password ~/.basespace | cut -f 2`
if [ -z "$BS_PASS" ]; then echo "ERROR: 'password' key not found in $BS_FILE"; exit 1; else echo "PASS=$BS_PASS";fi

BS_URL=$1
if [ -z "$BS_URL" ]; then 
    echo "SYNTAX: $0 BASE_SPACE_URL"
    exit 1
fi

#
# login
#
COOKIES_FILE=.basespace_cookies.$$.txt
POST_DATA="UserName=${BS_USER}&Password=${BS_PASS}"
LOGIN_URL="https://icom.illumina.com/login?service=basespace"
LOGIN_URL="https://accounts.illumina.com/?ForceSsl=True"
LOGIN_URL="https://accounts.illumina.com/auth/logon"

echo "Login with $POST_DATA"
wget \
    --verbose \
    --keep-session-cookies \
    --no-check-certificate \
    --save-cookies $COOKIES_FILE \
    --post-data="$POST_DATA" \
    "https://icom.illumina.com/login?service=basespace"
RC=$?
if [ $RC != 0 ]; then 
    echo "LOGIN FAILED: $RC"
    exit $RC
fi
chmod 700 $COOKIES_FILE


#
# pull URL
#
# Grab the link of the fastq files from basespace website 
# - my link was - https://basespace.illumina.com/runs/5750755/download/5484479/Fastq. 
# Now replace the link in following command, 

wget \
    --load-cookies $COOKIES_FILE \
    --no-check-certificate \
    --keep-session-cookies \
    -O basespace.$$.zip \
    "$BS_URL"
RC=$?
if [ $RC != 0 ]; then 
    echo "DOWNLOAD FAILED: $RC"
    exit $RC
fi

# clenaup cookies
rm $COOKIES_FILE

#
# 3. Extract the archive ( use 7zip as sometime unzip doesn,F"(Bt work)
# 
echo 7za x basespace.$$.zip
7za x basespace.$$.zip
RC=$?
if [ $RC != 0 ]; then 
    echo "EXTRACT FAILED: $RC"
    exit $RC
fi

 
