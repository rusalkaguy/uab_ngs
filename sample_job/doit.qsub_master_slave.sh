#!/bin/bash
#
# Run a set of samples through a script, perhaps filterd by a regex
# 
# doit [-f sample_list.txt] [grep_regex] [script_args]
#
LIST_FILE=doit.qsub_master_slave.sample_list.txt
FILT_REGEX="." 
# parse args
if [ "-i" == "$1" ]; then LIST_FILE=$2; shift 2; fi 
if [ "-e" == "$1" ]; then FILT_REGEX="$2"; shift 2; fi
 
cat $LIST_FILE \
    | grep -v "^#" \
    | grep $FILT_REGEX \
    | ~/uab_ngs/utils/xargs_array `dirname $0`/qsub_master_slave.sh $*
