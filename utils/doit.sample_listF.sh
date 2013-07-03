#!/bin/bash
#
# Run vicuna jobs on all known samples
#
#DEBUG=""; if [ "-debug" == "$1" ]; then DEBUG="-debug"; shift 1; fi
LOC=`dirname $0`
PAT="."; if [ -n "$1" ]; then PAT="$1"; shift; fi
 
cat ../sample_list \
    | grep -v "^#" \
    | grep $PAT \
    | ~/ics/xargs/xargs_array "$LOC/../bwa/doit.bwa_dirs.sh $*"
