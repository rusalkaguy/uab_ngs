#!/bin/bash
. ~/uab_ngs/uab_ngs_functions_v1.sh  # load shared run_cmd() & run_step()
#
# qsub_fastx_trimmer
# 
# 
# trim beginning or end of a fastq/fasta
# 
#$ -S /bin/bash
#$ -cwd
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -j y
# email at:  Begining, End, Abort, Suspend
#$ -m beas  
# stats to fit in verari nodes
#$ -l h_rt=48:00:00
#$ -l vf=900M

TASK_NAME=bwa_sam
CMD_LINE_PARAM_LIST=""
DERIVED_VAR_LIST="FASTX_IN FASTX_OUT FASTX_FIRST FASTX_LAST FASTX_GZIP FASTX_FLAGS CMD_LINE"

# load needed modules 
source /etc/profile.d/modules.sh
module load ngs-ccts/fastx-0.0.13

#====================================================================== 
# MASTER: submit-self on a head-node
#====================================================================== 
if [[ -z "$JOB_ID" || "$1" == "-inline" ]]; then
    

    # --------------------------
    # parse parameters
    # --------------------------

    # capture the original CMD_LINE
    if [ -z "$CMD_LINE" ]; then  CMD_LINE="$0 $*"; fi

    # check for debug or no qsub
    QSUB=`which qsub 2>/dev/null`
    if [[ -z "$QSUB" ]]; then 
	export JOB_ID=run_now
	echo "**** NO QSUB FOUND [NSLOTS=$NSLOTS] ****" 
    fi

    # parse params
    export FASTX_FLAGS="-Q33 "  # sanger/illumina encoding
    while [[ "$1" == -* ]]; do 
	#echo "PARSE FLAG: $1"
	if [[ "-debug" == "$1" || "-inline" == "$1" ]]; then
	    if [[ "-debug" == "$1" || "-inline" == "$1" ]]; then
		#export FASTX_FLAGS="$FASTX_FLAGS $1"
		shift 1
	    fi
	    export JOB_ID=run_now
	    echo "**** NO QSUB [NSLOTS=$NSLOTS] ****" 
	    continue
	fi
	## FASTX_TRIMMER native flags
	if [ "-i" == "$1" ]; then
	    export FASTX_IN=$2
	    #export FASTX_FLAGS="$FASTX_FLAGS $1 $2"
	    shift 2
	    continue
	fi
	if [ "-o" == "$1" ]; then
	    export FASTX_OUT=$2
	    #if [[ "$2" == */ ]]; then 
	    #export FASTX_FLAGS="$FASTX_FLAGS $1 $2"
	    #fi
	    shift 2
	    continue
	fi
	if [ "-f" == "$1" ]; then
	    export FASTX_FIRST=$2
	    export FASTX_FLAGS="$FASTX_FLAGS $1 $2"
	    shift 2
	    continue
	fi
	if [ "-l" == "$1" ]; then
	    export FASTX_LAST=$2
	    export FASTX_FLAGS="$FASTX_FLAGS $1 $2"
	    shift 2
	    continue
	fi
	if [ "-z" == "$1" ]; then
	    export FASTX_GZIP=true
	    export FASTX_FLAGS="$FASTX_FLAGS $1"
	    shift 1
	    continue
	fi
	# unknown flag
	echo "ERROR: unknown option: $1"
	exit 1
    done

    #
    # QSUB
    #
    SAMPLE_NAME=`basename $FASTX_IN`
    JOB_NAME="fastx_trimmer-$SAMPLE_NAME"
    if [ -z "$JOB_ID" ]; then
	echo -n "QSUB: $JOB_NAME "
	qsub -terse \
	    -N ${JOB_NAME} \
	    -M $USER@uab.edu \
	    $0
	RC=$?; if [ $? != 0 ]; then echo "ERROR: bad return code from QSUB=$RC"; exit $RC; fi
	exit 0
    else
	if [ $has_qsub==1 ]; then echo "skipped qsub (debug mode)"; fi
    fi
fi

# 
# SLAVE: do the actual work on a compute node
#
if [ -n "$JOB_ID"  ]; then
    echo "-- environment --"
    echo "SAMPLE_NAME=$SAMPLE_NAME"
    echo "JOB_NAME=$JOB_NAME"
    echo "JOB_ID=$JOB_ID"
    echo "NSLOTS=$NSLOTS"
    echo "-- cmd line params -- "
    print_cmd_line_params
    echo "-- derived params -- "
    print_derived_params
    
    echo "I'm a qsub slave: "

    #
    # figure out output name
    #
    IN_BASE_NAME=`basename $FASTX_IN .gz`
    IN_BASE_NAME=`basename $IN_BASE_NAME .fastq`
    if [[ $FASTX_OUT == */ ]]; then
	# if it's a dir, make sure it exists
	mkdir -p $FASTX_OUT
	FASTX_OUT=$FASTX_OUT$IN_BASE_NAME.trim_${FASTX_FIRST}-${FASTX_LAST}.fastq
    fi
    if [ -z $FASTX_OUT ]; then 
	# if not specified
	FASTX_OUT=`dirname $FASTX_IN`/${IN_BASE_NAME}.trim_${FASTX_FIRST}-${FASTX_LAST}.fastq
    fi
    if [ -n "$FASTX_GZIP" ]; then
	FASTX_OUT=${FASTX_OUT}.gz
    fi
    echo "FASTX_OUT=$FASTX_OUT"

    #
    # format flags
    #
    echo "trimming ${FASTX_IN}"
    run_step $IN_BASE_NAME $FASTX_OUT fastx_trimmer - fastx_trimmer -i $FASTX_IN -o $FASTX_OUT $FASTX_FLAGS

    exit 0
fi
