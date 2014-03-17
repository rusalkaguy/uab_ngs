#!/bin/bash
#$ -S /bin/sh
#$ -cwd
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -j y # merge stderr into stdout
# email at:  Begining, End, Abort, Suspend
#$ -m beas  
#
# *** DRMAA ****
#$ -l vf=2.9G -l h_vmem=3G  # observed max of 1.87G for chr1
#$ -l h_rt=1:00:00          # observed max of 30 min for chr1

TASK_NAME=cg_evid2bam
QSUB_DRMAA=
CMD_LINE_PARAM_LIST="TARGET"

# load modules
which cgatools 2>&1 > /dev/null; if [ $? != 0 ]; then 
    echo "****** TO FIX ******"
    echo "module load ngs-ccts/completegenomics-cgatools-1.7.1"
    exit 3 2>&1 > /dev/null
fi
which samtools 2>&1 > /dev/null; if [ $? != 0 ]; then 
    echo "****** TO FIX ******"
    echo "module load samtools"
    exit 3 2>&1 > /dev/null
fi

#
# MASTER or SLAVE?
#    
if [ -z "$JOB_ID" ]; then
    # capture the original CMD_LINE
    if [ -z "$CMD_LINE" ]; then  CMD_LINE="$0 $*"; fi

    # parse params
    FLAGS=
    while [[ "$1" == -* ]]; do 
	# check for debug or no qsub
	if [ "-debug" == "$1" ]; then
	    export FLAGS="$FLAGS -debug"
	    export JOB_ID=run_now
	    echo "**** DEBUG MODE ON ****" 
	    shift 1
	    continue
	fi
	# unknown flag
	echo "ERROR: unknown option: $1"
	exit 1
    done
    # cmd-line params
    for myvar in $CMD_LINE_PARAM_LIST ; do
	eval $myvar=$1
	export $myvar
	echo -n "$myvar	:"; eval echo \$$myvar
	if [[ -z "$1" && "$myvar" != "REV_FASTQ" ]] ; then
	    echo "$myvar	: MISSING"
	    SCRIPT_NAME=`basename $0`
	    echo ""
	    echo "ERROR: $SCRIPT_NAME [-debug] $CMD_LINE_PARAM_LIST"
	    echo ""
	    echo "TARGET may be: "
	    echo "   EVIDENCE directory"
	    echo "or "
	    echo "   a CG evidenceDnbs-chrN-XXX.tsv.bz2 file "
	    echo ""
	    exit 1
	fi
	shift
    done

    # setup
    if [ -d $TARGET ]; then 
	#
        # scan a directory for target files
	#
	TARG_FILES=`find $TARGET -name "evidenceDnbs-*-*.tsv.bz2"`
	for targ in $TARG_FILES; do
	    $0 $FLAGS $targ
	done
    else
	#
	# submit a file
	# 
	if [[ "$TARGET" != *evidenceDnbs-*.tsv.bz2 ]]; then 
	    echo "ERROR: $TARGET : not a evidenceDnbs-*.tsv.bz2 file"
	    exit 1 2>&1 > /dev/null
	fi

	if [ -z "$JOB_ID" ]; then
	    TARGET_NAME=`basename $TARGET .tsv.bz2`
	    QSUB_NAME=${TASK_NAME}.${TARGET_NAME}
	    qsub -terse \
		${QSUB_DRMAA} \
		-N ${QSUB_NAME} \
		-M $USER@uab.edu \
		$0
	    if [ $? != 0 ]; then echo "ERROR: bad return code from QSUB"; exit 1 2>&1 > /dev/null; fi
	    exit 0 2>&1 > /dev/null	
	else 
	    echo "[debug] skipped qsub"
	fi
    fi
fi

# 
# SLAVE: do the actual work on a compute node
#
if [ -n "$JOB_ID"  ]; then
    echo "-- environment --"
    echo "JOB_NAME: $JOB_NAME"
    echo "JOB_ID: $JOB_ID"
    echo "NSLOTS=$NSLOTS"
    echo "-- cmd line params -- "
    for myvar in $CMD_LINE_PARAM_LIST ; do
	echo -n "$myvar	:"; eval echo \$$myvar
    done
    echo "-- derived values line params -- "
    for myvar in $DERIVED_VAR_LIST ; do
	echo -n "$myvar	:"; eval echo \$$myvar
    done
    echo "I'm a qsub slave: "

    DIR=`dirname $TARGET`
    if [ $DIR = "." ]; then 
	DIR="../../BAM/EVIDENCE"
    else
	DIR1=`dirname $DIR`
	DIR2=`dirname $DIR1`
	DIR=$DIR2/BAM/EVIDENCE
    fi
    mkdir -p $DIR
    RC=$?; 
    if [ $RC != 0 ]; then
        echo "ERROR: Couldn't create $DIR: RC=$RC"
        exit 1 2>&1 > /dev/null;
    fi

    RESULT=`basename $TARGET .tsv.bz2`
    RESULT_PATH=$DIR/$RESULT

    touch $DIR/.$RESULT.0.start
    cgatools evidence2sam \
	--beta \
	--evidence-dnbs=$TARGET \
	--reference=${CG_REF}/ReferenceFiles/build37.crr | \
	samtools view -uS - | \
	samtools sort - $RESULT_PATH
    RC=$?; 
    if [ $RC != 0 ]; then 
	echo "ERROR: bad return code from evidence2sam or samtools: RC=$RC"
	exit 1 2>&1 > /dev/null; 
    fi
    samtools index ${RESULT_PATH}.bam ${RESULT_PATH}.bai
    RC=$?; 
    if [ $RC != 0 ]; then 
	echo "ERROR: bad return code from evidence2sam or samtools: RC=$RC"
	exit 1 2>&1 > /dev/null; 
    fi
    touch $DIR/.$RESULT.1.done

    exit 0 2>&1 > /dev/null;     
fi

# shouldn't get here
exit 666 2>&1 > /dev/null;     
