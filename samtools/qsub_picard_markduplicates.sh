#!/bin/bash
#######################################################################
#
# Self-submitting Qsub script.
# 
# qsub_master_slave.sh [-debug|-inline] WORK_DIR SAMPLE_NAME INFILE
#
# http://seqanswers.com/forums/showthread.php?t=9194
#
########################################################################
# libraries to load
. /etc/profile.d/modules.sh          # enable module loading
. ~/uab_ngs/uab_ngs_functions_v1.sh  # load shared run_cmd() & run_step()
#*** QSUB FLAGS ***
#$ -S /bin/bash
#$ -cwd # remember what dir I'm launched in 
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -m beas  #email at:  Begining, End, Abort, Suspend
# *** output logs ***
#$ -j y # merge stderr into stdout
#$ -o jobs/$JOB_NAME.$JOB_ID.out
#*** END QSUB ****
# run time
QSUB_TIME_DRMAA="-l h_rt=119:00:00 -l s_rt=120:55:00 " # 
RAM=12 # Gig
CPU=1 # Threads
if [ $CPU -gt 1 ]; then
    FRAM=$(( $RAM / $CPU ))
    export QSUB_MACH_DRMAA="-pe smp ${CPU} -l vf=${FRAM}G -l h_vmem={$FRAM}G" 
else
    export QSUB_MACH_DRMAA="-l vf=${RAM}G -l h_vmem=${RAM}G" # 1 CPUs
fi
RAM_1=$(( $RAM - 1 ))
JAVA_DRMAA="-Xms${RAM_1}g -Xmx${RAM_1}g" 

QSUB_TASK_NAME=PicMarkDup
# things we'll parse for
CMD_LINE_PARAM_LIST="WORK_DIR SAMPLE_NAME INFILE"
# things we'll report on 
DERIVED_VAR_LIST="CMD_LINE HOSTNAME PROJECT_DIR PIC_JAR SAMTOOLS_VER QSUB_MACH_DRMAA" 

# load needed modules 
module load ngs-ccts/samtools-0.1.19 # to do bam -> bai
PIC_JAR=/share/apps/galaxy/galaxy-tool-dependencies/picard/1.56.0/devteam/picard/e0232cbac965/jars/MarkDuplicates.jar

#====================================================================== 
# MASTER: submit-self on a head-node
#====================================================================== 
if [[ -z "$JOB_ID" || "$1" == "-inline" ]]; then
    
    # --------------------------
    # parse parameters
    # --------------------------

    # capture the original CMD_LINE
    if [ -z "$CMD_LINE" ]; then  CMD_LINE="$0 $*"; fi

    # check for debug or no qsub - so it will work outside the cluster
    QSUB=`which qsub 2>/dev/null`
    if [[ -z "$QSUB" ]]; then 
	export JOB_ID=run_now
	echo "**** NO QSUB FOUND ****" 
    fi

    # parse FLAGS on cmd-line
    while [[ "$1" == -* ]]; do 
	echo "PARSE FLAG: $1"
	if [[ "-debug" == "$1" || "-inline" == "$1" ]]; then
	    shift 1
	    export JOB_ID=run_now
	    echo "**** NO QSUB [NSLOTS=$NSLOTS] ****" 
	    continue
	fi
	# unknown flag
	echo "ERROR: unknown option: $1"
	qsub_exit 1
    done
    # parse PARAMS on cmd-line 
    for myvar in $CMD_LINE_PARAM_LIST ; do
	eval $myvar=$1
	export $myvar
	echo -n "Z: $myvar	:"; eval echo \$$myvar
	if [ -z "$1" ] ; then
	    echo "$myvar	: MISSING"
	    SCRIPT_NAME=`basename $0`
	    echo ""
	    echo "ERROR: $SCRIPT_NAME [-debug|-inline] ${CMD_LINE_PARAM_LIST}"
	    echo ""
	    echo ""
	    qsub_exit 1
	fi
	shift
    done

    # --------------------------
    # job SETUP
    # --------------------------
    export SAMTOOLS_VER=`samtools 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`

    # job dirs - enumerate and create
    # enumerate
    export DIR_LIST=WORK_DIR
    export DIR_LIST="$DIR_LIST JOB_DIR"; export JOB_DIR=${WORK_DIR}/jobs
    # create
    for dir in $DIR_LIST; do
	eval MDIR=\$$dir
	if [ ! -e ${MDIR} ]; then 
	    run_cmd - mkdir -p $MDIR
	fi
    done

    # --------------------------
    # QSUB myself to do the work on a node
    # --------------------------
    if [ -z "$JOB_ID" ]; then
	echo -n "${QSUB_TASK_NAME}:${SAMPLE_NAME}:QSUB:"
	QSUB_NAME=${QSUB_TASK_NAME}-${SAMPLE_NAME}
	if [[ $0 == /* ]]; then JOB_SCRIPT=$0; else JOB_SCRIPT=`pwd`/$0; fi
	pushd $WORK_DIR > /dev/null
	qsub -terse \
	    $QSUB_TIME_DRMAA $QSUB_MACH_DRMAA \
	    -N $QSUB_NAME \
	    -M $USER@uab.edu \
	    $JOB_SCRIPT
	popd > /dev/null
	if [ $? != 0 ]; then echo "ERROR: bad return code from QSUB"; qsub_exit 1; fi
    else
	echo "[debug] skipped qsub"
	$0
	RC=$?
	if [ $? != 0 ]; then echo "ERROR: bad return code from slave $0: $RC"; qsub_exit $RC; fi
    fi
    qsub_exit 0
fi

#====================================================================== 
# actual slave work
#====================================================================== 
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
    cd $WORK_DIR

    # STEP: write hello file
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    OUT_BAM=${SAMPLE_NAME}.mkdup.rm.bam
    OUT_BAI=${SAMPLE_NAME}.mkdup.rm.bai
    OUT_STATS=${SAMPLE_NAME}.mkdup.rm.txt
    OUT_FLAGSTATS=${SAMPLE_NAME}.mkdup.rm.flagstat
    run_step $SAMPLE_NAME $OUT_BAM PICARD_mark_dup - \
	java $JAVA_DRMAA \
	-Djava.io.tmpdir='/scratch/share/galaxy/temp' \
	-jar $PIC_JAR \
	VALIDATION_STRINGENCY=LENIENT \
	ASSUME_SORTED=true \
	INPUT=$INFILE \
	OUTPUT=$OUT_BAM \
	METRICS_FILE=$OUT_STATS \
	REMOVE_DUPLICATES=true \
	READ_NAME_REGEX="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100

    run_step $SAMPLE_NAME $OUT_BAI SAMTOOLS_index - \
	samtools index \
	${OUT_BAM} \
	${OUT_BAI}

    run_step $SAMPLE_NAME $OUT_FLAGSTATS SAMTOOLS_flagstat $OUT_FLAGSTATS \
	samtools flagstat \
	${OUT_BAM} \

    qsub_exit 0
fi
