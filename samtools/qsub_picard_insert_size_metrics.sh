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
. ~/uab_ngs/uab_ngs_functions_v2.sh  # load shared run_cmd() & run_step()
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
RAM=4 # Gig
CPU=1 # Threads
if [ $CPU -gt 1 ]; then
    FRAM=$(( $RAM / $CPU ))
    export QSUB_MACH_DRMAA="-pe smp ${CPU} -l vf=${FRAM}G -l h_vmem={$FRAM}G" 
else
    export QSUB_MACH_DRMAA="-l vf=${RAM}G -l h_vmem=${RAM}G" # 1 CPUs
fi
RAM_1=$(( $RAM - 1 ))
JAVA_DRMAA="-Xms${RAM_1}g -Xmx${RAM_1}g" 

QSUB_TASK_NAME=PicardCollectInsertSize
# things we'll parse for
CMD_LINE_PARAM_LIST="WORK_DIR SAMPLE_NAME INFILE"
# things we'll report on 
DERIVED_VAR_LIST="CMD_LINE HOSTNAME PROJECT_DIR PICARD_JAR PICARD_VER QSUB_MACH_DRMAA" 

# load needed modules 
module load ngs-ccts/samtools-0.1.19 # to do bam -> bai
module load python/python-2.7
module load R/R-3.0.1
if [ -z "$PICARD_JAR" ]; then PICARD_JAR=/share/apps/ngs-ccts/picard-tools/picard-tools-1.110/CollectInsertSizeMetrics.jar; fi

#====================================================================== 
# MASTER: submit-self on a head-node
#====================================================================== 
if [[ -z "$JOB_ID" || "$1" == "-inline" ]]; then
    
    # --------------------------
    # parse parameters
    # --------------------------
    parse_params $*


    # --------------------------
    # job SETUP
    # --------------------------
    export PICARD_VER=`basename $PICARD_JAR | perl -pe 's|.*/picard-tools-([0-9.]+)$/$1/;'`

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

    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    OUT_BASE=`echo ${INFILE} | perl -pe 's/\.[sb]am$//;'`
    OUT_FILE="${OUT_BASE}.fragstat"
    OUT_HIST="${OUT_BASE}.fraghist"
    run_step $SAMPLE_NAME $OUT_FILE PICARD_CollectInsertSizeMetrics - \
	java $JAVA_DRMAA \
	-Djava.io.tmpdir='/scratch/share/galaxy/temp' \
	-jar $PICARD_JAR \
	VALIDATION_STRINGENCY=LENIENT \
	ASSUME_SORTED=true \
	INPUT=$INFILE \
	OUTPUT=$OUT_FILE \
	HISTOGRAM_FILE=$OUT_HIST 

    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    # https://gist.github.com/davidliwei/2323462
    OUT2_FILE="${OUT_BASE}.insizestat"
    run_step $SAMPLE_NAME $OUT2_FILE DavidLiWei_getinsretsize.py - \
	samtools view $INFILE \
	\| ~/github/davdliwei.gist/getinsertsize.py - \
	\> $OUT2_FILE \

    qsub_exit 0
fi
