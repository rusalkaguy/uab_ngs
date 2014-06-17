#!/bin/bash
#######################################################################
#
# Merge & Trim overlapping reads
#
# Designed for datasets where fragment_size < read_length for many reads
# 
# qsub_master_slave.sh [-debug|-inline] WORK_DIR SAMPLE_NAME INFILE
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
QSUB_MACH_DRMAA="-pe smp 4 -l vf=1.9G -l h_vmem=2G" # 4 CPUs, 2G/CPU=8G total
#QSUB_MACH_DRMAA="-l vf=11.8G -l h_vmem=12G" # 1 CPUs, 12G

QSUB_TASK_NAME=fastq_merge_trim
# things we'll parse for
CMD_LINE_PARAM_LIST="WORK_DIR SAMPLE_NAME FWD_FQ REV_FQ"
# things we'll report on 
DERIVED_VAR_LIST="CMD_LINE HOSTNAME PROJECT_DIR" 

# load needed modules 
module load ngs-ccts/usearch/7.0.1001
#module load ngs-ccts/bwa.0.6.2


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
    #export SAMTOOLS_VER=`samtools 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`

    # job dirs - enumerate and create
    # enumerate
    export DIR_LIST=WORK_DIR
    export DIR_LIST="$DIR_LIST JOB_DIR"; export JOB_DIR=${WORK_DIR}/jobs
    export DIR_LIST="$DIR_LIST SAMPLE_DIR"; export SAMPLE_DIR="$WORK_DIR/$SAMPLE_NAME"
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
    if [ -z "$NSLOTS" ]; then NSLOTS=1; fi
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

    #
    # STEP: TRIMMOMATIC raw fastqs
    #
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    FWD_FQ_TRIM=$SAMPLE_DIR/$SAMPLE_NAME.R1.trim.paired.fastq
    REV_FQ_TRIM=$SAMPLE_DIR/$SAMPLE_NAME.R2.trim.paired.fastq
    FWD_FQ_TRIM_UNPAIR=$SAMPLE_DIR/$SAMPLE_NAME.R1.trim.unpair.fastq
    REV_FQ_TRIM_UNPAIR=$SAMPLE_DIR/$SAMPLE_NAME.R2.trim.unpair.fastq
    PE_FQ_TRIMMED_LOG=$SAMPLE_DIR/$SAMPLE_NAME.PE.trim.log
    run_step $SAMPLE_NAME $FWD_FQ_TRIM trimmomatic_PE - \
	java -Xmx6g \
	-jar /share/apps/ngs-ccts/Trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar \
	PE \
	-threads ${NSLOTS} \
	-trimlog $PE_FQ_TRIMMED_LOG \
	$FWD_FQ $REV_FQ\
	$FWD_FQ_TRIM $FWD_FQ_TRIM_UNPAIR \
	$REV_FQ_TRIM $REV_FQ_TRIM_UNPAIR \
	ILLUMINACLIP:/share/apps/ngs-ccts/Trimmomatic/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:3:50:7 \
	MINLEN:20


    FWD_FQ_TRIM_LHIST=$SAMPLE_DIR/$SAMPLE_NAME.R1.trim.paired.fastq.lhist
    run_step $SAMPLE_NAME $FWD_FQ_TRIM_LHIST seq_length_hist.fwd_fq_trimmed $FWD_FQ_TRIM_LHIST \
	~/uab_ngs/fastq/fastq_seq_len_hist.sh $FWD_FQ_TRIM

    REV_FQ_TRIM_LHIST=$SAMPLE_DIR/$SAMPLE_NAME.R2.trim.paired.fastq.lhist
    run_step $SAMPLE_NAME $REV_FQ_TRIM_LHIST seq_length_hist.rev_fq_trimmed $REV_FQ_TRIM_LHIST \
	~/uab_ngs/fastq/fastq_seq_len_hist.sh $REV_FQ_TRIM

    #
    # STEP: MERGE raw fastq (w/o staggers)
    #
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    FASTQ_MERGE_PP=$SAMPLE_DIR/$SAMPLE_NAME.merge.pp.fastq
    FASTQ_MERGE_PP_OUT=$SAMPLE_DIR/$SAMPLE_NAME.merge.pp.usearch.out
    FASTQ_MERGE_PP_ERR=$SAMPLE_DIR/$SAMPLE_NAME.merge.pp.usearch.err
    run_step $SAMPLE_NAME $FASTQ_MERGE_PP_OUT usearch_merge_proppairs - \
	usearch -threads $NSLOTS \
	        -fastq_mergepairs $FWD_FQ \
		-reverse $REV_FQ  \
                -fastq_truncqual 3 \
		-fastqout $FASTQ_MERGE_PP \
	        -log $FASTQ_MERGE_PP_OUT \
	        2> $FASTQ_MERGE_PP_ERR 

    FASTQ_MERGE_PP_LHIST=$SAMPLE_DIR/$SAMPLE_NAME.merge.pp.fastq.lhist
    run_step $SAMPLE_NAME $FASTQ_MERGE_PP_LHIST seq_length_hist_merge_proppairs $FASTQ_MERGE_PP_LHIST \
	~/uab_ngs/fastq/fastq_seq_len_hist.sh $FASTQ_PP_STAGS


    #
    # STEP: MERGE raw fastq  (WITH staggers)
    #
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    FASTQ_MERGE_WITH_STAGS=$SAMPLE_DIR/$SAMPLE_NAME.merge.stags.fastq
    FASTQ_MERGE_WITH_STAGS_OUT=$SAMPLE_DIR/$SAMPLE_NAME.merge.stags.usearch.out
    FASTQ_MERGE_WITH_STAGS_ERR=$SAMPLE_DIR/$SAMPLE_NAME.merge.stags.usearch.err
    run_step $SAMPLE_NAME $FASTQ_MERGE_WITH_STAGS_OUT usearch_merge_with_stags - \
	usearch -threads $NSLOTS \
	        -fastq_mergepairs $FWD_FQ \
		-reverse $REV_FQ  \
		-fastq_allowmergestagger \
                -fastq_truncqual 3 \
		-fastqout $FASTQ_MERGE_WITH_STAGS \
	        -log $FASTQ_MERGE_WITH_STAGS_OUT \
	        2> $FASTQ_MERGE_WITH_STAGS_ERR 

    FASTQ_MERGE_WITH_STAGS_LHIST=$SAMPLE_DIR/$SAMPLE_NAME.merge.stags.fastq.lhist
    run_step $SAMPLE_NAME $FASTQ_MERGE_WITH_STAGS_LHIST seq_length_hist_stags $FASTQ_MERGE_WITH_STAGS_LHIST \
	~/uab_ngs/fastq/fastq_seq_len_hist.sh $FASTQ_MERGE_WITH_STAGS


    #
    # STEP: MERGE trimmed fastq (w/o staggers)
    #
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    FASTQ_TRIM_MERGE_PP=$SAMPLE_DIR/$SAMPLE_NAME.trim.merge.pp.fastq
    FASTQ_TRIM_MERGE_PP_OUT=$SAMPLE_DIR/$SAMPLE_NAME.trim.merge.pp.usearch.out
    FASTQ_TRIM_MERGE_PP_ERR=$SAMPLE_DIR/$SAMPLE_NAME.trim.merge.pp.usearch.err
    run_step $SAMPLE_NAME $FASTQ_TRIM_MERGE_PP_OUT usearch_trim_merge_proppairs - \
	usearch -threads $NSLOTS \
	        -fastq_mergepairs $FWD_FQ_TRIM \
		-reverse $REV_FQ_TRIM  \
                -fastq_truncqual 3 \
		-fastqout $FASTQ_TRIM_MERGE_PP \
	        -log $FASTQ_TRIM_MERGE_PP_OUT \
	        2> $FASTQ_TRIM_MERGE_PP_ERR 

    FASTQ_TRIM_MERGE_PP_LHIST=$SAMPLE_DIR/$SAMPLE_NAME.trim.merge.pp.fastq.lhist
    run_step $SAMPLE_NAME $FASTQ_TRIM_MERGE_PP_LHIST seq_length_hist_trim_merge_proppairs $FASTQ_TRIM_MERGE_PP_LHIST \
	~/uab_ngs/fastq/fastq_seq_len_hist.sh $FASTQ_TRIM_MERGE_PP


    #
    # STEP: MERGE trimmed fastq  (WITH staggers)
    #
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    FASTQ_TRIM_MERGE_WITH_STAGS=$SAMPLE_DIR/$SAMPLE_NAME.trim.merge.stags.fastq
    FASTQ_TRIM_MERGE_WITH_STAGS_OUT=$SAMPLE_DIR/$SAMPLE_NAME.trim.merge.stags.usearch.out
    FASTQ_TRIM_MERGE_WITH_STAGS_ERR=$SAMPLE_DIR/$SAMPLE_NAME.trim.merge.stags.usearch.err
    run_step $SAMPLE_NAME $FASTQ_TRIM_MERGE_WITH_STAGS_OUT usearch_trim_merge_with_stags - \
	usearch -threads $NSLOTS \
	        -fastq_mergepairs $FWD_FQ_TRIM \
		-reverse $REV_FQ_TRIM  \
		-fastq_allowmergestagger \
                -fastq_truncqual 3 \
		-fastqout $FASTQ_MERGE_WITH_STAGS \
	        -log $FASTQ_MERGE_WITH_STAGS_OUT \
	        2> $FASTQ_TRIM_MERGE_WITH_STAGS_ERR 

    FASTQ_TRIM_MERGE_WITH_STAGS_LHIST=$SAMPLE_DIR/$SAMPLE_NAME.trim.merge.stags.fastq.lhist
    run_step $SAMPLE_NAME $FASTQ_TRIM_MERGE_WITH_STAGS_LHIST seq_length_hist_trim_merge_stags $FASTQ_TRIM_MERGE_WITH_STAGS_LHIST \
	~/uab_ngs/fastq/fastq_seq_len_hist.sh $FASTQ_TRIM_MERGE_WITH_STAGS


    qsub_exit 0
fi
