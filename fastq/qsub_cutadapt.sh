#!/bin/bash
#######################################################################
#
# Self-submitting script for 
#
#     cutadapat SEE: https://github.com/marcelm/cutadapt
# 
# FLAGS: 
#    [-debug|-inline] run on current node
#    [-clean] delete intermediate files
#    [-pe smp #] run with # cores
#
########################################################################
. /etc/profile.d/modules.sh          # enable module loading
. ~/uab_ngs/uab_ngs_functions_v2.sh  # load shared run_cmd() & run_step()
#
# 
# WARNING: !!!!!NOT RE-ENTRANT!!!!!
#
#$ -S /bin/sh
#$ -cwd
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -j y # merge stderr into stdout
# email at:  Begining, End, Abort, Suspend
#$ -m beas  
# ** RUN TIME ** 
## -l h_rt=119:00:00 -l s_rt=120:55:00 # 5 days
#$ -l h_rt=48:00:00 -l s_rt=47:55:00 # 2 days
#
# *** output logs ***
# -e jobs/$JOB_NAME.$JOB_ID.err # don't need because of -j y
#$ -o jobs/$JOB_NAME.$JOB_ID.out


TASK_NAME=cutadapt
CMD_LINE_PARAM_LIST="WORK_DIR RUN_NAME SAMPLE_NAME ADAPTER FWD_FASTQ REV_FASTQ "
DERIVED_VAR_LIST="CMD_LINE HOSTNAME  DONE_ONLY QSUB_PE_OVERRIDE"

QSUB_DRMAA="-l vf=1.9G -l h_vmem=2G"; if [ -z "$NSLOTS" ]; then export NSLOTS=1; fi
#QSUB_DRMAA="-pe smp 4 -l vf=1.9G -l h_vmem=2G"; if [ -z "$NSLOTS" ]; then export NSLOTS=4; fi

# load needed modules 
# we hit the exe directly
source /share/apps/ngs-ccts/QIIME-files-1.7/QIIME-1.7/activate.sh
export CUTADAPT_EXE=/share/apps/ngs-ccts/QIIME-files-1.7/QIIME-1.7/python-2.7.3-release/bin/cutadapt

#====================================================================== 
# MASTER: submit-self on a head-node
#====================================================================== 
if [[ -z "$JOB_ID" || "$1" == "-inline" ]]; then

    # --------------------------
    # parse parameters
    # --------------------------

    # capture the original CMD_LINE
    if [ -z "$CMD_LINE" ]; then  export CMD_LINE="$0 $*"; fi

    # check for debug or no qsub
    QSUB=`which qsub 2>/dev/null`
    if [ -z "$QSUB" ]; then 
	echo "**** no qsub [NSLOTS=$NSLOTS] ****"
    fi

    #HOSTNAME=DEBUG
    # parse params
    while [[ "$1" == -* ]]; do 
	echo "PARSE FLAG: $1"
	if [[ "-debug" == "$1" || "-inline" == "$1" ]]; then
	    export JOB_ID=run_now
	    echo "****  $1 [NSLOTS=$NSLOTS] ****" 
	    shift 1
	    continue
	fi
	# check for .done mode
	if [ "-done" == "$1" ]; then
	    export DONE_ONLY=yes
	    echo "**** .DONE MODE ON **** (skip step based ONLY on .done file) " 
	    shift 1
	    continue
	fi
	# check for -pe  over-ride
	if [ "-pe" == "$1" ]; then
	    export QSUB_DRMAA="$1 $2 $3 -l vf=1.9G -l h_vmem=2G "
	    export NSLOTS=$3
	    echo "**** -PE OVERRIDE=$QSUB_PE_OVERRIDE **** ) " 
	    shift 1
	    shift 1
	    shift 1
	    continue
	fi
	# check for -clean arguement
	if [[ "-clean" == "$1"  ]]; then
	    echo "**** -CLEAN **** " 
	    export CLEAN=true
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
	echo -n "Z: $myvar	:"; eval echo \$$myvar
	if [ -z "$1" ] ; then
	    echo "$myvar	: MISSING"
	    SCRIPT_NAME=`basename $0`
	    echo ""
	    echo "ERROR: $SCRIPT_NAME [-debug][-pe smp #] ${CMD_LINE_PARAM_LIST}"
	    echo ""
	    echo " "
	    echo ""
	    exit 1
	fi
	shift
    done

    # --------------------------
    # job setup
    # --------------------------

    export CUTADAPT_VER=`$CUTADAPT_EXE --version`
    export DERIVED_VAR_LIST="$DERIVED_VAR_LIST CUTADAPT_VER"

    if [ "$ADAPTER" == "Nextera" ]; then
	# from /share/apps/ngs-ccts/Trimmomatic/Trimmomatic-0.32/adapters/NexteraPE-PE.fa
	# we're loooking for short-fragment read-through, so it's the RC of the adaptor at the other end
	export REV_ADAPTER="Nextera_Trans1_rc=CTGTCTCTTATACACATCTGACGCTGCCGACGA"
	export FWD_ADAPTER="Nextera_Trans2_rc=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
    else 
	if [ "$ADAPTER" == "TruSeq3" ]; then
	# from /share/apps/ngs-ccts/Trimmomatic/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa
	# we're loooking for short-fragment read-through, so it's the RC of the adaptor at the other end
	export REV_ADAPTER="TruSeq3_PE1_rc=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
	export FWD_ADAPTER="TruSeq3_PE2_rc=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	else
	    echo "ERROR: unknown adapter: $ADAPTER"
	    qsub_exit 1
	fi
    fi

    # job dirs
    export DIR_LIST=
    export SAMPLE_DIR=${WORK_DIR}/${TASK_NAME}/${RUN_NAME}/${SAMPLE_NAME}	;export DIR_LIST="$DIR_LIST $SAMPLE_DIR"
    export JOB_DIR=${SAMPLE_DIR}/jobs		;export DIR_LIST="$DIR_LIST $JOB_DIR"
    export OLD_JOB_DIR=${SAMPLE_DIR}/old_jobs	;export DIR_LIST="$DIR_LIST $OLD_JOB_DIR"
    echo "DIR_LIST=$DIR_LIST"
    for dir in $DIR_LIST; do
	if [ ! -e "$dir" ]; then 
	    echo mkdir -p "$dir"
	    run_cmd - mkdir -p "$dir"
	fi
    done

    # --------------------------
    # qsub
    # --------------------------
    if [ -z "$JOB_ID" ]; then
	echo -n "${TASK_NAME}:${SAMPLE_NAME}:QSUB:"
	QSUB_NAME="${TASK_NAME}.${SAMPLE_NAME}"
	pushd ${SAMPLE_DIR}
	echo pushd ${SAMPLE_DIR}
	qsub -terse \
	    $QSUB_DRMAA \
	    -N $QSUB_NAME \
	    -M $USER@uab.edu \
	    $0
	popd
	if [ $? != 0 ]; then echo "ERROR: bad return code from QSUB"; exit 1; fi
	exit 0
    else
p	echo "[debug] skipped qsub"
    fi
fi

#====================================================================== 
# actual slave work
#====================================================================== 
if [ -n "$JOB_ID"  ]; then
    echo "-- environment --"
    echo "JOB_NAME: $JOB_NAME"
    echo "JOB_ID: $JOB_ID"
    echo "NSLOTS=$NSLOTS"
    echo "-- cmd line params -- "
    for myvar in $CMD_LINE_PARAM_LIST ; do
	echo -n "$myvar	:"; eval echo \$$myvar
    done
    echo "-- derrived values line params -- "
    for myvar in $DERIVED_VAR_LIST ; do
	echo -n "$myvar	:"; eval echo \$$myvar
    done
    
    echo "I'm a qsub slave: "
    cd ${WORK_DIR}
    # OUT file
    FWD_BASE=${SAMPLE_DIR}/${SAMPLE_NAME}_R1.cutadapt
    REV_BASE=${SAMPLE_DIR}/${SAMPLE_NAME}_R2.cutadapt
    FWD_TMP=${FWD_BASE}.tmp.fastq
    REV_TMP=${REV_BASE}.tmp.fastq
    #FWD_SHORT=${FWD_BASE}.short.fastq #	--too-short-output $FWD_SHORT \
    #REV_SHORT=${REV_BASE}.short.fastq #	--too-short-output $REV_SHORT \
    FWD_OUT=${FWD_BASE}.fastq
    REV_OUT=${REV_BASE}.fastq
    FWD_INFO=${FWD_BASE}.info
    REV_INFO=${REV_BASE}.info
    FWD_STATS=${FWD_BASE}.stats
    REV_STATS=${REV_BASE}.stats
    FWD_GZ=${FWD_OUT}.gz
    REV_GZ=${REV_OUT}.gz

    # overlap 10 (default 3)
    # error rate 0.2 (default lower)
    export CUTADAPT_PARAMS="\
	--overlap=5 \
	--error-rate=0.20 \
	--minimum-length=20 \
	--quality-cutoff=25 \
        "
    echo "CUTADAPT_PARAMS=$CUTADAPT_PARAMS"

    # FWD reads
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    run_step $SAMPLE_NAME $FWD_GZ "cutadapt_FWD_PE" $FWD_STATS \
	$CUTADAPT_EXE \
	$CUTADAPT_PARAMS \
	--paired-output=$REV_TMP \
	--info-file=$FWD_INFO \
	--output=$FWD_TMP \
	--adapter=$FWD_ADAPTER \
	$FWD_FASTQ \
	$REV_FASTQ
    cat $FWD_STATS

    # REV reads
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    run_step $SAMPLE_NAME $REV_GZ "cutadapt_REV_PE" $REV_STATS \
	$CUTADAPT_EXE \
	$CUTADAPT_PARAMS \
	--paired-output=$FWD_OUT \
	--info-file=$REV_INFO \
	--output=$REV_OUT \
	--adapter=$REV_ADAPTER \
	$REV_TMP \
	$FWD_TMP
    cat $REV_STATS

    # clean up temp files
    run_cmd - \
	rm -rf $FWD_TMP $REV_TMP

    # Compress output
    run_step $SAMPLE_NAME $FWD_GZ "gzip_fwd" - \
	gzip $FWD_OUT
    run_step $SAMPLE_NAME $REV_GZ "gzip_rev" - \
	gzip $REV_OUT


    # kick off GATK, if desired

    #
    # cleanup - if requested
    #
    if [ -n "$CLEAN" ]; then 
	run_cmd - \
	    rm -rf $FWD_INFO $REV_INFO $FWD_SHORT $REV_SHORT
    fi

    exit 0
fi
