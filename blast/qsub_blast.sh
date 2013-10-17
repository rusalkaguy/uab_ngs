#!/bin/bash
#######################################################################
#
# Self-submitting Qsub script.
# 
# qsub_master_slave.sh [-debug|-inline] WORK_DIR SAMPLE_NAME INFILE
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
#QSUB_MACH_DRMAA="-pe smp 4 -l vf=1.9G -l h_vmem=2G" # 4 CPUs, 2G/CPU=8G total
QSUB_MACH_DRMAA="-pe smp 8 -l vf=1.9G -l h_vmem=2G" # 4 CPUs, 2G/CPU=8G total
#QSUB_MACH_DRMAA="-l vf=11.8G -l h_vmem=12G" # 1 CPUs, 12G

QSUB_TASK_NAME=blastn
# things we'll parse for
CMD_LINE_PARAM_LIST="WORK_DIR SAMPLE_NAME INFILE DATABASE_PATH"
# things we'll report on 
DERIVED_VAR_LIST="CMD_LINE HOSTNAME PROJECT_DIR BLAST_VER BLAST_PATH" 

# load needed modules 
module load galaxy/galaxy-command-line
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
    export BLAST_VER=`blastn -version | head -n 1 | perl -pe 's/^[0-9]+([0-9.]+).*$/$1/;'`
    export BLAST_PATH=`which blastn`

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

    DB_NAME=`basename $DATABASE_PATH .fa`
    DATABASE_INDEX_PATH=${DATABASE_PATH}.nsq # Nuc: .nsq, .nin, .nhr #  http://www3.imperial.ac.uk/bioinfsupport/help/blast_searches
    OUT_FILE=${SAMPLE_NAME}.v.$DB_NAME.perc90.ungappedout
    
    # check if we need to do fastq_to_fasta conversion on source
    if [[ $INFILE == *.fastq ]]; then 
	INFILE_FASTQ=$INFILE
	INFILE=`dirname $INFILE`/`basename $INFILE .fastq`.fasta

        # STEP: fastq_to_fasta
        # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
#	run_step $SAMPLE_NAME $INFILE Fastq_to_Fasta - \
#	    fastq_to_fasta -v -i $INFILE_FASTQ -o $INFILE
	run_step $SAMPLE_NAME $INFILE Fastq_2_Fasta $INFILE \
	   ~/uab_ngs/blast/fastq_2_fasta.sh $INFILE_FASTQ 
    fi

    # STEP: Build Blast DB
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    run_step $SAMPLE_NAME $DATABASE_INDEX_PATH blast.formatdb - \
	formatdb -i ${DATABASE_PATH} -p F #-dbtype nucl

    # STEP: Blast
    #          outfmt=6: 12 column standard output
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    run_step $SAMPLE_NAME $OUT_FILE blastn - \
	blastn -db $DATABASE_PATH  \
	-query $INFILE \
	-out $OUT_FILE \
	-task blastn-short \
	-dust no \
	-perc_identity 90 \
	-ungapped \
	-outfmt 6 

    qsub_exit 0
fi
