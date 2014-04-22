#!/bin/bash
#######################################################################
#
# Self-submitting script for 
#
#     BWA (MEM) PE Illumna reads vs an un-indexed FASTA genome
# 
# FLAGS: 
#    [-debug|-inline] run on current node
#    [-clean] delete intermediate files
#    [-pe smp #] run with # cores
#    [-bwa exe_path] use that BWA executable
#    [-no_vcf] don't run samtools mpileup
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
#$ -l h_rt=119:00:00
#$ -l s_rt=120:55:00
#
# *** output logs ***
#$ -e jobs/$JOB_NAME.$JOB_ID.err
#$ -o jobs/$JOB_NAME.$JOB_ID.out


TASK_NAME=bwa_mem
CMD_LINE_PARAM_LIST="WORK_DIR SAMPLE_NAME FWD_FASTQ REV_FASTQ REF_FASTA"
DERIVED_VAR_LIST="CMD_LINE HOSTNAME PROJECT_DIR REV_FASTQ DONE_ONLY QSUB_PE_OVERRIDE"

QSUB_DRMAA="-pe smp 4 -l vf=1.9G -l h_vmem=2G" 
if [ -z "$NSLOTS" ]; then export NSLOTS=4; fi  # for debugging

# load needed modules 
# we hit the exe directly
module load ngs-ccts/bwa/0.7.7
module load ngs-ccts/samtools/0.1.19
if [ -z "$PICARD_JAR" ]; then PICARD_JAR=/share/apps/ngs-ccts/picard-tools/picard-tools-1.110/CollectInsertSizeMetrics.jar; fi
module load R/R-3.0.1  # for Picard to make PDF

export DERIVED_VAR_LIST="${DERIVED_VAR_LIST} FASTQ_DIR JOB_DIR"
export DERIVED_VAR_LIST="${DERIVED_VAR_LIST} BWA_VER BWA"
export DERIVED_VAR_LIST="${DERIVED_VAR_LIST} SAMTOOLS_VER SAMTOOLS"

#====================================================================== 
# MASTER: submit-self on a head-node
#====================================================================== 
if [[ -z "$JOB_ID" || "$1" == "-inline" ]]; then

    # capture the original CMD_LINE
    if [ -z "$CMD_LINE" ]; then  CMD_LINE="$0 $*"; fi
    
    # --------------------------
    # parse parameters
    # --------------------------

    # capture the original CMD_LINE
    if [ -z "$CMD_LINE" ]; then  CMD_LINE="$0 $*"; fi

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
	# check for -bwa  over-ride
	if [ "-bwa" == "$1" ]; then
	    export BWA="$2"
	    echo "**** -BWA OVERRIDE=$BWA **** " 
	    shift 1
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
	# check for -no_pileup  over-ride
	if [[ "-no_pileup" == "$1" || "-no_vcf" == "$1" ]]; then
	    echo "**** -NO_PILEUP/-NO_VCF **** " 
	    export NO_PILEUP=true
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
	    echo "ERROR: $SCRIPT_NAME [-debug] [-bwa exe_path] [-pe smp #] ${CMD_LINE_PARAM_LIST}"
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

    #
    # BWA
    #
    if [ -z "$BWA" ]; then export BWA=`which bwa`; fi
    export BWA_DIR=`dirname $BWA`
    # get BWA abbreviations 
    export BWA_VER=`$BWA 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`
    export SAMTOOLS_VER=`samtools 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`

    export INDEX_ABBREV=`basename $REF_FASTA .fa`
    echo "INDEX_ABBREV=$INDEX_ABBREV"
    
    # output file
    export OUT_SEP=-
    export OUT_NAME=${SAMPLE_NAME}${OUT_SEP}${INDEX_ABBREV}${OUT_SEP}bwa${BWA_VER}
    export OUT_DIR=. # relative to WORKDIR

    # job dirs
    export DIR_LIST=
    export JOB_DIR=${WORK_DIR}/jobs		;export DIR_LIST="$DIR_LIST $JOB_DIR"
    export OLD_JOB_DIR=${WORK_DIR}/old_jobs	;export DIR_LIST="$DIR_LIST $OLD_JOB_DIR"
    for dir in $DIR_LIST; do
	if [ ! -e "$dir" ]; then 
	    run_cmd - mkdir -p "$dir"
	fi
    done

    # --------------------------
    # qsub
    # --------------------------
    if [ -z "$JOB_ID" ]; then
	echo -n "${TASK_NAME}:${SAMPLE_NAME}:QSUB:"
	QSUB_NAME="${TASK_NAME}.${SAMPLE_NAME}.${INDEX_ABBREV}"
	pushd ${WORK_DIR}
	qsub -terse \
	    $QSUB_DRMAA \
	    -N $QSUB_NAME \
	    -M $USER@uab.edu \
	    -e ${JOB_DIR}/${QSUB_NAME}.$$.err.txt \
	    -o ${JOB_DIR}/${QSUB_NAME}.$$.out.txt \
	    $0
	popd
	if [ $? != 0 ]; then echo "ERROR: bad return code from QSUB"; exit 1; fi
	exit 0
    else
	echo "[debug] skipped qsub"
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

    # BWA Index of Reference, if needed
    export BWA_REF_INDEXED=${REF_FASTA}.bwa${BWA_VER}
    if [ ! -e "${BWA_REF_INDEXED}.bwt" ]; then
	# decide which indexing option
	REF_SIZE=`wc -c ${REF_FASTA} | cut -d " " -f 1`
	export BWA_INDEX_TYPE='is'
	BWA_MAX_IS=$(( 2**30 ))
	if [ $REF_SIZE -gt $BWA_MAX_IS ]; then
	    export BWA_INDEX_TYPE='bwtsw'
	fi
	echo "REF_SIZE_=$REF_SIZE"
	echo "BWA_MAX_IS=$BWA_MAX_IS"
	echo "BWA_INDEX_TYPE=$BWA_INDEX_TYPE"
    fi
    echo "BWA_REF_INDEXED=$BWA_REF_INDEXED"
    run_step $SAMPLE_NAME $BWA_REF_INDEXED.bwt BWA_index_ref_fasta - \
	bwa index -a $BWA_INDEX_TYPE -p $BWA_REF_INDEXED $REF_FASTA    

    # BWA mem
    SAM=${OUT_DIR}/${OUT_NAME}.sam.gz
    UBAM=${OUT_DIR}/${OUT_NAME}.unsorted.bam
    BAM_BASE=${OUT_DIR}/${OUT_NAME}
    BAM=${BAM_BASE}.bam

    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    echo "FWD_FASTQ=$FWD_FASTQ"
    # -R \"$RG_NAME\" \  # readgroup
    run_step $SAMPLE_NAME $SAM "BWA_mem" $SAM \
	$BWA mem -t $NSLOTS \
	-M -C \
	$BWA_REF_INDEXED \
	$FWD_FASTQ $REV_FASTQ \
	\| gzip 


    # .SAM.gz to sorted .BAM
    run_step $SAMPLE_NAME $UBAM "samtools_view" $UBAM \
	zcat $SAM \|  \
	sed "'s/\([1-9]\:N\:[0-9][0-9]*\:[A-Z]*\)/BC\:Z\:\1/'" \| \
	samtools view -bS -@ $NSLOTS -
    run_step $SAMPLE_NAME $BAM "samtools_sort" $BAM \
	  samtools sort -@ $NSLOTS $UBAM $BAM_BASE

    # SAMTools index
    BAI=${OUT_DIR}/${OUT_NAME}.bai
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    run_step $SAMPLE_NAME $BAI "samtools_index" - \
	samtools index ${BAM} ${BAI}

    # SAMTools flagstat
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    FLAGSTAT_OUT=${OUT_DIR}/${OUT_NAME}.flagstat
    run_step $SAMPLE_NAME $FLAGSTAT_OUT "samtools_index" $FLAGSTAT_OUT \
	samtools flagstat ${BAM}

    # PICARD insertSizeMetrics
    JAVA_RAM=$(( ($NSLOTS * 2 ) - 1 ))
    JAVA_DRMAA="-Xms${JAVA_RAM}g -Xmx${JAVA_RAM}g" 
    OUT_FILE="${BAM_BASE}.fragstat"
    OUT_HIST="${BAM_BASE}.fraghist.pdf"
    run_step $SAMPLE_NAME $OUT_FILE PICARD_CollectInsertSizeMetrics - \
	java $JAVA_DRMAA \
	-Djava.io.tmpdir='/scratch/share/galaxy/temp' \
	-jar $PICARD_JAR \
	VALIDATION_STRINGENCY=LENIENT \
	ASSUME_SORTED=true \
	INPUT=$BAM \
	OUTPUT=$OUT_FILE \
	HISTOGRAM_FILE=$OUT_HIST 


    if [ -z "$NO_PILEUP" ]; then
	# SAMTOOLs PILEUP -> VCF
	SAMTOOLS_REF_INDEXED=${REF_FASTA}.fai
	run_step $SAMPLE_NAME $SAMTOOLS_REF_INDEXED SAMTOOLS_index_ref -  \
	    samtools faidx ${REF_FASTA}
	BCF_RAW=${BAM}.raw.bcf
	run_step $SAMPLE_NAME $BCF_RAW SAMTOOLS_mpileup_vcf $BCF_RAW  \
	    samtools mpileup -C50 -g -f ${REF_FASTA} $BAM 
	VCF_OUT=${BAM}.vcf
	run_step $SAMPLE_NAME $VCF_OUT SAMTOOLS_mpileup_vcf $VCF_OUT  \
	    bcftools view -cv $BCF_RAW
    fi

    # kick off GATK, if desired

    #
    # cleanup - if requested
    #
    if [ -n "$CLEAN" ]; then 
	run_cmd - \
	    rm -rf $SAM $UBAM $BCF_RAW 
    fi

    exit 0
fi
