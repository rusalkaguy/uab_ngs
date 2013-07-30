#!/bin/bash
#
# Take read-group information & fastqs as inputs
# Self-submit
#
# designed for large sample counts with multiple read sets per sample
#
#$ -S /bin/sh
#$ -cwd
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -j y # merge stderr into stdout
# email at:  Begining, End, Abort, Suspend
#$ -m beas  
#
# *** DRMAA resources for BWA ALN ****
# worked for 85% of kimbley samples
# -pe smp 4
# -l vf=1.9G -l h_vmem=2G
#
# *** DRMAA resources for BWA SAMPE ONLY ****
# -l vf=3.9G -l h_vmem=4G
#
# ** RUN TIME ** 
#$ -l h_rt=119:00:00
#$ -l s_rt=120:55:00
#
# *** output logs ***
#$ -e jobs/$JOB_NAME.$JOB_ID.err
#$ -o jobs/$JOB_NAME.$JOB_ID.out

TASK_NAME=bwa_sam
CMD_LINE_PARAM_LIST="SAMPLE_NAME RG_ID RG_LIB FWD_FASTQ"
DERIVED_VAR_LIST="CMD_LINE HOSTNAME PROJECT_DIR REV_FASTQ DONE_ONLY BWA BWA_VER SAMTOOLS_VER INDEX INDEX_ABBREV OUT_DIR OUT_NAME RG_NAME QSUB_PE_OVERRIDE"

export PROJECT_DIR=/scratch/user/${USER}/kimberly_lupus

QSUB_DRMAA="-pe smp 4 -l vf=1.9G -l h_vmem=2G" # worked for all but 11 kimberly samples - those crash
#QSUB_DRMAA="-l vf=3.9G -l h_vmem=4G" # DRMAA resources for BWA SAMPE ONLY 
# change default because parameter based over-ride didn't work and out of time to debug
#QSUB_DRMAA="-pe smp 8 -l vf=1.9G -l h_vmem=2G" # try bigger for final 11 kimberly


#
# MASTER: submit-self on a head-node
#
if [ -z "$JOB_ID" ]; then
    # capture the original CMD_LINE
    if [ -z "$CMD_LINE" ]; then  CMD_LINE="$0 $*"; fi

    # parse params
    while [[ "$1" == -* ]]; do 
	# check for debug or no qsub
	if [ "-debug" == "$1" ]; then
	    export JOB_ID=run_now
	    echo "**** DEBUG MODE ON ****" 
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
	# check for -bwa  over-ride
	if [ "-bwa" == "$1" ]; then
	    export BWA="$2"
	    echo "**** -BWA OVERRIDE=$BWA **** " 
	    shift 1
	    shift 1
	    continue
	fi
	# check for -index  over-ride
	if [ "-index" == "$1" ]; then
	    export INDEX="$2"
	    echo "**** -INDEX OVERRIDE=$INDEX **** " 
	    shift 1
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
	    echo "ERROR: $SCRIPT_NAME [-debug] RG_sample_name RG_id RG_lib FWD.fastq[i.gz] REV.fastq[i.gz]"
	    echo ""
	    echo "Read Group fields"
	    echo "  RG_id = flow_cell.LANE#"
	    echo "  RG_lib = library/prep/batch name"
	    echo "  RG_sample_name = sample name"
	    echo "If there are multiple preps of a sample, "
	    echo "they should have different RG_id's or RG_libs'"
	    echo ""
	    echo "FWD/REV read FASTQ files may (optionally) be compressed."
	    echo "If .fastqi suffix is used, then -I (Illumina 1.3 Phred scores)"
	    echo "  will be passed to BWA"
	    echo ""
	    echo "Output files will be named: "
	    echo " RG_sample_name-RG_lib-RG_id.[whatever]"
	    echo "so avoid '-'s in your names, stick to underscores and periods!"
	    echo " "
	    echo ""
	    exit 1
	fi
	shift
    done

    # compute REV_FASTQ filename
    if [ -z "$REV_FASTQ" ]; then
	export REV_FASTQ=`echo $FWD_FASTQ | perl -pe 's/R1\./R3./;s/_1\.ILM/_3.ILM/;'`
    fi

    # BWA and INDEX versions
    if [ -z "$BWA" ]; then export BWA="/share/apps/ngs-ccts/bwa-0.6.2/bwa"; fi
    #if [ -z "$BWA" ]; then export BWA="/share/apps/ngs-ccts/bwa-0.5.9/bwa"; fi
    export BWA_DIR=`dirname $BWA`
    if [ -z "$INDEX" ]; then export INDEX="/scratch/share/public_datasets/ngs/genomes_handbuilt/dkcrossm/ucsc.hg19/"`basename $BWA_DIR`"/ucsc.hg19.fa"; fi
    echo "INDEX=$INDEX"
    # get BWA/INDEX abbreviations for use in filename
    export BWA_VER=`$BWA 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`
    export SAMTOOLS_VER=`samtools 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`
    export INDEX_ABBREV=`basename $INDEX .fa | sed -e 's/[.-]/_/g;'`

    # create output file name
    export OUT_DIR=${PROJECT_DIR}/bwa/${RG_LIB}/${SAMPLE_NAME}
    mkdir -p ${OUT_DIR} 
    export OUT_SEP=-
    export OUT_NAME=${SAMPLE_NAME}${OUT_SEP}${RG_LIB}${OUT_SEP}${RG_ID}${OUT_SEP}${INDEX_ABBREV}${OUT_SEP}bwa${BWA_VER}
    # create job output dir
    JOB_DIR=${OUT_DIR}/jobs
    OJOB_DIR=${OUT_DIR}/old_jobs
    mkdir -p ${JOB_DIR} ${OJOB_DIR}
    mv ${JOB_DIR}/${TASK_NAME}-${OUT_NAME}* ${OJOB_DIR}

    # create Read Group name (TAB separated)
    export RG_NAME="@RG\tID:${RG_ID}\tLB:${RG_LIB}\tPL:ILLUMINA\tSM:${SAMPLE_NAME}\tCN:GTAC"


    # qsub
    if [ -z "$JOB_ID" ]; then
	echo -n "${OUT_NAME}:QSUB:"
	QSUB_NAME=${TASK_NAME}-${OUT_NAME}
	pushd ${OUT_DIR}
	qsub -terse \
	    ${QSUB_DRMAA} \
	    -N ${QSUB_NAME} \
	    -M $USER@uab.edu \
	    $0
	popd 
	if [ $? != 0 ]; then echo "ERROR: bad return code from QSUB"; exit 1; fi
	exit 0
    else
	echo "[debug] skipped qsub"
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
    echo "-- derrived values line params -- "
    for myvar in $DERIVED_VAR_LIST ; do
	echo -n "$myvar	:"; eval echo \$$myvar
    done
    
    echo "I'm a qsub slave: "
    #CMD="/home/curtish/ics/log_processes.sh ${PROC_LOG} "
    #echo CMD: $CMD
    #$CMD > /dev/null &

    echo "cd $OUT_DIR"
    cd ${OUT_DIR}

    # BWA aln
    echo;echo `date`"	TS	Starting BWA aln 	${OUT_NAME}"
    FWD_ENCODING=`~/ics/fastq/fastqFormatDetectM.pl -q $FWD_FASTQ`
    REV_ENCODING=`~/ics/fastq/fastqFormatDetectM.pl -q $REV_FASTQ`
    #if [[ $FWD_FASTQ == *.fastqi* || $FWD_FASTQ == *.ILM* ]]; then 
    if [ "$FWD_ENCODING" == "64" ]; then
	FWD_FMT_FLAG=-I
	echo "FWD_FASTQ: using illumina 1.3 Phred scaling (filename contains .fastqi)"
    fi
    #if [[ $REV_FASTQ == *.fastqi* || $REV_FASTQ == *.ILM* ]]; then 
    if [ "$REV_ENCODING" == "64" ]; then
	REV_FMT_FLAG=-I
	echo "REV_FASTQ: using illumina 1.3 Phred scaling (filename contains .fastqi)"
    fi
    FWD_SAI=${OUT_DIR}/${OUT_NAME}-R1.sai
    TEMP_FILE_LIST="${TEMP_FILE_LIST} ${FWD_SAI}"
    STEP_TARGET=$FWD_SAI
    STEP_NAME="BWA aln FWD"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="$BWA aln $INDEX ${FWD_FMT_FLAG} -t $NSLOTS -f ${FWD_SAI} ${FWD_FASTQ}"
	echo $CMD
	$CMD
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi
    REV_SAI=${OUT_DIR}/${OUT_NAME}-R2.sai
    TEMP_FILE_LIST="${TEMP_FILE_LIST} ${REV_SAI}"
    STEP_TARGET=$REV_SAI
    STEP_NAME="BWA aln REV"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="$BWA aln $INDEX ${REV_FMT_FLAG} -t $NSLOTS  -f ${REV_SAI} ${REV_FASTQ}"
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi
    
    # BWA sampe
    SAM_ALIGN=${OUT_DIR}/${OUT_NAME}.aligned.unsorted.sam
    TEMP_FILE_LIST="${TEMP_FILE_LIST} ${SAM_ALIGN}"
    STEP_TARGET=$SAM_ALIGN
    STEP_NAME="BWA sampe"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
	echo "STEP_TARGET: ${STEP_TARGET}"
	\ls -lst ${STEP_TARGET} ${STEP_TARGET}.done
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="$BWA sampe -r $RG_NAME -f $SAM_ALIGN $INDEX $FWD_SAI $REV_SAI $FWD_FASTQ $REV_FASTQ"
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    # SAMTools view
    BAM_ALIGN=${OUT_DIR}/${OUT_NAME}.aligned.unsorted.bam
    STEP_TARGET=$BAM_ALIGN
    STEP_NAME="samtools view"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="samtools view -bS -o ${BAM_ALIGN} ${SAM_ALIGN}"
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    # SAMTools sort
    BAM_SORTED_BASE=${OUT_DIR}/${OUT_NAME}
    BAM_SORTED=${BAM_SORTED_BASE}.bam
    STEP_TARGET=$BAM_SORTED
    STEP_NAME="samtools sort"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="samtools sort ${BAM_ALIGN} ${BAM_SORTED_BASE}"  
	echo $CMD
	$CMD
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    # SAMTools index
    BAI_SORTED=${OUT_DIR}/${OUT_NAME}.bai
    STEP_TARGET=$BAI_SORTED
    STEP_NAME="samtools index"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="samtools index ${BAM_SORTED} ${BAI_SORTED}"
	echo $CMD
	$CMD
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    # SAMTools flagstat
    FLAGSTAT_OUT=${OUT_DIR}/${OUT_NAME}.flagstat
    STEP_TARGET=$FLAGSTAT_OUT
    STEP_NAME="samtools flagstat"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="samtools flagstat ${BAM_SORTED}"
	echo $CMD
	echo "${OUT_NAME} sample_name" > ${FLAGSTAT_OUT}
	$CMD >> ${FLAGSTAT_OUT}
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    # Remove intermediate files
    echo;echo `date`"  TS      START removing intermediate files  ${OUT_NAME}"
    CMD="echo DEBUG ECHO rm -f ${TEMP_FILE_LIST}"
    echo $CMD
    $CMD
    RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR	${OUT_NAME}"; exit $RC; fi
    echo;echo `date`"  TS      DONE removing intermediate files  ${OUT_NAME}"

    # kick off GATK
#    STEP_TARGET=
    STEP_NAME="start_gatk"
#    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
#	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
#    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="/home/$USER/ics/consults/kimberly/ics223/qsub_gatk.sh -bwa ${BWA} -index ${INDEX} ${SAMPLE_NAME} ${RG_ID} ${RG_LIB}"
	echo $CMD
	(JOB_ID=; $CMD)
 	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
#	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
#    fi

    # done
    echo;echo `date`"	TS	BWA_SAM analysis completed on $SAMPLE	${OUT_NAME}"

    exit 0
fi
