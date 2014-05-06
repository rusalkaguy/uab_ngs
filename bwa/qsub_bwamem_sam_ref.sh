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
## -l h_rt=119:00:00 -l s_rt=120:55:00 # 5 days
## -l h_rt=48:00:00 -l s_rt=47:55:00 # 2 days
#$ -l h_rt=12:00:00 -l s_rt=11:55:00 # 1/2 day
#
# *** output logs ***
#$ -e jobs/$JOB_NAME.$JOB_ID.err
#$ -o jobs/$JOB_NAME.$JOB_ID.out


TASK_NAME=bwa_mem
CMD_LINE_PARAM_LIST="WORK_DIR SAMPLE_NAME FWD_FASTQ REV_FASTQ REF_FASTA"
DERIVED_VAR_LIST="CMD_LINE HOSTNAME PROJECT_DIR REV_FASTQ DONE_ONLY QSUB_PE_OVERRIDE"

QSUB_DRMAA="-pe smp 4 -l vf=1.9G -l h_vmem=2G" 
if [ -z "$NSLOTS" ]; then export NSLOTS=4; fi  # for debugging

export TEMP=/scratch/user/$USER/tmp
# load needed modules 
# we hit the exe directly
module load ngs-ccts/bwa/0.7.7
module load ngs-ccts/samtools/0.1.19
if [ -z "$PICARD_DIR" ]; then PICARD_DIR=/share/apps/ngs-ccts/picard-tools/picard-tools-1.110; fi
module load R/R-3.0.1  # for Picard to make PDF
if [ -z "$IGVTOOLS_JAR" ]; then IGVTOOLS=/share/apps/ngs-ccts/IGVTools/IGVTools_2.3.32/igvtools.jar; fi
if [ -z "$GATK_JAR" ]; then GATK_JAR=/share/apps/ngs-ccts/GenomeAnalysisTK/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar; fi

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
    run_step $SAMPLE_NAME $SAM "BWA_mem" $SAM \
	$BWA mem -t $NSLOTS \
	-M -C \
	-R "'@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}'" \
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
    JAVA_DRMAA="-Xms${JAVA_RAM}g -Xmx${JAVA_RAM}g  -Djava.io.tmpdir=$TEMP" 
    OUT_FILE="${BAM_BASE}.fragstat"
    OUT_HIST="${BAM_BASE}.fraghist.pdf"
    run_step $SAMPLE_NAME $OUT_FILE PICARD_CollectInsertSizeMetrics - \
	java $JAVA_DRMAA \
	-Djava.io.tmpdir='/scratch/share/galaxy/temp' \
	-jar ${PICARD_DIR}/CollectInsertSizeMetrics.jar \
	VALIDATION_STRINGENCY=LENIENT \
	ASSUME_SORTED=true \
	INPUT=$BAM \
	OUTPUT=$OUT_FILE \
	HISTOGRAM_FILE=$OUT_HIST 

    # PICARD de-dup
    BAM_DEDUP=${BAM_BASE}.dedup.bam
    BAI_DEDUP=${BAM_BASE}.dedup.bai
    DEDUP_MATRICS=${BAM_BASE}.dedup.metrics.txt
    run_step $SAMPLE_NAME $BAM_DEDUP PICARD_dedup - \
	java $JAVA_DRMAA \
	-jar ${PICARD_DIR}/MarkDuplicates.jar \
	VALIDATION_STRINGENCY=SILENT \
	CREATE_INDEX=True \
	TMP_DIR=$TEMP \
	MAX_RECORDS_IN_RAM=3000000 \
	INPUT=${BAM} \
	OUTPUT=${BAM_DEDUP} \
	METRICS_FILE=${DEDUP_METRICS} \
	ASSUME_SORTED=True
    # SAMTools index
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    run_step $SAMPLE_NAME $BAI_DEDUP "samtools_index" - \
	samtools index ${BAM_DEDUP} ${BAI_DEDUP}

    # IGVTools Coverage
    COV_W=10 # window 
    COV_REF=~/igv/genomes/hcmvMerlin.genome
    if [[ "$INDEX_ABBREV" == *TR-BAC* ]]; then COV_REF=~/igv/genomes/hcmvTR-BAC.genome; fi
    OUT_COV="${BAM_BASE}.cov.w${COV_W}.wig"
    run_step $SAMPLE_NAME $OUT_COV IGVTOOLS_count - \
	java $JAVA_DRMAA \
	-Djava.awt.headless=true -jar $IGVTOOLS \
	count \
	-z 10 -w $COV_W -e 0 \
	${BAM} \
	${OUT_COV} \
	${COV_REF}

    # GATK DepthOfCoverage (requires RG)
    OUT_DOC=${BAM_BASE}.gatk_doc
    run_step $SAMPLE_NAME $OUT_DOC GATK_DepthOfCoverage - \
	java $JAVA_DRMAA \
	-jar $GATK_JAR \
	-T DepthOfCoverage \
	-R ${REF_FASTA} \
	-I ${BAM} \
	-ct 10 -ct 100 -ct 250 -ct 500 -ct 1000 \
	--out $OUT_DOC

    # GATK re-align targets
    GATK_REALIGN_INTERVALS=${BAM_BASE}.dedup.intervals
    run_step $SAMPLE_NAME $GATK_REALIGN_INTERVALS GATK_realign_create_targets - \
	java $JAVA_DRMAA \
	-jar ${GATK_JAR} \
	-nt $NSLOTS \
	-T RealignerTargetCreator \
        -R ${REF_FASTA} \
	-I ${BAM_DEDUP} \
	-o ${GATK_REALIGN_INTERVALS} 
    BAM_REALIGN=${BAM_BASE}.dedup.realign.bam
    BAI_REALIGN=${BAM_BASE}.dedup.realign.bai
    run_step $SAMPLE_NAME $BAM_REALIGN GATK_realign - \
	java $JAVA_DRMAA \
	-jar ${GATK_JAR} \
	-T IndelRealigner \
        -R ${REF_FASTA} \
	-I ${BAM_DEDUP} \
	-targetIntervals ${GATK_REALIGN_INTERVALS} \
	-o ${BAM_REALIGN}
    # SAMTools index
    run_step $SAMPLE_NAME $BAI_REALIGN "samtools_index" - \
	samtools index ${BAM_REALIGN} ${BAI_REALIGN}

    if [ 1 == 0 ]; then
        # GATK recalibrator
	RECAL_BQSR=${BAM_BASE}.dedup.realign.bqsr
	BAM_RECAL=${BAM_BASE}.dedup.realign.recal.bam
	BAI_RECAL=${BAM_BASE}.dedup.realign.recal.bai
	run_step $SAMPLE_NAME $RECAL_BQSR GATK_recalibrator - \
	    java $JAVA_DRMAA \
	    -jar ${GATK_JAR} \
	    -T BaseRecalibrator \
	    -cov QualityScoreCovariate \
	    -cov CycleCovariate \
	    -cov ContextCovariate \
	    -R ${REF_FASTA} \
	    -I ${BAM_REALIGN} \
	    -o ${BAM_RECAL}
	run_step $SAMPLE_NAME $RECAL_BQSR GATK_printReads - \
	    java $JAVA_DRMAA \
	    -jar ${GATK_JAR} \
	    -T PrintReads \
	    -R ${REF_FASTA} \
	    -I ${BAM_REALIGN} \
	    -BQSR $RECAL_BQSR \
	    -o ${BAM_RECAL}
	# SAMTools index
	run_step $SAMPLE_NAME $BAI_RECAL "samtools_index" - \
	    samtools index ${BAM_RECAL} ${BAI_RECAL}
    fi

    # GATK UnifiedGenotype (VCF)
    UG_P1_VCF=${BAM_BASE}.p1.vcf
    # can also use -nct : better because less mem? 
    # --output_mode EMIT_ALL_SITES \
    run_step $SAMPLE_NAME $UG_P1_VCF GATK_UG - \
	java $JAVA_DRMAA \
	-jar ${GATK_JAR} \
	-nt ${NSLOTS} \
	-T UnifiedGenotyper \
	-R ${REF_FASTA} \
	-ploidy 1 \
	-glm BOTH \
	-stand_call_conf 20 \
	-stand_emit_conf 20 \
	-I ${BAM_REALIGN} \
	-o ${UG_P1_VCF}

    UG_P10_VCF=${BAM_BASE}.p10.vcf
    UG_P10AF_VCF=${BAM_BASE}.p10.af.vcf
    # can also use -nct : better because less mem? 
    run_step $SAMPLE_NAME $UG_P10_VCF GATK_UG - \
	java $JAVA_DRMAA \
	-jar ${GATK_JAR} \
	-nt ${NSLOTS} \
	-T UnifiedGenotyper \
	-R ${REF_FASTA} \
	-ploidy 10 \
	-stand_emit_conf 20 \
	-glm BOTH \
	-I ${BAM_REALIGN} \
	-o ${UG_P10_VCF}
    # filter to just partialy penatrant variants
    run_step $SAMPLE_NAME $UG_P10AF_VCF grep_AF_not_1  $UG_P10AF_VCF \
        grep -v "AF=1.00" $UG_P10_VCF 

    UG_P20_VCF=${BAM_BASE}.p20.vcf
    UG_P20AF_VCF=${BAM_BASE}.p20.af.vcf
    # can also use -nct : better because less mem? 
    run_step $SAMPLE_NAME $UG_P20_VCF GATK_UG - \
	java $JAVA_DRMAA \
	-jar ${GATK_JAR} \
	-nt ${NSLOTS} \
	-T UnifiedGenotyper \
	-R ${REF_FASTA} \
	-ploidy 20 \
	-stand_emit_conf 20 \
	-glm BOTH \
	-I ${BAM_REALIGN} \
	-o ${UG_P20_VCF}
    # filter to just partialy penatrant variants
    run_step $SAMPLE_NAME $UG_P20AF_VCF grep_AF_not_1 $UG_P20AF_VCF \
        grep -v "AF=1.00" $UG_P20_VCF 
	    

    # GATK ReferenceMaker (fasta)
    CON_NAME=${SAMPLE_NAME}.consensus
    CON_FA=${OUT_DIR}/${CON_NAME}.fa
    run_step $SAMPLE_NAME $CON_FA GATK_ReferenceMaker - \
	java $JAVA_DRMAA \
	-jar ${GATK_JAR} \
	-T FastaAlternateReferenceMaker \
	-R ${REF_FASTA} \
	--variant ${UG_P1_VCF} \
	-o ${CON_FA}
    # BWA index CONSENSUS
    CON_SIZE=`wc -c ${CON_FA} | cut -d " " -f 1`
    export CON_BWA_INDEX_TYPE='is'
    CON_BWA_MAX_IS=$(( 2**30 ))
    if [ $CON_SIZE -gt $CON_BWA_MAX_IS ]; then
	export CON_BWA_INDEX_TYPE='bwtsw'
    fi
    CONSENSUS_BWA_IDX=${OUT_DIR}/${CON_NAME}.fa.bwa${BWA_VER}.bwt
    run_step $SAMPLE_NAME $CONSENSUS_BWA_IDX.bwt BWA_index_CON_fa - \
	bwa index -a $CON_BWA_INDEX_TYPE -p $CONSENSUS_BWA_IDX $CON_FA    

    # 
    # re-align reads to consensus
    #
    CON_SAM=${OUT_DIR}/${CON_NAME}.sam.tgz
    CON_UBAM=${OUT_DIR}/${CON_NAME}.unsorted.bam
    CON_BAM_BASE=${OUT_DIR}/${CON_NAME}
    CON_BAM=${CON_BAM_BASE}.bam
    CON_BAI=${CON_BAM_BASE}.bai

    run_step $SAMPLE_NAME $CON_SAM "BWA_mem_vs_CONSENSUS" $CON_SAM \
	$BWA mem -t $NSLOTS \
	-M -C \
	-R "'@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}'" \
	$CONSENSUS_BWA_IDX \
	$FWD_FASTQ $REV_FASTQ \
	\| gzip 


    # .SAM.gz to sorted .BAM
    run_step $SAMPLE_NAME $CON_UBAM "samtools_view_CONSENSUS" $CON_UBAM \
	zcat $CON_SAM \|  \
	sed "'s/\([1-9]\:N\:[0-9][0-9]*\:[A-Z]*\)/BC\:Z\:\1/'" \| \
	samtools view -bS -@ $NSLOTS -
    run_step $SAMPLE_NAME $CON_BAM "samtools_sort_CON" $CON_BAM \
	  samtools sort -@ $NSLOTS $CON_UBAM $CON_BAM_BASE


    # SAMTools index
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    run_step $SAMPLE_NAME $CON_BAI "samtools_index" - \
	samtools index ${CON_BAM} ${CON_BAI}

    # SAMTools flagstat
    # run_step(SAMPLE_NAME,TARGET,STEP_NAME,STDOUT,cmd...)
    CON_FLAGSTAT_OUT=${OUT_DIR}/${CON_NAME}.flagstat
    run_step $SAMPLE_NAME $CON_FLAGSTAT_OUT "samtools_flagstat_CON" $CON_FLAGSTAT_OUT \
	samtools flagstat ${CON_BAM}



    #
    # cleanup - if requested
    #
    if [ -n "$CLEAN" ]; then 
	run_cmd - \
	    rm -rf $SAM $UBAM $BAM $BAI \
	    $CON_SAM $CON_UBAM \
	    $BAM_DEDUP  $BAI_DEDUP \

    fi

    exit 0
fi
