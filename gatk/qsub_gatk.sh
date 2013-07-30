#!/bin/bash
#
# Take a BWA/samtools .BAM file, and run GATK re-align and SNP caller
# 
#$ -S /bin/sh
#$ -cwd
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -j y # merge stderr into stdout
# email at:  Begining, End, Abort, Suspend
#$ -m beas  
#
# *** DRMAA resources for GATK ****
#$ -pe smp 4 -l vf=5.9G -l h_vmem=6G
#
## avoid SSG nodes, which cause SNPeff failure
# -q all.q,sipsey.q
#
# *** DRMAA resources for BWA SAMPE ONLY ****
# -l vf=3.9G -l h_vmem=4G
#
# *** DRMAA resources for GATK VariantFiltration ONLY ****
# -l vf=3.9G -l h_vmem=4G
#
# ** RUN TIME ** 
#$ -l h_rt=335:00:00
#$ -l s_rt=336:55:00
#
# *** output logs ***
#$ -e jobs/$JOB_NAME.$JOB_ID.err
#$ -o jobs/$JOB_NAME.$JOB_ID.out

TASK_NAME=gatk
CMD_LINE_PARAM_LIST="PROJECT_DIR BAM_SUFFIX REGIONS.BED SLICE_SIZE"
DERIVED_VAR_LIST="CMD_LINE HOSTNAME SLICE_START SLICE_END"


#
# MASTER: submit-self on a head-node
#
if [ -z "$JOB_ID" ]; then

    # capture the original CMD_LINE
    if [ -z "$CMD_LINE" ]; then  export CMD_LINE="$0 $*"; fi

    # parse params
    while [[ "$1" == -* ]]; do 
	echo "FLAG: $1"
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
	    export REF_FASTA="$2"
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
	if [ -z "$1" ] ; then
	    echo "$myvar	: MISSING"
	    SCRIPT_NAME=`basename $0`
	    echo ""
	    echo "ERROR: $SCRIPT_NAME [-debug] $CMD_LINE_PARAM_LIST"
	    echo ""
	    echo ""
	    exit 1
	fi
	shift
    done

    # Prog & INDEX versions
    if [ -z "$BWA" ]; then export BWA="/share/apps/ngs-ccts/bwa-0.6.2/bwa"; fi
    if [ -z "${REF_FASTA}" ]; then export REF_FASTA="/scratch/share/public_datasets/ngs/genomes_handbuilt/dkcrossm/ucsc.hg19/bwa/ucsc.hg19.fa"; fi
    export PICARD_DIR="/share/apps/ngs-ccts/picard-tools-1.83/"
    export TEMP=/scratch/user/${USER}/tmp
    mkdir -p ${TEMP}
    export PUBLIC_VARIANT_DIR="/scratch/share/public_datasets/ngs/databases/gatk_bundle/2.3"
    export GATK_DIR="/share/apps/ngs-ccts/GenomeAnalysisTK-2.4-9-g532efad"
    export SNPEFF_DIR="/share/apps/ngs-ccts/snpEff_3_1"

    # get VERSIONS abbreviations for use in filename
    export BWA_VER=`$BWA 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`
    export SAMTOOLS_VER=`samtools 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`
    export REF_ABBREV=`basename ${REF_FASTA} .fa | sed -e 's/[.-]/_/g;'`

    # create BWA output file name (BAM input)
    export BWA_BAM_DIR=${PROJECT_DIR}/gatk/
    mkdir -p ${BWA_BAM_DIR} 
    export OUT_NAME=${SAMPLE_NAME}${OUT_SEP}${RG_LIB}${OUT_SEP}${RG_ID}${OUT_SEP}${REF_ABBREV}${OUT_SEP}bwa${BWA_VER}

    # create PICARD/GATK output filenames (BAM output, VCF)
    export GATK_OUT_DIR=${PROJECT_DIR}/gatk/hc${RG_LIB}/${SAMPLE_NAME}
    mkdir -p ${GATK_OUT_DIR}
    if [ -z "$GATK_RECAL_MODE" ]; then export GATK_RECAL_MODE=recal_const; fi
    export GATK_INTERVAL_NAME=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}.intervals
    export GATK_REALIGN_BAM=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}realigned.bam
    export GATK_RECAL_INIT_GRP=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}recal_initial_data.grp 
    export GATK_RECAL_BAM=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}recal.bam
    export GATK_RECAL_FINAL_GRP=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}recal_final_data.grp 
    export GATK_RECAL_REDUCED_BAM=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}recal_reduced.bam 
    export GATK_RAW_VCF=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}raw_snps_indels.vcf
    export GATK_RAW_ANNO_VCF=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}raw_snps_indels_anno.vcf
    export GATK_RECAL_SNPS=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps.recal
    export GATK_RECAL_INDEL=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}indels.recal
    export GATK_RECAL_SNPS_TRANCHES=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps.recal.tranches
    export GATK_RECAL_INDEL_TRANCHES=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}indels.recal.tranches
    export GATK_RECAL_SNPS_RSCRIPT=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps.recal.plots.R
    export GATK_RECAL_INDEL_RSCRIPT=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}indels.recal.plots.R
    export GATK_SEMIFINAL_VCF=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps_recal_filtered_ts97.vcf
    export GATK_FINAL_VCF=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps_indels_recal_filtered_ts97.vcf
    export GATK_FILT_SNPS_VCF=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}var_filt${OUT_SEP}snps_filtered.vcf
    export GATK_FILT_FINAL_VCF=${GATK_OUT_DIR}/${OUT_NAME}${OUT_SEP}gatk${GATK_VER}${OUT_SEP}var_filt${OUT_SEP}snps_indels_filtered.vcf

    # create job output dir
    JOB_DIR=${GATK_OUT_DIR}/jobs
    OJOB_DIR=${GATK_OUT_DIR}/old_jobs
    mkdir -p ${JOB_DIR} ${OJOB_DIR}
    mv ${JOB_DIR}/${TASK_NAME}-${OUT_NAME}* ${OJOB_DIR} > /dev/null 2>&1

    # qsub
    if [ -z "$JOB_ID" ]; then
	echo -n "${OUT_NAME}:QSUB:"
	pushd ${GATK_OUT_DIR} > /dev/null
	QSUB_NAME=${TASK_NAME}-${OUT_NAME}
	qsub -terse \
	    -N $QSUB_NAME \
	    -M $USER@uab.edu \
	    $0
	popd > /dev/null
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
    echo "ulimit: "`ulimit`
    echo "-- cmd line params -- "
    for myvar in $CMD_LINE_PARAM_LIST ; do
	echo -n "$myvar	:"; eval echo \$$myvar
    done
    echo "-- derrived values line params -- "
    for myvar in $DERIVED_VAR_LIST ; do
	echo -n "$myvar	:"; eval echo \$$myvar
    done
    echo "-- printenv -- "
    printenv
    echo "-- /printenv -- "

    echo "I'm a qsub slave: "
    #CMD="/home/curtish/ics/log_processes.sh ${PROC_LOG} "
    #echo CMD: $CMD
    #$CMD > /dev/null &

    #echo "cd $OUT_DIR"
    #cd ${OUT_DIR}

    #
    # check if BWA/samtools finished
    #
    BWA_DONE=${BWA_BAM_NAME}.done
    echo "CHECK: $BWA_DONE"
    if [ ! -e "${BWA_DONE}" ]; then
	echo; echo `date`"	TS	ERROR BWA_NOT_DONE	${OUT_NAME}";
	exit 1
    fi

    #
    # picard mark dulicates
    #

    # PICARD mark dup
    STEP_TARGET=${DEDUP_BAM_NAME}
    STEP_NAME="Picard MarkDup"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="java \
  	-Xmx22g \
		-Djava.io.tmpdir=$TEMP \
	-jar ${PICARD_DIR}MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT CREATE_INDEX=True TMP_DIR=$TEMP MAX_RECORDS_IN_RAM=3000000 INPUT=${BWA_BAM_NAME} OUTPUT=${DEDUP_BAM_NAME} METRICS_FILE=${DEDUP_BAM_NAME}_metrics.txt ASSUME_SORTED=True"
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    # PICARD sort
    STEP_TARGET=${DEDUP_SORT_BAM_NAME}
    STEP_NAME="Picard Sort"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="java \
  	-Xmx22g \
		-Djava.io.tmpdir=$TEMP \
        -jar ${PICARD_DIR}ReorderSam.jar CREATE_INDEX=True VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=280000000 REFERENCE=${REF_FASTA} INPUT=${DEDUP_BAM_NAME}  TMP_DIR=$TEMP  OUTPUT=${DEDUP_SORT_BAM_NAME}"
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    #
    # GATK REALIGN
    #

    STEP_TARGET=${GATK_INTERVAL_NAME}
    STEP_NAME="GATK RealignerTargetCreator"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="java -Xmx22g -Djava.io.tmpdir=$TEMP \
	-jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	-nt 4 \
	-T RealignerTargetCreator \
        -R ${REF_FASTA} \
	-known:mills,VCF $PUBLIC_VARIANT_DIR/Mills_and_1000G_gold_standard.indels.hg19.vcf \
	-known:1000g,VCF $PUBLIC_VARIANT_DIR/1000G_phase1.indels.hg19.vcf \
	-o ${GATK_INTERVAL_NAME} \
	-I ${DEDUP_SORT_BAM_NAME}"
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    STEP_TARGET=${GATK_REALIGN_BAM}
    STEP_NAME="GATK IndelRealigner"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="java -Xmx22g -Djava.io.tmpdir=$TEMP \
	-jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	-T IndelRealigner \
        -I ${DEDUP_SORT_BAM_NAME} \
        -R ${REF_FASTA} \
	-known:mills,VCF $PUBLIC_VARIANT_DIR/Mills_and_1000G_gold_standard.indels.hg19.vcf \
	-known:1000g,VCF $PUBLIC_VARIANT_DIR/1000G_phase1.indels.hg19.vcf \
	-o ${GATK_REALIGN_BAM} \
	-targetIntervals ${GATK_INTERVAL_NAME} \
        "
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    STEP_TARGET=${GATK_RECAL_INIT_GRP}
    STEP_NAME="GATK BaseRecalibrator"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="\
            java -Xmx22g -Djava.io.tmpdir=$TEMP \
	    -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	    -T BaseRecalibrator \
	    -R ${REF_FASTA} \
	    -I ${GATK_REALIGN_BAM} \
	    -knownSites:dbsnp,VCF $PUBLIC_VARIANT_DIR/dbsnp_137.hg19.vcf \
	    -knownSites:mills,VCF $PUBLIC_VARIANT_DIR/Mills_and_1000G_gold_standard.indels.hg19.vcf \
	    -knownSites:1000g,VCF $PUBLIC_VARIANT_DIR/1000G_phase1.indels.hg19.vcf \
	    -cov ReadGroupCovariate \
	    -cov QualityScoreCovariate \
	    -cov CycleCovariate \
	    -cov ContextCovariate \
	    -o ${GATK_RECAL_INIT_GRP} \
        "          
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    STEP_TARGET=${GATK_RECAL_BAM}
    STEP_NAME="GATK PrintReads"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD="\
            java -Xmx22g -Djava.io.tmpdir=$TEMP \
	    -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	    -T PrintReads \
            -R ${REF_FASTA} \
	    -I ${GATK_REALIGN_BAM} \
	    -BQSR ${GATK_RECAL_INIT_GRP} \
    	    -o ${GATK_RECAL_BAM} \
        " 
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    #
    # Samtools index on the recalibrated bam file???
    #

    #
    # ReduceReads  (not used down stream)
    #
    STEP_TARGET=${GATK_RECAL_REDUCED_BAM}
    STEP_NAME="GATK ReduceReads"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD=" \
	    java -Xmx22g -Djava.io.tmpdir=$TEMP \
	    -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	    -T ReduceReads \
            --dont_compress_read_names \
            -R ${REF_FASTA} \
	    -I ${GATK_RECAL_BAM} \
	    -o ${STEP_TARGET} \
        "
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    #
    # GATK Recalibration
    #
    STEP_TARGET=${GATK_RECAL_FINAL_GRP}
    STEP_NAME="GATK Final BaseRecalibrator"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	CMD=" \
	    java -Xmx22g -Djava.io.tmpdir=$TEMP \
	    -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	    -T BaseRecalibrator \
	    -R ${REF_FASTA} \
	    -I ${GATK_RECAL_BAM} \
	    -knownSites:dbsnp,VCF $PUBLIC_VARIANT_DIR/dbsnp_137.hg19.vcf \
	    -knownSites:mills,VCF $PUBLIC_VARIANT_DIR/Mills_and_1000G_gold_standard.indels.hg19.vcf \
	    -knownSites:1000g,VCF $PUBLIC_VARIANT_DIR/1000G_phase1.indels.hg19.vcf \
	    -cov ReadGroupCovariate \
	    -cov QualityScoreCovariate \
	    -cov CycleCovariate \
	    -cov ContextCovariate \
	    -o ${GATK_RECAL_FINAL_GRP} \
        "
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    # 
    # GATK HaplotypeCaller (single thread) 
    #
    STEP_TARGET=${GATK_RAW_VCF}
    STEP_NAME="GATK HaplotypeCaller"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	# GATK2.4-9-g532efad  : HaplotypeCaller supports NEITHER -nt nor -nct
	CMD=" \
        java -Xmx22g -Djava.io.tmpdir=$TEMP \
	-jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R ${REF_FASTA} \
	-I ${GATK_RECAL_BAM} \
	--dbsnp $PUBLIC_VARIANT_DIR/dbsnp_137.hg19.vcf \
	-o ${STEP_TARGET} \
        "
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi

    STEP_TARGET=${GATK_RAW_ANNO_VCF}
    STEP_NAME="SNPSift"
    if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
    else
	echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	# note: using Xmx22g to match GATK, but not needed
	# WARNING: DO NOT USE -Xmx - if you do both java -Xmx, qsub -l h_vmem, and "annotate" then snpEff will fail:
	# "I/O problem while mapping file '/scratch/share/public_datasets/ngs/databases/gatk_bundle/2.3/dbsnp_137.hg19.vcf'"
	# annotate | annMem
	# dbsnp_137.hg19 uses (annotate: ~12.1G; annMem: ~12.8G)
	CMD=" \
            java \
            \
            -Djava.io.tmpdir=$TEMP \
	    -jar ${SNPEFF_DIR}/SnpSift.jar \
	    annMem \
	    -id ${PUBLIC_VARIANT_DIR}/dbsnp_137.hg19.vcf \
	    ${GATK_RAW_VCF} \
        "
	echo $CMD "> ${STEP_TARGET}"
	$CMD > ${STEP_TARGET}
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi


    #
    # VariantRecalibrator needs to be run differenntly for "small samples"
    # http://gatkforums.broadinstitute.org/discussion/1927/variantrecalibrator-error-stack-trace
    #
    if [ 0 == 1 ]; then 
	STEP_TARGET=${GATK_RECAL_FILE}
	STEP_NAME="GATK VariantRecalibrator"
	if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	    echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
	else
	    echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	    CMD=" \
	    java -Xmx22g -Djava.io.tmpdir=$TEMP \
	    -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	    -T VariantRecalibrator \
	    -R ${REF_FASTA} \
	    -input ${GATK_RAW_ANNO_VCF} \
	    --maxGaussians 4 \
	    -percentBad 0.05 \
	    -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ${PUBLIC_VARIANT_DIR}/hapmap_3.3.hg19.vcf \
	    -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ${PUBLIC_VARIANT_DIR}/1000G_omni2.5.hg19.vcf \
	    -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=6.0 ${PUBLIC_VARIANT_DIR}/dbsnp_137.hg19.vcf \
	    -resource:mills,VCF,known=false,training=true,truth=true,prior=12.0 ${PUBLIC_VARIANT_DIR}/Mills_and_1000G_gold_standard.indels.hg19.vcf \
	    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP -an ClippingRankSum \
	    -mode BOTH \
	    -recalFile ${GATK_RECAL_FILE} \
	    -tranchesFile ${GATK_RECAL_TRANCHES} \
	    -rscriptFile ${GATK_RECAL_RSCRIPT} \
	    "
	    echo $CMD
	    $CMD 
	    RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	    touch ${STEP_TARGET}.done
	    echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
	fi

	# Apply Recalibration to SNPs/indels
	STEP_TARGET=${GATK_FINAL_VCF}
	STEP_NAME="GATK ApplyRecalibration"
	if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	    echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
	else
	    echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	    CMD=" \
	    java -Xmx22g -Djava.io.tmpdir=$TEMP \
	    -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	    -T ApplyRecalibration \
	    -R ${REF_FASTA} \
	    -input ${GATK_RAW_ANNO_VCF} \
	    -recalFile ${GATK_RECAL_FILE} \
	    -tranchesFile ${GATK_RECAL_TRANCHES} \
	    --ts_filter_level 97.0 \
	    -mode BOTH \
	    -o ${STEP_TARGET} \
	    "
	    echo $CMD
	    $CMD 
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

	# done
	echo;echo `date`"	TS	Picard/GATK analysis completed on $SAMPLE	${OUT_NAME}"
    fi

    #
    # SMALL SAMPLE VariantRecalibrator
    # based on comments here 
    # http://gatkforums.broadinstitute.org/discussion/1927/variantrecalibrator-error-stack-trace
    # See more at: http://gatkforums.broadinstitute.org/discussion/1186/best-practice-variant-detection-with-the-gatk-v4-for-release-2-0#sthash.07RzWE0U.dpuf
    #
    if [ 1 == 1 ]; then 
	# Filter SNPs
	STEP_TARGET=${GATK_FILT_SNPS_VCF}
	STEP_NAME="GATK VariantFiltration SNPS"
	if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	    echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
	else
	    echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	    CMD=' 
	    java -Xmx3g -Djava.io.tmpdir='$TEMP' 
	    -jar '${GATK_DIR}'/GenomeAnalysisTK.jar 
	    -T VariantFiltration 
	    -R '${REF_FASTA}' 
	    --variant '${GATK_RAW_ANNO_VCF}' 
	    --filterName "gatkExomeSnpsQDlt2" --filterExpression "QD < 2.0 && vc.isSNP()" 
	    --filterName "gatkExomeSnpsMQlt40" --filterExpression "MQ < 40.0 && vc.isSNP()" 
	    --filterName "gatkExomeSnpsFSgt60" --filterExpression "FS > 60.0 && vc.isSNP()" 
	    --filterName "gatkExomeSnpsHSgt13" --filterExpression "HaplotypeScore > 13.0 && vc.isSNP()" 
	    --filterName "gatkExomeSnpsMQRSltNeg12p5" --filterExpression "MQRankSum < -12.5 && vc.isSNP()" 
	    --filterName "gatkExomeSnpsRPRSltNeg8" --filterExpression "ReadPosRankSum < -8.0 && vc.isSNP()" 
            --out '${STEP_TARGET}' 
	    '
	    echo $CMD
	    eval $CMD 
	    RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	    touch ${STEP_TARGET}.done
	    echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
	fi

	# filter INDELs
	# omitted "InbreadingCoeff < -0.8" as I'm doing this per-sample.
	STEP_TARGET=${GATK_FILT_FINAL_VCF}
	STEP_NAME="GATK VariantFiltration INDELS"
	if [[ ( -e "${STEP_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${STEP_TARGET}.done" ) ]]; then
	    echo;echo `date`"	TS	SKIP ${STEP_NAME}	${OUT_NAME}"; 
	else
	    echo;echo `date`"	TS	START ${STEP_NAME}	${OUT_NAME}"; 
	    CMD=" \
	    java -Xmx3g -Djava.io.tmpdir=$TEMP \
	    -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	    -T VariantFiltration \
	    -R ${REF_FASTA} \
	    --variant ${GATK_FILT_SNPS_VCF} \
	    --filterName \"gatkExomeIndelQDlt2\" --filterExpression \"QD < 2.0 && not vc.isSNP()\" \
	    --filterName \"gatkExomeIndelRPRSltNeg20\" --filterExpression \"ReadPosRankSum < -20.0 && not vc.isSNP()\" \
	    --filterName \"gatkExomeIndelFSgt200\" --filterExpression \"FS > 200.0 && not vc.isSNP()\" \
            --out ${STEP_TARGET} \
	    "
	    echo $CMD
	    eval $CMD 
	    RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC; fi
	    touch ${STEP_TARGET}.done
	    echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
	fi
    fi
    
    exit 0
fi
