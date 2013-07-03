#!/bin/bash
#
# Find all .bams and run them through GATK HC, slicing by probe region
# 
#$ -S /bin/sh
#$ -cwd
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -j y
# email at:  Begining, End, Abort, Suspend
#$ -m beas  
#
## avoid SSG nodes, which cause SNPeff failure
#$ -q all.q,sipsey.q
#
# *** DRMAA resources for GATK HAPLOTYPE CALLER(single threaded) ****
# -l vf=24G -l h_vmem=23G
#
# ** RUN TIME ** 
#$ -l h_rt=335:00:00
#$ -l s_rt=336:55:00
#
# *** output logs ***
#$ -e jobs/$JOB_NAME.$JOB_ID.err
#$ -o jobs/$JOB_NAME.$JOB_ID.out

TASK_NAME=gatkHC
CMD_LINE_PARAM_LIST="BAM_DIR BAM_PATTERN"
DERIVED_VAR_LIST="HOSTNAME CMD_LINE RUN_TAG SAFE_DIR SAFE_PAT OUT_NAME OUT_DIR GATK_RAW_VCF BAM_FLIST GATK_VER INTERVAL_BED INTERVAL_PAD SLICE_MODE SLICE QSUB_RAM QSUB_PE NSLOTS" 

export QSUB_RAM=24
export QSUB_PE=;
export PROJECT_DIR=/scratch/user/${USER}/kimberly
if [ -z "$INTERVAL_BED" ]; then export INTERVAL_BED=/home/curtish/ics/consults/kimberly/ics223/ESRD_SLE_TDS_Files/hg19_MyGen_SLE_TMA_CHIP.bed; fi
if [ -z "$INTERVAL_PAD" ]; then export INTERVAL_PAD=0; fi
if [ -z "$SLICE_MODE" ]; then export SLICE_MODE=all; fi

# DEFAULT PATHS
export BWA="/share/apps/ngs-ccts/bwa-0.6.2/bwa"
#export BWA="/share/apps/ngs-ccts/bwa-0.5.9/bwa"
export REF_FASTA="/scratch/share/public_datasets/ngs/genomes_handbuilt/dkcrossm/ucsc.hg19/bwa/ucsc.hg19.fa"
export PICARD_DIR="/share/apps/ngs-ccts/picard-tools-1.83/"
export TEMP=/scratch/user/${USER}/tmp
export GATK_DIR="/share/apps/ngs-ccts/GenomeAnalysisTK-2.5-2-gf57256b"
export QSUB=qsub

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
	# check for TEST (echo instead of qsub)
	if [ "-test" == "$1" ]; then
	    export QSUB="echo qsub"
	    echo "**** TEST MODE ON ****" 
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
	    export QSUB_PE="$1 $2 $3"
	    export NSLOTS=$3
	    echo "**** -PE OVERRIDE=$QSUB_PE **** ) " 
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
	# check for -run TAG - a name for this 
	if [ "-tag" == "$1" ]; then
	    export RUN_TAG="-$2"
	    echo "**** -tag=$RUN_TAG **** " 
	    shift 1
	    shift 1
	    continue
	fi
	# check for -slice mode - run one node per chromosome
	if [ "-slice" == "$1" ]; then
	    export SLICE_MODE="$2"
	    echo "**** -slice=$SLICE_MODE **** " 
	    shift 1
	    shift 1
	    continue
	fi
	# check for -slice mode - run one node per chromosome
	if [ "-pad" == "$1" ]; then
	    export INTERVAL_PAD="$2"
	    echo "**** -pad=$INTERVAL_PAD **** " 
	    shift 1
	    shift 1
	    continue
	fi
	# check for -slice mode - run one node per chromosome
	if [ "-ram" == "$1" ]; then
	    export QSUB_RAM="$2"
	    echo "**** -ram=$QSUB_RAM **** " 
	    shift 1
	    shift 1
	    continue
	fi
	# unknown flag
	echo "ERROR: unknown option: $1"
	exit 1 > /dev/null
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
	    echo "ERROR: $SCRIPT_NAME [-debug] ${CMD_LINE_PARAM_LIST}"
	    echo " "
	    echo "OPTIONAL FLAGS "
	    echo " -debug            [default: ] - run w/o qsub"
	    echo " -done             [default: $DONE_ONLY] - check only .done files"
	    echo " -tag   RUN_TAG    [default: $RUN_TAG] - add tag to dir/run name"
	    echo " -slice SLICE_MODE [default: $SLICE_MODE] - slice by chr per node"
	    echo " -pad   INTERVAL_PAD [default: $INTERVAL_PAD] - bp to scan outside target intervals"
	    echo " -pe    smp N      [default: $PE_DRMAA] - qsub with parallel env"
	    echo " -bwa   EXE_PATH   [default: $BWA] - path to exe"
	    echo " -index REF_FASTA  [default: $REF_FASTA]"
	    echo " -ram   QSUB_RAM   [default: $QSUB_RAM] - in Gig"
	    echo " "
	    echo " "
	    echo " "
	    echo ""
	    exit 1 > /dev/null
	fi
	shift
    done

    # Prog & INDEX versions
    mkdir -p ${TEMP}
    export REF_ABBREV=`basename ${REF_FASTA} .fa | sed -e 's/[.-]/_/g;'`
    export GATK_VER=`basename ${GATK_DIR} | cut -d - -f 2-3`
    export GATK_SHORT_VER=`basename ${GATK_DIR} | cut -d - -f 2`
    export PUBLIC_VARIANT_DIR="/scratch/share/public_datasets/ngs/databases/gatk_bundle/${GATK_SHORT_VER}/hg19"  # could use REF_ABBREV?
    #export SNPEFF_DIR="/share/apps/ngs-ccts/snpEff_3_1"
    export SNPEFF_DIR="/share/apps/ngs-ccts/snpEff_3_2a"

    # get VERSIONS abbreviations for use in filename
    export BWA_VER=`$BWA 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`
    export SAMTOOLS_VER=`samtools 2>&1 | grep "^Version" | cut -d " " -f 2  | sed -e 's/[.-]/_/g;'`

    # create BWA output file name (BAM input)
    export SAFE_PAT=`echo -n "${BAM_PATTERN}" | perl -pe 's/[^a-z0-9._-]/./gi;'`
    echo "SAFE_PAT=$SAFE_PAT"
    export SAFE_DIR=`basename $BAM_DIR` 
    echo "SAFE_DIR=$SAFE_DIR"
    export SLICE_TAG=; if [ "all" != "$SLICE_MODE" ]; then export SLICE_TAG=-slice_$SLICE_MODE; fi
    export PAD_TAG=; if [ "0" != "$INTERVAL_PAD" ]; then export PAD_TAG=-pad_$INTERVAL_PAD; fi
    export OUT_NAME=${SAFE_DIR}-${SAFE_PAT}${SLICE_TAG}${PAD_TAG}${RUN_TAG}
    export OUT_DIR=$PROJECT_DIR/gatk/hc/${OUT_NAME}
    mkdir -p ${OUT_DIR} 
    
    #
    # create groups of input bams per cohort
    #
    export BAM_ALL_LIST=${OUT_DIR}/bams-${SAFE_DIR}-${SAFE_PAT}-all.list
    echo "BAM_ALL_LIST=$BAM_ALL_LIST"
    #CMD="lfs find $BAM_DIR -name \"*${BAM_PATTERN}*.bam\" > $BAM_ALL_LIST"
    #echo $CMD

    FIND_CMD="find"; which lfs 2>/dev/null >/dev/null; if [ $? == 0 ]; then FIND_CMD="lfs find"; fi
    $FIND_CMD $BAM_DIR -name "*${BAM_PATTERN}*.bam" > $BAM_ALL_LIST

    export BAM_HEALTHY_TAG=healthy
    export BAM_HEALTHY_LIST=${OUT_DIR}/bams-${SAFE_DIR}-${SAFE_PAT}-${BAM_HEALTHY_TAG}.list
    grep -i healthy $BAM_ALL_LIST > $BAM_HEALTHY_LIST

    export BAM_SLE_ONLY_TAG=sle_only
    export BAM_SLE_ONLY_LIST=${OUT_DIR}/bams-${SAFE_DIR}-${SAFE_PAT}-${BAM_SLE_ONLY_TAG}.list
    grep -i non_nephritis_  $BAM_ALL_LIST > $BAM_SLE_ONLY_LIST

    export BAM_ESRD_TAG=esrd
    export BAM_ESRD_LIST=${OUT_DIR}/bams-${SAFE_DIR}-${SAFE_PAT}-${BAM_ESRD_TAG}.list
    grep -i esrd_sle_  $BAM_ALL_LIST > $BAM_ESRD_LIST

    export BAM_HAPMAP_TAG=hapmap
    export BAM_HAPMAP_LIST=${OUT_DIR}/bams-${SAFE_DIR}-${SAFE_PAT}-${BAM_HAPMAP_TAG}.list
    grep -i hapmap_  $BAM_ALL_LIST > $BAM_HAPMAP_LIST

    echo -n "BAMs found matching ${BAM_PATTERN}"
    wc -l \
	$BAM_ALL_LIST 
    wc -l \
	$BAM_HEALTHY_LIST \
	$BAM_SLE_ONLY_LIST \
	$BAM_ESRD_LIST \
	$BAM_HAPMAP_LIST
	
    #
    # create region slices for paralellization
    #
    SLICE_LIST=all
    if [ "chr" == "$SLICE_MODE" ]; then 
	echo "SLICE_MODE=${SLICE_MODE} ****** spliting by chromosome "
	export GATK_SLICE_BASE=${OUT_DIR}/slice-${OUT_NAME}
	# get list of uniq chromosomes
	export SLICE_LIST=`cut -f 1 ${INTERVAL_BED} | sort | uniq`
	for SLICE in $SLICE_LIST; do
	    grep "^$SLICE[[:space:]]" ${INTERVAL_BED} > ${GATK_SLICE_BASE}-${SLICE}.bed 
	done
	echo -n "Slices "`\ls -1 ${GATK_SLICE_BASE}-*.bed | wc -l` 
	wc -l ${GATK_SLICE_BASE}-*.bed
    fi

    # 
    # ram 
    # 
    export QSUB_VF=$QSUB_RAM
    #export QSUB_VMEM=$(($QSUB_RAM - 1))".5"
    export QSUB_VMEM=$QSUB_RAM
    export JVM_RAM=$(($QSUB_RAM - 1))

    # create PICARD/GATK output filenames (BAM output, VCF)
    if [ -z "$GATK_RECAL_MODE" ]; then export GATK_RECAL_MODE=recal_const; fi
    export OUT_SEP=-

    # create job output dir
    JOB_DIR=${OUT_DIR}/jobs
    OJOB_DIR=${OUT_DIR}/old_jobs
    mkdir -p ${JOB_DIR} ${OJOB_DIR}
    mv ${JOB_DIR}/${TASK_NAME}-${OUT_NAME}* ${OJOB_DIR} > /dev/null 2>&1

    # qsub
#    if [ -z "$JOB_ID" ]; then
	pushd ${OUT_DIR} > /dev/null
	# iterate over slices (or just all)
	for SLICE in $SLICE_LIST; do 
	    if [ "$SLICE" == "all" ]; then 
		echo -n "${OUT_NAME}:QSUB:"
		export JOB_NAME=${TASK_NAME}-${OUT_NAME}
	    else
		echo -n "${OUT_NAME}-${SLICE}:QSUB:"
		export INTERVAL_BED=${GATK_SLICE_BASE}-${SLICE}.bed
		export JOB_NAME=${TASK_NAME}-${OUT_NAME}-${SLICE}
	    fi
	    # per-slice VCF/filenames
	    export GATK_RAW_VCF=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}raw_snps_indels${OUT_SEP}${SLICE}.vcf
	    export GATK_RAW_ANNO_VCF=${OUT_DIR}gatk${GATK_VER}${OUT_SEP}raw_snps_indels_anno${OUT_SEP}${SLICE}.vcf
	    export GATK_RECAL_SNPS=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps${OUT_SEP}${SLICE}.recal
	    export GATK_RECAL_INDEL=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}indels${OUT_SEP}${SLICE}.recal
	    export GATK_RECAL_SNPS_TRANCHES=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps.recal${OUT_SEP}${SLICE}.tranches
	    export GATK_RECAL_INDEL_TRANCHES=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}indels.recal${OUT_SEP}${SLICE}.tranches
	    export GATK_RECAL_SNPS_RSCRIPT=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps.recal.plots${OUT_SEP}${SLICE}.R
	    export GATK_RECAL_INDEL_RSCRIPT=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}indels.recal.plots${OUT_SEP}${SLICE}.R
	    export GATK_SEMIFINAL_VCF=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps_recal_filtered_ts97${OUT_SEP}${SLICE}.vcf
	    export GATK_FINAL_VCF=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}${GATK_RECAL_MODE}${OUT_SEP}snps_indels_recal_filtered_ts97${OUT_SEP}${SLICE}.vcf
	    export GATK_FILT_SNPS_VCF=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}var_filt${OUT_SEP}snps_filtered${OUT_SEP}${SLICE}.vcf
	    export GATK_FILT_FINAL_VCF=${OUT_DIR}/gatk${GATK_VER}${OUT_SEP}var_filt${OUT_SEP}snps_indels_filtered${OUT_SEP}${SLICE}.vcf
	    if [ "$JOB_ID" != "run_now" ]; then
		# qsub!
		MEM_DRMAA="-l vf=${QSUB_VF}G -l h_vmem=${QSUB_VMEM}G" 
		MEM_DRMAA="-l vf=${QSUB_VF}G"
		${QSUB} -terse \
		    $PE_DRMAA $MEM_DRMAA \
		    -N $JOB_NAME \
		    -M $USER@uab.edu \
		    -v SLICE=$SLICE \
		    -v INTERVAL_BED=$INTERVAL_BED \
		    -v GATK_RAW_VCF=$GATK_RAW_VCF \
		    $0
	    else
		echo $JOB_ID "(in-line, no qsub)"
		export SLICE=$SLICE
		export INTERVAL_BED=$INTERVAL_BED
		export GATK_RAW_VCF=$GATK_RAW_VCF
		$0  # call self as slave
	    fi
	done
	popd > /dev/null
	if [ $? != 0 ]; then echo "ERROR: bad return code from QSUB"; exit 1 > /dev/null; fi
	exit 0 > /dev/null
#    else
#	echo "[debug] skipped qsub"
#    fi
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
    echo "CMD: printenv > jobs/$JOB_NAME.$JOB_ID.env"
    printenv > jobs/$JOB_NAME.$JOB_ID.env

    echo "I'm a qsub slave: "

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
        java -Xmx${JVM_RAM}g -Djava.io.tmpdir=$TEMP \
	-jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R ${REF_FASTA} \
	--dbsnp $PUBLIC_VARIANT_DIR/dbsnp_137.hg19.vcf \
        --intervals $INTERVAL_BED \
        --interval_padding $INTERVAL_PAD \
	-o ${STEP_TARGET} \
	--input_file:$BAM_HEALTHY_TAG ${BAM_HEALTHY_LIST} \
	--input_file:$BAM_SLE_ONLY_TAG ${BAM_SLE_ONLY_LIST} \
	--input_file:$BAM_ESRD_TAG ${BAM_ESRD_LIST} \
	--input_file:$BAM_HAPMAP_TAG ${BAM_HAPMAP_LIST} \
        "
	echo $CMD
	$CMD 
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC>/dev/null; fi
	touch ${STEP_TARGET}.done
	echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
    fi
    
    # rest needs more coding
    #exit 0 > /dev/null

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
            java -Xmx${JVM_RAM}g \
            \
            -Djava.io.tmpdir=$TEMP \
	    -jar ${SNPEFF_DIR}/SnpSift.jar \
	    annMem \
	    -id ${PUBLIC_VARIANT_DIR}/dbsnp_137.hg19.vcf \
	    ${GATK_RAW_VCF} \
        "
	echo $CMD "> ${STEP_TARGET}"
	$CMD > ${STEP_TARGET}
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC > /dev/null; fi
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
	    RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC > /dev/null; fi
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
	    RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC > /dev/null; fi
	    touch ${STEP_TARGET}.done
	    echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
	fi

	# Remove intermediate files
	echo;echo `date`"  TS      START removing intermediate files  ${OUT_NAME}"
	CMD="echo DEBUG ECHO rm -f ${TEMP_FILE_LIST}"
	echo $CMD
	$CMD
	RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR	${OUT_NAME}"; exit $RC > /dev/null; fi
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
	    RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC > /dev/null; fi
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
	    RC=$?; if [ $RC != 0 ]; then echo;echo `date`"	TS	ERROR ${STEP_NAME}	${OUT_NAME}"; exit $RC > /dev/null; fi
	    touch ${STEP_TARGET}.done
	    echo;echo `date`"	TS	DONE ${STEP_NAME}	${OUT_NAME}"; 
	fi
    fi
    
    exit 0 > /dev/null
fi
