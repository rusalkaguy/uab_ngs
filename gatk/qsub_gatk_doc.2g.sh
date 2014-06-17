#!/bin/bash
#
# Compute Depth Of Coverge from BAM file 
# http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_coverage_DepthOfCoverage.html
# 
# libraries to load
. /etc/profile.d/modules.sh          # enable module loading
. ~/uab_ngs/uab_ngs_functions_v1.sh  # load shared run_cmd() & run_step()
# *** QSUB FLAGS ***
#$ -S /bin/bash
#$ -cwd
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -j y # merge stderr into stdout
# email at:  Begining, End, Abort, Suspend
#$ -m beas  
#
# *** DRMAA resources for GATK ****
# DepthOfCoverage IS SINGLE THREADED (-nt > 1 crashes)
# -pe smp 4 -l vf=1.9G -l h_vmem=1.9G
#$ -l vf=2G -l h_vmem=2G
if [ -z "$NSLOTS" ]; then NSLOTS=1; fi
JAVA_XMX=1500m
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
#*** END QSUB ****

TASK_NAME=gatk_doc
CMD_LINE_PARAM_LIST="PROJECT_DIR BAM GENE_LIST"
DERIVED_VAR_LIST="CMD_LINE HOSTNAME CT NSLOTS JAVA_XMX"

export CT="-ct 10 -ct 50 -ct 100 -ct 200 -ct 500 -ct 1000 -ct 5000 -ct 100000 -ct 500000"

# modules
module load java/1.7.0_25
export PATH=/share/apps/java/1.7.0_25/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/1.7.0_25/lib:$LD_LIBRARY_PATH
module load ngs-ccts/samtools-0.1.19
export PATH=/share/apps/ngs-ccts/samtools-0.1.19:$PATH
export PATH=/share/apps/ngs-ccts/samtools-0.1.19/bcftools:$PATH
export PATH=/share/apps/ngs-ccts/samtools-0.1.19/misc:$PATH
export GATK_DIR="/share/apps/ngs-ccts/GenomeAnalysisTK-2.8-1-g932cd3a"
export PUBLIC_VARIANT_DIR="/scratch/share/public_datasets/ngs/databases/gatk_bundle/2.3"
export PICARD_DIR="/share/apps/ngs-ccts/picard-tools-1.87"

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
	    echo "**** -PE OVERRIDE=$QSUB_DRMAA **** ) " 
	    shift 1
	    shift 1
	    shift 1
	    continue
	fi
	# unknown flag
	echo "ERROR: unknown option: $1"
	qsub_exit 1
    done
    # cmd-line params
    for myvar in $CMD_LINE_PARAM_LIST ; do
	eval $myvar="$1"
	export $myvar
	echo -n "$myvar	:"; eval echo \$$myvar
	if [ -z "$1" ] ; then
	    echo "$myvar	: MISSING"
	    SCRIPT_NAME=`basename $0`
	    echo ""
	    echo "ERROR: $SCRIPT_NAME [-debug] $CMD_LINE_PARAM_LIST"
	    echo ""
	    echo ""
	    qsub_exit 1
	fi
	shift
    done

    # Prog & INDEX versions
    export TEMP=/scratch/user/${USER}/tmp
    mkdir -p ${TEMP}

    # what version of find to use.
    FIND_CMD="find"; which lfs 2>/dev/null >/dev/null; if [ $? == 0 ]; then FIND_CMD="lfs find"; fi

    # create job output dir
    GATK_OUT_DIR=.
    JOB_DIR=${GATK_OUT_DIR}/jobs
    OJOB_DIR=${JOB_DIR}/old
    mkdir -p ${JOB_DIR} ${OJOB_DIR}
    #mv ${JOB_DIR}/${TASK_NAME}-${OUT_NAME}* ${OJOB_DIR} > /dev/null 2>&1

    # check for multiple targets (BAMs)
    export BAM_LIST=$BAM
    if [ -d "$BAM" ]; then
	echo "SCANNING DIRECTORY $BAM"
	export BAM_LIST=`$FIND_CMD $BAM -name "*sort.bam"`
	echo BAM_LIST=$BAM_LIST
	JOB_ID=
	for bam in $BAM_LIST; do
	    echo EXE $0 $PROJECT_DIR "$bam" "$GENE_LIST"
	    $0 $PROJECT_DIR "$bam" "$GENE_LIST"
	    RC=$?
	    if [ $RC != 0 ]; then qsub_exit $RC; fi
	done
	qsub_exit $RC
    fi

    # qsub
    export OUT_FNAME=`dirname $BAM`/`basename $BAM .bam`.doc.txt
    export OUT_NAME=`echo $OUT_FNAME | sed -e 's|/|-|g'`
    if [ -z "$JOB_ID" ]; then
	echo -n "${OUT_NAME}:QSUB:"
	pushd ${GATK_OUT_DIR} > /dev/null
	QSUB_NAME=${TASK_NAME}-${OUT_NAME}
	qsub -terse \
	    -N $QSUB_NAME \
	    -M $USER@uab.edu \
	    $0
	popd > /dev/null
	if [ $? != 0 ]; then echo "ERROR: bad return code from QSUB"; qsub_exit 1; fi
	qsub_exit 0
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
    #printenv
    echo "-- /printenv -- "

    echo "I'm a qsub slave: "

    #
    # find matching fasta
    #
    BFILE=`basename $BAM`
    if [[ $BAM = *_alignVsAssembly_sort.bam ]]; then
	REF_FASTA=`dirname $BAM`/`basename $BAM _alignVsAssembly_sort.bam`_fixed_assembly_w60.fa
	GENE_LIST=
    fi
    if [[ $BAM = *_alignVsRef_sort.bam ]]; then 
	REF_FASTA=~/ics/consults/ross/ics63/ref_genomes/NC_006273.2_HCMV_Merlin.fa
	GENE_LIST="~/ics/consults/ross/ics63/ref_genomes/NC_006273.2_HCMV_Merlin.gb"
    fi
    GENE_LIST_ARG=
    if [ -n "$GENE_LIST" ]; then GENE_LIST_ARG="-geneList $GENE_LIST"; fi

    # 
    # PICARD/SAMTOOLS index FASTA
    #
    STEP_TARGET=`dirname ${REF_FASTA}`/`basename ${REF_FASTA} .fa`.dict
    run_step $OUT_NAME $STEP_TARGET "PICARD_DICT" - \
	java -Xmx${JAVA_XMX} -jar $PICARD_DIR/CreateSequenceDictionary.jar \
	R= ${REF_FASTA} \
	O= ${STEP_TARGET}

    STEP_TARGET=${REF_FASTA}.fai
    run_step $OUT_NAME $STEP_TARGET "SAMTOOLS_INDEX" - \
	samtools faidx ${REF_FASTA}
    
    #
    # GATK DepthOfCoverage
    #
    STEP_TARGET=`dirname $BAM`/`basename $BAM .bam`.doc.txt
    STEP_NAME="GATK RealignerTargetCreator"
    run_step $OUT_NAME $STEP_TARGET "GATK_DOC" - \
	java -Xmx$JAVA_XMX -Djava.io.tmpdir=$TEMP \
	-jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	-nt $NSLOTS \
	-T DepthOfCoverage \
        -R ${REF_FASTA} \
	-o $STEP_TARGET \
        -I $BAM \
        $GENE_LIST_ARG \
	$CT 
    fi

    
    qsub_exit 0
fi
