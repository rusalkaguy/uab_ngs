#!/bin/bash
#
# GATK CombineVariants 
# http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html
#

#if [[ -z $NSLOTS ]]; then
#	$NSLOTS=1
#fi

TARGETGENE="ITGAM"

module load java/jre1.7.0_51
GATK=/share/apps/ngs-ccts/GenomeAnalysisTK/GenomeAnalysisTK-3.0-0/GenomeAnalysisTK.jar 
REF=/scratch/share/public_datasets/ngs/databases/gatk_bundle/2.5/b37/human_g1k_v37.fasta
MERGE_OPT=PRIORITIZE # UNIQUIFY
# get list of cohorts
GROUP_LIST=`ls *.vcfBeta*vcf | cut -f 1 -d . | uniq`
echo GROUP_LIST=$GROUP_LIST

ALL_GROUP=`echo $GROUP_LIST | perl -pe 's/ /./g;'`
echo ALL_GROUP=$ALL_GROUP
SUFFIX="vcfBeta-ALL-ASM.$TARGETGENE.vcf"
mkdir -p MERGED JOBS

# for each group
for GROUP in $GROUP_LIST; do 

    GROUP_OUT="MERGED/$GROUP.$SUFFIX"

    # get list of samples
    SAMPLES=`ls $GROUP.vcfBeta-GS*.vcf $GROUP.vcfBeta-NA*.vcf 2>/dev/null`

    # build variant list
    VAR_LIST=
    PRIO_LIST=
    INDEX=0
    for S in $SAMPLES; do 
	INDEX=$(($INDEX+1))
	PRIO=v$INDEX
	VAR_LIST="$VAR_LIST --variant:$PRIO $S"
	if [ -n "$PRIO_LIST" ]; then PRIO_LIST="$PRIO_LIST,$PRIO"; else PRIO_LIST=$PRIO; fi
    done

    # run GATK merger
    echo "$GROUP has " `ls -1 $SAMPLES | wc -l` "samples"
    java -Xmx20g -Xms20g -jar $GATK \
	-nt 1 \
	-R $REF \
	-T CombineVariants \
	-o $GROUP_OUT \
	-genotypeMergeOptions $MERGE_OPT  \
	-priority $PRIO_LIST \
	$VAR_LIST \
	> JOBS/$GROUP.merge.so \
	2> JOBS/$GROUP.merge.se 
done

jobs 
wait

#
# merge groups
#

# build variant list
VAR_LIST=
    PRIO_LIST=
    INDEX=0
for S in $GROUP_LIST; do 
	INDEX=$(($INDEX+1))
	PRIO=v$INDEX
	VAR_LIST="$VAR_LIST --variant:$PRIO MERGED/$S.$SUFFIX"
	if [ -n "$PRIO_LIST" ]; then PRIO_LIST="$PRIO_LIST,$PRIO"; else PRIO_LIST=$PRIO; fi
done
ALL_GROUP_OUT=MERGED/$ALL_GROUP.$SUFFIX

# run GATK merger
echo "$ALL_GROUP has " `echo $GROUP_LIST | wc -w` "samples"
java -Xmx20g -Xms20g -jar $GATK \
	-nt 1 \
	-R $REF \
	-T CombineVariants \
	-o $ALL_GROUP_OUT \
	-genotypeMergeOptions $MERGE_OPT \
	-priority $PRIO_LIST \
	$VAR_LIST \
	> JOBS/$ALL_GROUP.merge.so \
	2> JOBS/$ALL_GROUP.merge.se 

jobs
wait
