#!/bin/bash
#
# GATK SelectVariants
# 
# create 1 VCF per group from a combined VCF using SAMPLE->GROUP list file
#
# http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html
#
VCF_IN=$1            ; if [ ! -e "$VCF_IN" ]; then           echo "ERROR: file not found VCF_IN=$VCF_IN"; fi
echo VCF_IN $VCF_IN
SAMPLE_GROUP_MAP=$2  ; if [ ! -e "$SAMPLE_GROUP_MAP" ]; then echo "ERROR: file not found SAMPLE_GROUP_MAP=$SAMPLE_GROUP_MAP"; fi
echo SAMPLE_GROUP_MAP=$SAMPLE_GROUP_MAP 
shift 2; GROUP_LIST=$*    
echo GROUP_LIST=$GROUP_LIST 

if [[ -z "$VCF_IN" || -z "$SAMPLE_GROUP_MAP" ]]; then 
    echo "SYNTAX: $0 VCF_IN SAMPLE_GROUP_MAP [group1 [group2...]]"
    exit 1;
fi

# GATK SETUP
module load java/jre1.7.0_51
GATK=/share/apps/ngs-ccts/GenomeAnalysisTK/GenomeAnalysisTK-3.0-0/GenomeAnalysisTK.jar 
IGV_TOOLS=/share/apps/ngs-ccts/IGVTools/IGVTools_2.3.32/igvtools.jar
REF=`grep "^##reference=file:" $VCF_IN | cut -d / -f 3-`
echo REF=$REF
if [ ! -e "$REF" ]; then
    echo "ERROR: REF .fasta doesn't exist: $REF"
    exit 1
fi
#REF=/scratch/share/public_datasets/ngs/databases/gatk_bundle/2.5/b37/human_g1k_v37.fasta

# threading
if [ -z "$NSLOTS" ]; then NSLOTS=1; fi 
if [ -z "$RAM" ]; then RAM=6; fi 

# output dir
mkdir -p JOBS SPLIT

# check VCF is index
if [[ "$VCF_IN" == *.vcf && ! -e "${VCF_IN}.idx" ]]; then
    echo INDEXING $VCF_IN
    java -Xmx${RAM}g -Xms${RAM}g -jar $IGV_TOOLS \
	index \
	$VCF_IN
    RC=$?
    if [ "$RC" != 0 ];then 
	echo "ERROR: GATK returned $RC"
	exit $RC
    fi
fi
    
# get list of groups
if [ -z "$GROUP_LIST" ]; then 
    echo "Loading group list from $SAMPLE_GROUP_MAP"
    GROUP_LIST=`awk '(/^#/){next;}{print $2}' $SAMPLE_GROUP_MAP | sort | uniq | xargs echo`
    echo "Found "`echo $GROUP_LIST | wc -w`" groups: $GROUP_LIST"
fi

# for each group
for GROUP in $GROUP_LIST; do 
    echo GROUP=$GROUP

    VCF_OUT="SPLIT/"`basename $VCF_IN .vcf`".$GROUP.vcf"

    # get list of samples
    SAMPLE_NAMES=`awk --assign "GROUP=$GROUP" '(/^#/){next;}(GROUP==$2){print $1}' $SAMPLE_GROUP_MAP | sort | uniq | xargs echo`
    #echo SAMPLE_NAMES=$SAMPLE_NAMES
    # build sample name flags
    ARG_LIST=
    for SNAME in $SAMPLE_NAMES; do 
	ARG_LIST="$ARG_LIST -sn $SNAME"
    done

    # run GATK SelectVariants
    # --excludeNonVariants (-env) : remove variants not observed in selected samples
    echo "$GROUP has " `echo $SAMPLE_NAMES | wc -w` "samples"
    java -Xmx${RAM}g -Xms${RAM}g -jar $GATK \
	-nt $NSLOTS \
	-R $REF \
	-T SelectVariants \
	--variant $VCF_IN \
	-o $VCF_OUT \
	--excludeNonVariants \
	$ARG_LIST \
	> JOBS/$GROUP.split.so \
	2> JOBS/$GROUP.split.se 
    RC=$?
    grep . JOBS/$GROUP.split.s?
    if [ "$RC" != 0 ];then 
	echo "ERROR: GATK returned $RC"
	exit $RC
    fi
	
done

jobs 
wait

exit 0

