#!/bin/bash
######################################################################
# 
# produce several TFAM files (one per case/control) from a .VCF file
# by 
# ANNOVAR-2013-09-11/convert2annovar.sh -format vcf4old -include -comment 
#
# combined Cohorts:
# ASW                HEALTHY
# CEU                HEALTHY
# YRI                HEALTHY
# WUSTL.HEALTHY      HEALTHY
# WUSTL.SLE          SLE
# WUSTL.ESRD         ESRD+SLE
# Lupus              ESRD+SLE
# RA                 RA
######################################################################

# ARGS
IN=$1  # vcf4old.inc.com/foo.avinput
OUT_BASE=`basename $1 .avinput`
PHENO_CODING=$2
if [ ! -e "$IN" ]; then 
    echo "ERROR: can't read input .vcf file: $IN"; 
    exit 1 2>&1 > /dev/null
fi
if [ ! -e "$PHENO_CODING" ]; then 
    echo "ERROR: can't read input coding2group2sample file: $PHENO_CODING"; 
    echo "ERROR: first col=sample name"
    echo "ERROR: second col=cohort/phenotype name"
    exit 1 2>&1 > /dev/null
fi
echo "PHENO_CODING=$PHENO_CODING"; 
if [[ -z "$PHENO_CODING" || ! -e "$PHENO_CODING" ]]; then echo "!!MISSING!!"; exit 1; fi
PHENO_DEF=`grep "^##PHENO_DEF=" $PHENO_CODING | perl -pe 's/.*=([^ \t#]+).*/$1/;'`
echo "PHENO_DEF=$PHENO_DEF"; if [ ! -e $PHENO_DEF ]; then echo "!!MISSING!!"; exit 1; fi
if [[ -z "$PHENO_DEF" || ! -e "$PHENO_DEF" ]]; then echo "!!MISSING!!"; exit 1; fi
SAMPLE_DEF=`grep "^##SAMPLE_DEF=" $PHENO_DEF | perl -pe 's/.*=([^ \t#]+).*/$1/;'`
echo "SAMPLE_DEF=$SAMPLE_DEF"; 
if [[ -z "$SAMPLE_DEF" || ! -e "$SAMPLE_DEF" ]]; then echo "!!MISSING!!"; exit 1; fi

# get sample list
if [ `grep -c "^#CHROM" $IN` != 1 ]; then
    echo "ERROR: .vcf file must have exactly one #CHROM header with sample names: $IN"
    exit 1 2>&1 > /dev/null
fi
    
SAMPLE_LIST=`grep "#CHROM" $IN | cut -f 10-`
SAMPLE_NUM=`echo $SAMPLE_LIST | wc -w`
echo "Found $SAMPLE_NUM samples"

# for each phenotype definition
COMPARE_LIST=`sort $PHENO_CODING | awk '(/^#/){next;}($1!=LVAL){print $1; LVAL=$1}'`
for phenotype in $COMPARE_LIST; do

    echo "PHENOTYPE: $phenotype"
    OUT=$OUT_BASE.$phenotype.tfam
    echo "    Writing $OUT"

    # one line per sample, 6 columns: 
    # 1.   Family ID
    # 2.   Individual ID
    # 3.   Paternal ID
    # 4.   Maternal ID
    # 5.   Sex (1=male; 2=female; other=unknown)
    # 6.   Phenotype
    #    Phenotype can be quantitative or affection status
    #            affection = { -9 missing,  0 missing,  1 unaffected,  2 affected }
    # constants
    STATUS_MISSING=0
    STATUS_UNAFFECTED=1
    STATUS_AFFECTED=2
    # counters
    STATUS_0=0
    STATUS_1=0
    STATUS_2=0

    for sample in $SAMPLE_LIST; do 
	INDIVIDUAL_ID=$sample
	FAMILY_ID=$sample # every individual is a founder
	PATERNAL_ID=0 # no family info available
	MATERNAL_ID=0 # no family info available
	SEX_ID=2 # female
	PHENOTYPE_STATUS=0 # missing

	# look up sample group/cohort
	SAMPLE_GROUP=`egrep "^$sample\s" $SAMPLE_DEF|cut -f 2`
	if [ -z "$SAMPLE_GROUP" ]; then echo "ERROR: MISSING $sample NOT FOUND IN $SAMPLE_DEF"; exit 1; fi
	# look up group/cohort phenotype
	SAMPLE_PHENO=`egrep "^$SAMPLE_GROUP\s" $PHENO_DEF|cut -f 2`
	if [ -z "$SAMPLE_PHENO" ]; then echo "ERROR: MISSING $SAMPLE_GROUP NOT FOUND IN $PHENO_DEF"; exit 1; fi
	# look up phenotype coding for this comparison & this group and phenotype
	PHENO_CODE=`awk '($1=="'$phenotype'"&&($2=="'$SAMPLE_GROUP'"||$3=="'$SAMPLE_PHENO'")){print $4}' $PHENO_CODING`
	if [ -n "$PHENO_CODE" ]; then
	    PHENOTYPE_STATUS=$PHENO_CODE
	fi

	# count phenotype status
	eval 'STATUS_'$PHENOTYPE_STATUS'=$(($STATUS_'$PHENOTYPE_STATUS' + 1))'

	# debug
	echo "     $PHENOTYPE_STATUS $sample"

	# output
	TAB="	"
	echo -n "$FAMILY_ID$TAB" >> $OUT
	echo -n "$INDIVIDUAL_ID$TAB" >> $OUT
	echo -n "$PATERNAL_ID$TAB" >> $OUT
	echo -n "$MATERNAL_ID$TAB" >> $OUT
	echo -n "$SEX_ID$TAB" >> $OUT
	echo "$PHENOTYPE_STATUS$TAB" >> $OUT
    done;
    echo "    Mapped sample->phenotype_status: "
    echo "        0 mising     $STATUS_0"
    echo "        1 unaffected $STATUS_1"
    echo "        2 affected   $STATUS_2"
    echo "          total      "$(($STATUS_0 + $STATUS_1 + $STATUS_2)) 
done

    