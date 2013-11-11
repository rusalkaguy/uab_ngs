#!/bin/bash
######################################################################
# 
# produce several TFAM files (one per case/control) from a .VCF file
# by 
# ANNOVAR-2013-09-11/convert2annovar.sh -format vcf4old -include -comment 
#  with or without -allallele (the uniq fixes that)
######################################################################

# ARGS
IN=$1
OUT_BASE=`basename $1 .avinput`
if [ ! -e "$IN" ]; then 
    echo "ERROR: can't read input .vcf file: $IN"; 
    exit 1 2>&1 > /dev/null
fi

# get sample list
if [ `grep -c "^#CHROM" $IN` != 1 ]; then
    echo "ERROR: .vcf file must have exactly one #CHROM header with sample names: $IN"
    exit 1 2>&1 > /dev/null
fi
    
SAMPLE_LIST=`grep "#CHROM" $IN | cut -f 10-`
SAMPLE_NUM=`echo $SAMPLE_LIST | wc -w`
echo "Found $SAMPLE_NUM samples"

# for each phenotype definition
for phenotype in sleVnorm esrdVnorm; do

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
    STATUS_0=0
    STATUS_1=0
    STATUS_2=0

    for sample in $SAMPLE_LIST; do 
	INDIVIDUAL_ID=$sample
	FAMILY_ID=$sample # every individual is a founder
	PATERNAL_ID=0 # no family info available
	MATERNAL_ID=0 # no family info available
	SEX_ID=2 # female
	
	# phenotype parsed from sample name
	PHENOTYPE_STATUS=1 # normal/unaffected

	if [[ "$sample" == *HapMap* || "$sample" == *NEG* || "$sample" == *Yung* ]]; then
	    # HapMap samples are outside our dataset
	    PHENOTYPE_STATUS=0
	    #echo "MISSING $sample"
	else
	    # which phenotype are we looking at
	    if [ "$phenotype" == "sleVnorm" ]; then
		if [[ "$sample" == *SLE* || "$sample" == *NON_NEPHRITIS* ]]; then
		    PHENOTYPE_STATUS=2
		fi
	    fi	    
	    if [ "$phenotype" == "esrdVnorm" ]; then
		if [[ "$sample" == *ESRD* ]]; then
		    PHENOTYPE_STATUS=2
		fi
	    fi	    
	fi
	 
	# count phenotype status
	eval 'STATUS_'$PHENOTYPE_STATUS'=$(($STATUS_'$PHENOTYPE_STATUS' + 1))'

	# debug
	#echo "     $PHENOTYPE_STATUS $sample"

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

    