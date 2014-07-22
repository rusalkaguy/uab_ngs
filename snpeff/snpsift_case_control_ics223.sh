#!/bin/bash
echo -n "START "; date
# ******** inputs ***********
VCF_IN=$PWD/hc_pad200_hg19.vcf
NAME=HC200 
ANNOVAR_IN=vcf4old.inc.com/hc_pad200_hg19.reduce.step4.varlist
if [ -n "$1" ]; then VCF_IN="$1"; fi
if [ -n "$2" ]; then NAME="$2"; fi
if [ -n "$3" ]; then ANNOVAR_IN="$3"; fi
echo "VCF_IN=$VCF_IN"; if [ ! -e $VCF_IN ]; then echo "!!MISSING!!"; exit 1; fi
echo "NAME=$NAME"
echo "ANNOVAR_IN=$ANNOVAR_IN"; if [ ! -e $ANNOVAR_IN ]; then echo "!!MISSING!!"; exit 1; fi

# ********* derived names ***************
IN_BASE=`basename $VCF_IN .avinput`
OUT_VCF=`basename $VCF_IN .avinput`.caseControl.vcf
OUT_TXT=`basename $VCF_IN .avinput`.caseControl.txt
echo "OUT_TXT=$OUT_TXT"

echo "***************************************************************"
echo 'to clean dir: rm *.*.vcf *.tfam */*.vcf '$OUT_TXT $OUT_VCF
#rm *.*.vcf *.tfam */*.vcf $OUT_TXT $OUT_VCF


echo "***************************************************************"
echo "# extract subset VCF from ANNOVAR REDUCEd list"
ANNOVAR_VCF=${ANNOVAR_IN}.vcf
if [[ ! -e $ANNOVAR_VCF || "$ANNOVAR_IN" -nt "$ANNOVAR_VCF" || "$0" -nt "$ANNOVAR_VCF" ]]; then 
    echo "grep ^# $VCF_IN > $ANNOVAR_VCF"
    grep "^#" $VCF_IN > $ANNOVAR_VCF
    RC=$?; if [ $RC != 0 ]; then echo "ERROR: RC=$RC"; exit $RC; fi
    echo -n "# "; wc -l $ANNOVAR_VCF
    echo "grep -v ^# $ANNOVAR_IN | cut -f 6- >> $ANNOVAR_VCF"
    grep -v "^#" $ANNOVAR_IN | cut -f 6- >> $ANNOVAR_VCF
    RC=$?; if [ $RC != 0 ]; then echo "ERROR: RC=$RC"; exit $RC; fi
    echo -n "# "; wc -l $ANNOVAR_VCF
else
    echo "SKIP"
fi
echo `grep -vc "^#" $ANNOVAR_VCF`" variants in $ANNOVAR_VCF"

echo "***************************************************************"
echo "generate .tfam files from VCF $IN"
INPUT=$VCF_IN
SCRIPT=~/uab_ngs/snpeff/vcf2tped_ics223.sh
OUTPUT=$IN_BASE.sleVnorm.tfam
if [[ ! -e "$OUTPUT" || "$INPUT" -nt "$OUTPUT" || "$SCRIPT" -nt "$OUTPUT" ]]; then 
    $SCRIPT $INPUT
    RC=$?; if [ $RC != 0 ]; then echo "ERROR: RC=$RC"; exit $RC; fi
else
    echo SKIP
fi

echo "***************************************************************"
echo "Generate phenotype counts in VCF"
PHE_IN=$ANNOVAR_VCF
# for each phenotype definition
for phenotype in sleVnorm esrdVnorm esrdVsle esrdVsleNorm; do

    echo "PHENOTYPE: $phenotype"
    PHE_OUT=$IN_BASE.$phenotype.vcf

    #**************************************************
    echo "caseControl first phenotype: $phenotype"
    INPUT=$PHE_IN
    INPUT2=$IN_BASE.${phenotype}.tfam
    OUTPUT=$PHE_OUT
    if [[ ! -e "$OUTPUT" || "$INPUT" -nt "$OUTPUT" || "$INPUT2" -nt "$OUTPUT" ]]; then 
	CMD="java -jar /share/apps/ngs-ccts/snpEff_3_3/SnpSift.jar caseControl \
	    -tfam $INPUT2 \
	    -name _${NAME}_$phenotype \
	    $INPUT"
	echo "$CMD > $PHE_OUT"
	$CMD > $PHE_OUT
	RC=$?
	if [ $RC != 0 ]; then echo "ERROR: RC=$RC"; exit $RC; fi
    else
	echo SKIP
    fi

    # use as input for next step
    PHE_IN=$PHE_OUT
done

echo "***************************************************************"
echo "linking $OUT_VCF"
CMD="ln -fs $PHE_OUT $OUT_VCF"
$CMD
RC=$?; if [ $RC != 0 ]; then echo "ERROR: RC=$RC"; fi
echo ""



echo "***************************************************************"
echo "Excel Extract" #**************************************************

echo "Creating $OUT_TXT"
INPUT=$OUT_VCF
SCRIPT=$0
OUTPUT=$OUT_TXT
if [[ ! -e "$OUTPUT" || "$INPUT" -nt "$OUTPUT" || "$SCRIPT" -nt "$OUTPUT" ]]; then 

    # extract field list
    echo "Extract field list"
    COL_LIST=`grep -v "^#" $OUT_VCF | grep Controls | head -n 1 \
              | cut -f 8 \
              | tr ";" "\n" \
              | egrep "^(Cases|Control|CC_)" \
              | cut -d = -f 1 \
              | uniq`
    echo "    "`echo $COL_LIST | wc -w`" fields ("`echo $COL_LIST | tr "	" "\n"| cut -d _ -f 1 | sort | uniq | perl -pe 's/\n//;'`")"
    echo "    $COL_LIST"
    echo "Subscript field list"
    for f in $COL_LIST; do
	if [[ "$f" == Cases* || "$f" == Controls* ]]; then
	    COL_LIST_SUB="$COL_LIST_SUB ${f}[0] ${f}[1] ${f}[2]"
	    MOD_LIST="$MOD_LIST $f"
	else
	    COL_LIST_SUB="$COL_LIST_SUB $f"
	fi
    done
    echo "run SnpSift.jar/extractFields"
    CMD="java -jar /share/apps/ngs-ccts/snpEff_3_3/SnpSift.jar extractFields \
	$OUT_VCF \
	CHROM POS REF ALT QUAL FILTER AD DP GQ GT PL AC AF AN BaseQRankSum ClippingRankSum DB DP DS FS HaplotypeScore InbreedingCoeff MLEAC MLEAF MQ MQ0 MQRankSum QD ReadPosRankSum $COL_LIST_SUB"
    echo "$CMD | [column name remap] > $OUT_TXT"
    $CMD | \
	perl -pe 's/((Cases|Controls)_\S+)\[0\]/$1_HOM/g;s/((Cases|Controls)_\S+)\[1\]/$1_HET/g;s/((Cases|Controls)_\S+)\[2\]/$1_TOT/g;' \
	> $OUT_TXT
    RC=$?; if [ $RC != 0 ]; then echo "ERROR: RC=$RC"; fi
else
    echo "SKIP"
fi

echo -n "COMPLETE "; date
