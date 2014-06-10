#!/bin/bash
# 
# run full variant filter/annotate protocol on a VCF
# 1. Annovar/Filter
# 3. Annovar/Annotate
# 2. SNPSift (case/control count/pvalue)
#
#
VCF_IN=$1
ABBREV_IN=$2
echo `# CMD: basename $0`" $*"

DIR_BASENAME=vcf4old.inc.com/`basename $VCF_IN .vcf`
AVINPUT=${DIR_BASENAME}.avinput
FILTERED_AVINPUT=${DIR_BASENAME}.reduce.step2.varlist  # step2 = # of steps in doit.variant_reductions.sh - need a ln -s for final to be less fragile
FULL_MULTIANNO_TXT=${DIR_BASENAME}.full.hg19_multianno.txt  
MULTIANNO_TXT=${FILTERED_AVINPUT}.hg19_multianno.txt  

CASECONTROL_VCF=`basename $VCF_IN .vcf`.caseControl.vcf
CASECONTROL_TXT=`basename $VCF_IN .vcf`.caseControl.txt
OUT_FIELDS=`basename $VCF_IN .vcf`.caseControl.annovar.fields.txt
OUT_DATA=`basename $VCF_IN .vcf`.caseControl.annovar.data.txt
OUT_DATA_COLCLEAN=`basename $VCF_IN .vcf`.caseControl.annovar.data.cols.txt
OUT_DATA_HFILT=`basename $VCF_IN .vcf`.caseControl.annovar.data.cols.filt.txt
OUT_DATA_HFILT_STATS=`basename $VCF_IN .vcf`.caseControl.annovar.data.cols.filt.stats
OUT_DATA_AAF05=`basename $VCF_IN .vcf`.caseControl.annovar.data.cols.filt.aaf05.txt
OUT_DATA_AAF05_STATS=`basename $OUT_DATA_AAF05 .txt`.stats
OUT_DATA_AAF03=`basename $VCF_IN .vcf`.caseControl.annovar.data.cols.filt.aaf03.txt
OUT_DATA_AAF03_STATS=`basename $OUT_DATA_AAF03 .txt`.stats
OUT_DATA_AAF01=`basename $VCF_IN .vcf`.caseControl.annovar.data.cols.filt.aaf01.txt
OUT_DATA_AAF01_STATS=`basename $OUT_DATA_AAF01 .txt`.stats
OUT_DATA_NOVEL=`basename $VCF_IN .vcf`.caseControl.annovar.data.cols.filt.aaf00.txt
OUT_DATA_NOVEL_STATS=`basename $OUT_DATA_NOVEL .txt`.stats

if [[ -z "$VCF_IN" || ! -e "$VCF_IN" || -z $ABBREV_IN ]]; then
    echo "USAGE: "`basename $0`" input.VCF abbrev"
    exit 1
fi

echo "#"
echo "# ANNOVAR: convert VCF to AVINPUT"
echo "#"
INPUT=$VCF_IN
OUTPUT=$AVINPUT
SCRIPT=~/uab_ngs/annovar/convert2annovar.vcf4old.sh
if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; then echo "SKIP"; else
    date
    $SCRIPT $VCF_IN $AVINPUT
    RC=$?; date
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"
echo ""

# turn off filtering
FILTERED_AVINPUT=$AVINPUT
#echo "#"
#echo "# ANNOVAR: filter out uninteresting variants"
#echo "#"
#INPUT=$AVINPUT
#OUTPUT=$FILTERED_AVINPUT
#SCRIPT=~/uab_ngs/annovar/doit.variants_reduction.sh
#if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; thenecho "SKIP"; else
#    date
#    $SCRIPT $AVINPUT
#    RC=$?; date
#    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
#fi
#echo "# OUTPUT=$OUTPUT"
#echo ""

echo "#"
echo "# ANNOVAR: annotate remaining variants"
echo "#"
INPUT=$FILTERED_AVINPUT
OUTPUT=$MULTIANNO_TXT
SCRIPT=~/uab_ngs/annovar/doit.table_annovar.sh
if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; then echo "SKIP"; else
    date
    $SCRIPT $INPUT $OUTPUT
    RC=$?; date
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"

echo "#"
echo "# ANNOVAR: annotate ALL variants"
echo "#"
INPUT=$AVINPUT
OUTPUT=$FULL_MULTIANNO_TXT
SCRIPT=~/uab_ngs/annovar/doit.table_annovar.sh
if [[ 1 == 1 || -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; then echo "SKIP"; else
    date
    $SCRIPT $INPUT $OUTPUT
    RC=$?; date
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"

echo "#"
echo "# SnpSift - compute variant counts/cohort & p-values"
echo "#"
INPUT=$AVINPUT
OUTPUT=$CASECONTROL_VCF
OUTPUT2=$CASECONTROL_TXT
SCRIPT=~/uab_ngs/snpeff/snpsift_case_control_ics223.sh
if [ -e "$OUTPUT" ]; then echo "SKIP"; else
    date
    $SCRIPT $INPUT $ABBREV_IN $FILTERED_AVINPUT
    RC=$?; date
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"
echo "# OUTPUT2=$OUTPUT2"
echo ""

echo "#"
echo "# Merge outputs"
echo "#"
INPUT=$CASECONTROL_TXT 
INPUT2=$MULTIANNO_TXT
OUTPUT=$OUT_DATA
OUTPUT2=$OUT_FIELDS
SCRIPT=~/uab_ngs/annovar/merge.snpeff-case-control.table_annovar.sh
if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT"  && "$OUTPUT" -nt "$INPUT2" ]]; then echo "SKIP"; else
    date; 
    $SCRIPT $INPUT $INPUT2
    RC=$?; date
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"
echo "# OUTPUT2=$OUTPUT2"
echo ""

echo ""
echo "#"
echo "# ANNOVAR: cleanup compound columns"
echo "#"
INPUT=$OUT_DATA
OUTPUT=$OUT_DATA_COLCLEAN
SCRIPT=~/uab_ngs/annovar/annovar_multianno_split_columns.pl
if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; then echo "SKIP"; else
    date
    cat $INPUT | $SCRIPT > $OUTPUT
    RC=$?; date
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"
echo ""

echo ""
echo "#"
echo "# ANNOVAR: hand-filter annotated variants"
echo "#"
INPUT=$OUT_DATA_COLCLEAN
OUTPUT=$OUT_DATA_HFILT
OUTPUT2=$OUT_DATA_HFILT_STATS
SCRIPT=~/uab_ngs/annovar/annovar_multianno_filter.pl
if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; then echo "SKIP"; else
    date
    cat $INPUT | $SCRIPT > $OUTPUT 2> $OUTPUT2
    RC=$?; date
    cat $OUTPUT2
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"
echo ""

echo ""
echo "#"
echo "# ANNOVAR: hand-filter 5% variants"
echo "#"
INPUT=$OUT_DATA_HFILT
OUTPUT=$OUT_DATA_AAF05
OUTPUT2=$OUT_DATA_AAF05_STATS
SCRIPT=~/uab_ngs/annovar/filter_named_column.pl
if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; then echo "SKIP"; else
    date
    cat $INPUT | $SCRIPT max_aaf "<=" 0.05 > $OUTPUT 2> $OUTPUT2
    RC=$?; date
    cat $OUTPUT2
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"
echo ""

echo ""
echo "#"
echo "# ANNOVAR: hand-filter 3% variants"
echo "#"
INPUT=$OUT_DATA_HFILT
OUTPUT=$OUT_DATA_AAF03
OUTPUT2=$OUT_DATA_AAF03_STATS
SCRIPT=~/uab_ngs/annovar/filter_named_column.pl
if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; then echo "SKIP"; else
    date
    cat $INPUT | $SCRIPT max_aaf "<=" 0.03 > $OUTPUT 2> $OUTPUT2
    RC=$?; date
    cat $OUTPUT2
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"
echo ""

echo ""
echo "#"
echo "# ANNOVAR: hand-filter 1% variants"
echo "#"
INPUT=$OUT_DATA_HFILT
OUTPUT=$OUT_DATA_AAF01
OUTPUT2=$OUT_DATA_AAF01_STATS
SCRIPT=~/uab_ngs/annovar/filter_named_column.pl
if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; then echo "SKIP"; else
    date
    cat $INPUT | $SCRIPT max_aaf "<=" 0.01 > $OUTPUT 2> $OUTPUT2
    RC=$?; date
    cat $OUTPUT2
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"
echo ""

echo ""
echo "#"
echo "# ANNOVAR: hand-filter novel variants"
echo "#"
INPUT=$OUT_DATA_HFILT
OUTPUT=$OUT_DATA_NOVEL
OUTPUT2=$OUT_DATA_NOVEL_STATS
SCRIPT=~/uab_ngs/annovar/filter_named_column.pl
if [[ -e "$OUTPUT" && "$OUTPUT" -nt "$SCRIPT" && "$OUTPUT" -nt "$INPUT" ]]; then echo "SKIP"; else
    date
    cat $INPUT | $SCRIPT is_novel "eq" 1 > $OUTPUT 2> $OUTPUT2
    RC=$?; date
    if [ $RC != 0 ]; then echo "FAILED: RC=$RC"; exit $RC; fi
fi
echo "# OUTPUT=$OUTPUT"
echo ""
echo "#"
echo "# FINAL OUTPUT"
echo "#"
grep -vc "^#" \
    $VCF_IN $FILTERED_AVINPUT \
    $CASECONTROL_TXT \
    $MULTIANNO_TXT $MULTIANNO_COLCLEAN_TXT $MULTIANNO_HFILT_TXT \
    $OUT_FIELDS \
    $OUT_DATA $OUT_DATA_COLCLEAN $OUT_DATA_HFILT \
    $OUT_DATA_AAF05 $OUT_DATA_AAF03 $OUT_DATA_AAF01 \
    $OUT_DATA_NOVEL \
    






