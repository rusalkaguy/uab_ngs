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

AVINPUT=vcf4old.inc.com/`basename $VCF_IN .vcf`.avinput
FILTERED_AVINPUT=vcf4old.inc.com/`basename $VCF_IN .vcf`.reduce.step4.varlist
MULTIANNO_TXT=vcf4old.inc.com/`basename $VCF_IN .vcf`.reduce.step4.varlist.hg19_multianno.txt
CASECONTROL_VCF=`basename $VCF_IN .vcf`.caseControl.vcf
CASECONTROL_TXT=`basename $VCF_IN .vcf`.caseControl.txt
OUT_DATA=`basename $VCF_IN .vcf`.caseControl.annovar.data.txt
OUT_FIELDS=`basename $VCF_IN .vcf`.caseControl.annovar.fields.txt

if [[ -z "$VCF_IN" || ! -e "$VCF_IN" || -z $ABBREV_IN ]]; then
    echo "USAGE: "`basename $0`" input.VCF abbrev"
    exit 1
fi

echo "#"
echo "# ANNOVAR: convert VCF to AVINPUT"
echo "#"
OUTPUT=$AVINPUT
if [ -e "$OUTPUT" ]; then echo "SKIP"; else
    ~/uab_ngs/annovar/convert2annovar.vcf4old.sh $VCF_IN $AVINPUT
fi
echo "# OUTPUT=$OUTPUT"
echo ""


echo "#"
echo "# ANNOVAR: filter out uninteresting variants"
echo "#"
OUTPUT=$FILTERED_AVINPUT
if [ -e "$OUTPUT" ]; then echo "SKIP"; else
    ~/uab_ngs/annovar/doit.variants_reduction.sh $AVINPUT
fi
echo "# OUTPUT=$OUTPUT"
echo ""

echo "#"
echo "# ANNOVAR: annotate remaining variants"
echo "#"
OUTPUT=$MULTIANNO_TXT
if [ -e "$OUTPUT" ]; then echo "SKIP"; else
    ~/uab_ngs/annovar/doit.table_annovar.sh $FILTERED_AVINPUT
fi
echo "# OUTPUT=$OUTPUT"
echo ""

echo "#"
echo "# SnpSift - compute variant counts/cohort & p-values"
echo "#"
OUTPUT=$CASECONTROL_VCF
OUTPUT2=$CASECONTROL_TXT
if [ -e "$OUTPUT" ]; then echo "SKIP"; else
    ~/uab_ngs/snpeff/snpsift_case_control_ics223.sh  \
     $AVINPUT $ABBREV_IN $FILTERED_AVINPUT
fi
echo "# OUTPUT=$OUTPUT"
echo "# OUTPUT2=$OUTPUT2"
echo ""

echo "#"
echo "# Merge outputs"
echo "#"
OUTPUT=$OUT_DATA
OUTPUT2=$OUT_FIELDS
if [ -e "$OUTPUT" ]; then echo "SKIP"; else
    ~/uab_ngs/annovar/merge.snpeff-case-control.table_annovar.sh \
    $CASECONTROL_TXT \
    $MULTIANNO_TXT 
fi
echo "# OUTPUT=$OUTPUT"
echo "# OUTPUT2=$OUTPUT2"
echo ""
echo "#"
echo "# FINAL OUTPUT"
echo "#"
grep -vc "^#" $VCF_IN $FILTERED_AVINPUT $CASECONTROL_TXT $MULTIANNO_TXT $OUT_DATA $OUT_FIELDS






