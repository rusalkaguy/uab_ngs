#!/bin/bash
echo -n "START "; date
IN_SNPSIFT_CC_TXT=$1
IN_SNPSIFT_CC_VCF=`dirname $IN_SNPSIFT_CC_TXT`/`basename $IN_SNPSIFT_CC_TXT .txt`.vcf
IN_TABLE_ANNOVAR=$2


OUT_DATA_USORT=`basename $IN_SNPSIFT_CC_TXT .txt`.annovar.data.unsorted
OUT_DATA=`basename $IN_SNPSIFT_CC_TXT .txt`.annovar.data.txt
OUT_FIELDS=`basename $IN_SNPSIFT_CC_TXT .txt`.annovar.fields.txt

if [ ! -e "$IN_SNPSIFT_CC_TXT" ]; then echo "ERROR: can't find ARG1: IN_SNPSIFT_CC_TXT: $IN_SNPSIFT_CC_TXT"; exit 1; fi
if [ ! -e "$IN_SNPSIFT_CC_VCF" ]; then echo "ERROR: can't find       IN_SNPSIFT_CC_VCF: $IN_SNPSIFT_CC_VCF"; exit 1; fi
if [ ! -e "$IN_TABLE_ANNOVAR" ];  then echo "ERROR: can't find ARG2: IN_TABLE_ANNOVAR:  $IN_TABLE_ANNOVAR"; exit 1; fi

IN_SCCT_VARS=`grep -c "^c" $IN_SNPSIFT_CC_TXT`
echo "$IN_SCCT_VARS variants in $IN_SNPSIFT_CC_TXT"
IN_SCCV_VARS=`grep -c "^c" $IN_SNPSIFT_CC_VCF`
echo "$IN_SCCV_VARS variants in $IN_SNPSIFT_CC_VCF"
IN_TA_VARS=`grep -c "^c" $IN_TABLE_ANNOVAR`
echo "$IN_TA_VARS variants in $IN_TABLE_ANNOVAR"

if [ $IN_SCCT_VARS != $IN_TA_VARS ]; then
    echo "ERROR: variant counts do not match!"
    exit 1
fi

echo "pasting data"
paste \
    <(grep -v "^##" $IN_SNPSIFT_CC_TXT| cut -f 1-) \
    <(grep -v "^#" $IN_TABLE_ANNOVAR | cut -f 1-) \
    > $OUT_DATA_USORT

(head -n 1 $OUT_DATA_USORT; \
    grep -v "^#" $OUT_DATA_USORT \
	| ./uab_ngs/linux_plus/sort_chr_pos.sh \
) > $OUT_DATA

echo ""
echo "OUT_DATA + header " `wc -l $OUT_DATA`

echo "extracting header"
echo "SRC	ID	Number	Type	Description" > $OUT_FIELDS
egrep "^##(FORMAT|INFO)=" $IN_SNPSIFT_CC_VCF \
    | perl -pe 's/^.*<ID=(.*),Number=(.*),Type=(.*),Description="(.*)">.*$/join("\t",$1,$2,$3,$4)/e;' \
    >> $OUT_FIELDS
echo "OUT_FIELDS " `wc -l $OUT_FIELDS`
 
    

