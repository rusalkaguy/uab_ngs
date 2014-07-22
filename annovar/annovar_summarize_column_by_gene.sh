#!/bin/bash
#
# summarize variant count per gene, with gene posisiton
#
SAMPLE_NAME=1
if [ -n "$1" ]; then SAMPLE_NAME=$1; fi
COLUMN="AC"
if [ -n "$2" ]; then COLUMN=$2; fi
T=`mktemp`
cat - > $T
COLNUM=`grep "#CHR" $T | ~/uab_ngs/linux_plus/transpose | grep -n $COLUMN | cut -d : -f 1`
if [ -z "$COLNUM" ]; then 
    echo "ERROR: column not found: $COLUMN in "`grep "^#CHR" $T`
    exit 1
fi
echo "[$SAMPLE_NAME] COLNUM=$COLNUM  ($COLUMN)"
# select chr,start & Gene.refGene
paste <(cut -f 1,2,76 $T) <(cut -f $COLNUM $T) \
    | sort -k3,3 -k1.4,1.10n -k2,2n \
    | awk -v "SAMPLE_NAME=$SAMPLE_NAME" \
    ' \
        BEGIN{OFS="\t";chr="";loc=0;ct=0;gene=""} \
        {sub(/[(;].*/,"",$3);} \
        (gene!=$3 && gene != "" && "#CHROM"==chr){print chr,loc,gene,SAMPLE_NAME;ct=0;} \
        (gene!=$3 && gene != "" && "#CHROM"!=chr){print chr,loc,gene,ct;ct=0;} \
        {chr=$1;loc=$2;gene=$3;ct+=$4;}' \
    | sort -k1.4,1.10n -k2,2n


