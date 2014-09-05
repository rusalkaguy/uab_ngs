#!/bin/bash
#
# summarize variant count per gene, with gene posisiton
#
SRC=$1
SAMPLE_NAME=1
if [ -n "$2" ]; then SAMPLE_NAME=$2; fi

# select chr,start & Gene.refGene
REF_GENE_COL=`head -n 1 $SRC | transpose | grep -n "Gene.refGene" | cut -d : -f 1`
cut -f 1,2,$REF_GENE_COL $SRC \
    | sort -k3,3 -k1.4,1.10n -k2,2n \
    | awk -v "SAMPLE_NAME=$SAMPLE_NAME" \
    ' \
        BEGIN{OFS="\t";chr="";loc=0;ct=0;gene=""} \
        {sub(/[(;].*/,"",$3);} \
        (gene!=$3 && gene != "" && "#CHROM"==chr){print chr,loc,gene,SAMPLE_NAME;ct=0;} \
        (gene!=$3 && gene != "" && "#CHROM"!=chr){print chr,loc,gene,ct;ct=0;} \
        {chr=$1;loc=$2;gene=$3;ct++;}' \
    | sort -k1.4,1.10n -k2,2n


