#!/bin/bash
#
# summarize variant count per gene, with gene posisiton
#
SAMPLE_NAME=1
if [ -n "$1" ]; then SAMPLE_NAME=$1; fi
OPERATION='if($4!="0"){ct++}' # count rows
if [ -n "$2" ]; then
    if [[ "$2" != "count" && "$2" != "sum" ]]; then 
	echo "ERROR: operation must be 'count' or 'sum', not '$2'"
	echo "SYNTAX: $0 abbreviation operation column_name"
	exit 1
    fi
    if [ "$2" == "sum" ]; then 
	OPERATION='ct+=$4'
    fi
fi    
COLUMN="AC"
if [ -n "$3" ]; then COLUMN=$3; fi
T=`mktemp`
cat - > $T

# look up column numbers
COLNUM=`head -n 1000 $T | grep "#CHR" | ./uab_ngs/linux_plus/transpose | grep -wn "$COLUMN" | cut -d : -f 1`
if [ -z "$COLNUM" ]; then 
    echo "ERROR: column not found: $COLUMN in "`grep "^#CHR" $T`
    exit 1
fi
GENENAME_CN=`head -n 1000 $T | grep "#CHR" | ./uab_ngs/linux_plus/transpose | grep -wn "Gene.refGene" | cut -d : -f 1`
if [ -z "$GENENAME_CN" ]; then 
    echo "ERROR: column not found: Gene.refGene in "`grep "^#CHR" $T`
    exit 1
fi
# select chr,start & Gene.refGene
# (strip off qualified gene names before the sort)
paste <(cut -f 1,2,$GENENAME_CN $T) <(cut -f $COLNUM $T) \
    | awk '{sub(/[(;].*/,"",$3);print $0}' \
    | sort -k3,3 -k1.4,1.10n -k2,2n \
    | awk -v "SAMPLE_NAME=$SAMPLE_NAME" \
    ' \
        BEGIN{OFS="\t";chr="";loc=0;ct=0;gene=""} \
        {sub(/[(;].*/,"",$3);} \
        (gene!=$3 && gene != "" && "#CHROM"==chr){print chr,loc,gene,SAMPLE_NAME;ct=0;} \
        (gene!=$3 && gene != "" && "#CHROM"!=chr){print chr,loc,gene,ct;ct=0;} \
        {chr=$1;loc=$2;gene=$3;'"$OPERATION"';} \
        END{print chr,loc,gene,ct}' \
    | sort -k1.4,1.10n -k2,2n


