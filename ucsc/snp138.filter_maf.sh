#!/bin/sh
IN=$1
if [[ -z "$IN" || ! -e "$IN" ]]; then "IN=$IN ERROR:MISING"; exit 1; fi
AF=$2
if [[ -z "$AF" ]]; then AF=20; fi
echo AF=$AF
OUT=`basename $IN .gz`.aaf$AF
CAT=cat
CAT_OUT=cat
if [[ "$IN" == *.gz ]]; then CAT=zcat; CAT_OUT="bgzip -c"; OUT="$OUT.gz"; fi

# column name to index
COL=`$CAT $IN | head | grep "^#chrom" | ~/uab_ngs/linux_plus/transpose | grep -n "alleleFreqs" | cut -d : -f 1`
echo COL=$COL

module load ngs-ccts/tabix/0.2.6

# filter only records with a MAJOR allele frequency of 80% (thus a aaf of 20%)
( $CAT $IN | head -n 1000 | grep "#" | grep -v "^#chr" ; \
    $CAT $IN | \
    awk 'BEGIN{FS="\t";OFS="\t"}("#chrom"==$1){print $0,"MajAlleleFeq"}(/^#/){next;}{split($'$COL',mafs,","); n=asort(mafs);if(mafs[n]<(1.0-0.'$AF')){print $0,mafs[n]}}' \
    ) \
    | $CAT_OUT \
    > $OUT

#    perl -pe '@lines = split("\t",chomp($_));$_=$lines['$COL'-1];'


#  chomp($_);$_ = $_ . "=>" . join("|",reverse(sort(split(m/,/,$_)))) . "\n";' | head
