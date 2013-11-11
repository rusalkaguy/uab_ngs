#!/bin/bash
######################################################################
# 
# recover a VCF file from an avinput file created
# by 
# ANNOVAR-2013-09-11/convert2annovar.sh -format vcf4old -include -comment 
#  with or without -allallele (the uniq fixes that)
######################################################################

# ARGS
IN=$1
OUT=`basename $1 .avinput`.rev.vcf
if [ ! -e "$IN" ]; then 
    echo "ERROR: can't read input avinput file: $IN"; 
    exit 1 2>&1 > /dev/null
fi

# headers - remove trailing stuff
grep "^##" $IN | cut -f 1 > $OUT
# "#CHROM" header - don't remove trailing stuff
grep "^#[^#]" $IN  >> $OUT
# variants - strip first 5 columns
grep -v "^#" $IN |cut -f 6- | uniq >> $OUT

echo "wrote "`grep -cv "^#" $OUT`" variants to $OUT"

