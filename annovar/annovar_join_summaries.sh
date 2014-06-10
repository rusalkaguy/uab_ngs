#!/bin/bash
#
# join a set of summary files in order
# 
# file schema (TSV): chr loc gene count
#
# SYNTAX: $0 file1 file2 [file3...]
#

# get list of genes
MASTER=`mktemp`
cut -f 1,2,3 $* \
    |sort -k3,3 -k1.4,1.10n -k2,2n  \
    | uniq \
    | awk '{print $0,"0"}' \
    > $MASTER


LEFT_OUT="1.1 1.2 0 1.4"
RIGHT_OUT="2.4"
NEXT_OUT=4
# while there are files, join to the master
while [ -n "$1" ]; do
    OUT_FILE=`mktemp`
    #echo "Join [$NEXT_OUT] $1"
    join -e 0 -o $LEFT_OUT $RIGHT_OUT -a1 -a2 -1 3 -2 3 \
	<(sort -k3,3 $MASTER) \
	<(sort -k3,3 $1) \
	> $OUT_FILE

    # move to next file
    mv $OUT_FILE $MASTER
    NEXT_OUT=$(($NEXT_OUT +1))
    LEFT_OUT="$LEFT_OUT 1.$NEXT_OUT"
    shift
done

# output result -w/o the 0 count from the original master
sort -k1.4,1.10n -k2,2n $MASTER \
    | cut -d " " --output-delimiter="	" -f 1-3,5-

# clean up
rm $MASTER


