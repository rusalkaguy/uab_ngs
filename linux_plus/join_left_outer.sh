#!/bin/bash
#
# outer_join key_columns out_columns [flags]
#
# 
TAB="	"; KEY_SEP=:
KEY_COLUMNS=$1
OUT_COLUMNS=$2
FILE1=$3
FILE2=$4
OUTFILE=$5

if [[ -z "$KEY_COLUMNS" || -z "$OUT_COLUMNS" || -z "$FILE1" || -z "$FILE2" || -z "$OUTFILE" ]]; then
    echo "ERROR: "`basename $0`" key_columns_csv out_columns_csv src1 src2 outfile"
    exit 1
fi
if [ ! -e "$FILE1" ]; then echo "ERORR: ${FILE1}: doesn't exist"; exit 1; fi
wc -l $FILE1
if [ ! -e "$FILE2" ]; then echo "ERORR: ${FILE2}: doesn't exist"; exit 1; fi
wc -l $FILE2

#
# process input file to build keys, subset columns
#
FILE1_KEYED=.`basename $FILE1 .txt`.keyed
FILE2_KEYED=.`basename $FILE2 .txt`.keyed

paste -d "$TAB" \
    <(cut -f $KEY_COLUMNS -d "$TAB" --output-delimiter "$KEY_SEP" $FILE1) \
    <(cut -f $OUT_COLUMNS -d "$TAB" --output-delimiter "$TAB"     $FILE1) \
    | `dirname $0`/sort_commented.sh \
    >  $FILE1_KEYED
wc -l $FILE1_KEYED

paste -d "$TAB" \
    <(cut -f $KEY_COLUMNS -d "$TAB" --output-delimiter "$KEY_SEP" $FILE2) \
    <(cut -f $OUT_COLUMNS -d "$TAB" --output-delimiter "$TAB"     $FILE2) \
    | `dirname $0`/sort_commented.sh \
    >  $FILE2_KEYED
wc -l $FILE2_KEYED

#
# master key list
#
KEY_FILE=.key_list
cut -f 1 $FILE1_KEYED $FILE2_KEYED \
    | `dirname $0`/sort_commented.sh \
    | uniq \
    > $KEY_FILE
wc -l $KEY_FILE

# 
# build output list for join
#
JOIN_FORMAT1=
JOIN_FORMAT2=
INDEX=1
for C in `grep -v "^#" $FILE1_KEYED | head -n 1 | perl -pe 's/[^\t]*/1/g;'`; do
    JOIN_FORMAT1+=,1.$INDEX
    JOIN_FORMAT2+=,2.$INDEX
    INDEX=$(($INDEX+1))
done
echo "JOIN_FORMAT1=$JOIN_FORMAT1"
echo "JOIN_FORMAT2=$JOIN_FORMAT2"

#
# do the join
#
JOINFILE1=.join_file1
join -t "$TAB" -a 1 -a 2  -e "" -o 0$JOIN_FORMAT1$JOIN_FORMAT2 $FILE1_KEYED $FILE2_KEYED \
    > $JOINFILE1
wc -l $JOINFILE1

JOINFILE=$JOINFILE1
#
# expand the keys
#
paste -d "$TAB" \
    <(cut -f 1 $JOINFILE | cut -d "$KEY_SEP" --output-delimiter "$TAB" -f 1- ) \
    <(cut -f 2- $JOINFILE ) \
    > $OUTFILE
wc -l $OUTFILE
