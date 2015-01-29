#!/bin/bash
#
# JOIN multiple files on a single key extracting a single data column
#
# test case
#  for x in {1..9}; do echo "geneA  ${x}00" > p${x}.txt; done
#  join_multiple.sh p?.txt
#  cat final.txt | column -t
#   ID     p1.txt  p2.txt  p3.txt  p4.txt  p5.txt  p6.txt  p7.txt  p8.txt  p9.txt
#   geneA  100     200     300     400     500     600     700     800     900
#  
# arguments: list of files to join
FILELIST=$*
echo FILELIST=$FILELIST
# handle temp file cleanup 
# http://stackoverflow.com/questions/687014/removing-created-temp-files-in-unexpected-bash-exit
export MY_TEMP_DIR=`mktemp -d` # establish dir for all our tmp files
trap 'rm -rf $MY_TEMP_DIR' EXIT  # auto-delete all temp files

# output destination
FINAL=final.txt
# key and data column indicies for input files
COLS=1,2

# start with no previous file
PREV=
# iterate over all files
for CUR in $FILELIST; do 
    # create auto-delete temp file for intermediate result
    OUT=`TMPDIR="" mktemp -p $MY_TEMP_DIR` 
    
    # 
    # on first time through, just strip columns
    # we don't (yet) have a PREV to merge
    #
    #echo "PREV=$PREV"
    if [ -z "$PREV" ]; then 
	cut -f $COLS $CUR | sort > $OUT
	PREV=$OUT
	echo "first file $CUR"
	continue
    fi

    #
    # merge our running intersect list with CUR data file
    #
    join -j 1 \
	$PREV \
	<(cut -f $COLS $CUR | sort) \
	> $OUT
    echo "merged "`wc -l $CUR`
    #echo "merged $PREV and $CUR -> $OUT" # too many details
    echo "	current intersection "`wc -l $OUT`
    #grep . $OUT 

    # shift current intersection to be the PREV file in next merge
    rm $PREV # cleanup temp files we're done with
    PREV=$OUT
done

# compute headers
HEADERS="ID $*"

# pre-pend headers to result
(echo $HEADERS; cat $OUT) > $FINAL
rm $OUT
echo "output in $FINAL"
wc -l $FINAL