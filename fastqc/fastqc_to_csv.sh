#!/bin/bash
#
# Extract data from the fastqc_data.txt files in a galaxy (or non-galaxy) run of FASTQC
# so that data from multiple samples can be easily combined in Excel. 
#
if [[ "$1" == "-o" && -n "$2" ]]; then
    echo "OUTFILE not supported [yet]"
    exit 1
    OUTFILE="$2"
    shift
    shift
fi

if [[ "$1" == "" ]] ; then
    IN_FILES="-"
else
    IN_FILES="$*"
fi


for TARG in $IN_FILES; do

    if [ "$TARG" == "-" ]; then
	TARG=`mktemp`
	cat - > $TARG
    fi

    FILE=`grep "^Filename" $TARG | cut -f 2`
    #echo "FASTQ: $FILE"
    #
    # extract base quality
    #

    # filename
    echo -n "$FILE	" 
    # seq count
    grep "^Total Sequences" $TARG | perl -pe 's/\n//;' 
    grep "^Sequence length" $TARG | perl -pe 's/\n//;' 
    SEQ_LEN=`grep "^Sequence length" $TARG | cut -f 2`
    grep "^>>" $TARG | grep -v END_MODULE | perl -pe 's/\n/\t/;s/>>//g;' 

    HEAD_LINES=`echo "13+$SEQ_LEN"|bc -q`
    TAIL_LINES=`echo "1+$SEQ_LEN"|bc -q`
    # mean
    head -${HEAD_LINES} $TARG | tail -${TAIL_LINES} | cut -f 2 | perl -pe 's/\n/\t/;'
    # median
    head -${HEAD_LINES} $TARG | tail -${TAIL_LINES} | cut -f 3 | perl -pe 's/\n/\t/;' 
    # lower Quartile
    head -${HEAD_LINES} $TARG | tail -${TAIL_LINES} | cut -f 4 | perl -pe 's/\n/\t/;' 
    # upper Quartile
    head -${HEAD_LINES} $TARG | tail -${TAIL_LINES} | cut -f 5 | perl -pe 's/\n/\t/;' 

    echo ""

    # finished
    #echo "	done" 

done
