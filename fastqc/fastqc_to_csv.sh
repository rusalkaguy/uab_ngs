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

    TAB="	"
    EMODULE_STATS=`grep -n ">>END_MODULE" Ross-1_S1_L001_R1_001_fastqc/fastqc_data.txt  | head -n 1 | tail -n 1  | cut -d : -f 1`
    EMODULE_QUAL=`grep -n ">>END_MODULE" Ross-1_S1_L001_R1_001_fastqc/fastqc_data.txt  | head -n 2 | tail -n 1 | cut -d : -f 1`
    HEAD_LINES=$((${EMODULE_QUAL}-1))
    TAIL_LINES=$((${EMODULE_QUAL}-${EMODULE_STATS}-2-1))
    #
    # build header
    #
    if [ -z $HEADER_DONE ]; then 
	echo -n "Filename"
	echo -n "$TAB"
	echo -n "Total Sequences"
	echo -n "$TAB"
	echo -n "Sequence Length"
	echo -n "$TAB"
	echo -n "GC"
	echo -n "$TAB"
	grep "^>>" $TARG | grep -v END_MODULE | cut -f 1 | cut -c 3- | perl -pe 's/\n/\t/;' # section headers

        # fetch base range
        BASE_HEADER="\""`head -n ${HEAD_LINES} $TARG | tail -n ${TAIL_LINES} | cut -f 1 | perl -pe 's/\n/"\t"/;'`"\""
        #BASE_HEADER="'"`head -n ${HEAD_LINES} $TARG | tail -n ${TAIL_LINES} | cut -f 1 | perl -pe 's/\n/''\t/;'`
	
	# columnar sections
	for c in Mean Median "Lower Quartile" "Upper Quartile"; do
	    echo -n "$c"
	    echo -n "$TAB"
	    echo -n "$BASE_HEADER"
	done
	
	# only emit the header once
	HEADER_DONE=1
	echo ""
    fi
    
    # 
    # emit the actual values
    #
    echo -n `grep "^Filename" $TARG | cut -f 2`"$TAB"
    echo -n `grep "^Total Sequences" $TARG | cut -f 2`"$TAB"
    echo -n `grep "^Sequence length" $TARG | cut -f 2`"$TAB"
    echo -n `grep '^%GC' $TARG | cut -f 2`"$TAB"
    grep "^>>" $TARG | grep -v END_MODULE | cut -f 2 | perl -pe 's/\n/\t/;' # section results (pass/fail)


    # mean
    echo -n "Mean"; echo -n "$TAB"
    head -${HEAD_LINES} $TARG | tail -${TAIL_LINES} | cut -f 2 | perl -pe 's/\n/\t/;'
    # median
    echo -n "Median"; echo -n "$TAB"
    head -${HEAD_LINES} $TARG | tail -${TAIL_LINES} | cut -f 3 | perl -pe 's/\n/\t/;' 
    # lower Quartile
    echo -n "Lower Quartile"; echo -n "$TAB"
    head -${HEAD_LINES} $TARG | tail -${TAIL_LINES} | cut -f 4 | perl -pe 's/\n/\t/;' 
    # upper Quartile
    echo -n "Upper Quartile"; echo -n "$TAB"
    head -${HEAD_LINES} $TARG | tail -${TAIL_LINES} | cut -f 5 | perl -pe 's/\n/\t/;' 

    echo ""

    # finished
    #echo "	done" 

done
