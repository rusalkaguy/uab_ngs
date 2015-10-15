#!/bin/bash
#
# CORRECTLY sort a file that starts with col1=chr*, col2=pos
# to match GATK canonical sort order
# http://www.broadinstitute.org/gatk/guide/article?id=1204
#

# get leading flags 
PREFIX=chr
SORT_FLAGS=
while [[ -n "$1" && "$1" != "-" && ("$1" == -* || ! -e "$1") ]]; do 
    if [ "$1" == "-chr" ]; then 
	PREFIX=
	shift
	#echo "# assume no chr prefix"
    else
	SORT_FLAGS+="$1 "; shift; echo "SORT_FLAGS=$SORT_FLAGS"; 
    fi
done

# the rest are source files
# all normal file names
if [ -z "$1" ]; then
    SRCS=-
else
    SRCS="$*"
fi


# stdin is a special case, since we scan the input files more than once
SAFE_SRCS=
STDIN=
# prefix handling
if [ -n "$PREFIX" ]; then 
    PEND=$((`echo $PREFIX | wc -c` + 1))
else
    PEND=1
fi
#echo "PEND=$PEND"
for SRC in $SRCS; do
    if [ "$SRC" == "-" ]; then
	if [ -z "$STDIN" ]; then
	    STDIN=`mktemp`
	    cat - > $STDIN
	fi
	SAFE_SRCS+="$STDIN "
    else
	SAFE_SRCS+="$SRC "
    fi
done

# do the actual work
grep -h "^#"            $SAFE_SRCS ; RC=$?; if [ $RC != 0 ]; then exit $RC; fi
grep -h "^$PREFIX[0-9]"     $SAFE_SRCS | sort -k1.$PEND,1.100n -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then exit $RC; fi
grep -h "^$PREFIX[XY]"      $SAFE_SRCS | sort -k1.$PEND,1.5    -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then exit $RC; fi
grep -h "^$PREFIX[M]"       $SAFE_SRCS | sort -k1.$PEND,1.5    -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then exit $RC; fi
grep -h "^$PREFIX[^0-9XYM]" $SAFE_SRCS | sort -k1.$PEND,1.100  -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then exit $RC; fi
exit 0
