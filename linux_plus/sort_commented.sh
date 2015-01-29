#!/bin/bash
#
# sort a file, without touching any header lines that start with #
#
#

# get leading flags 
SORT_FLAGS=
while [[ -n "$1" && "$1" != "-" && ("$1" == -* || ! -e "$1") ]]; do 
    SORT_FLAGS+="$1 "; shift; echo "SORT_FLAGS=$SORT_FLAGS"; 
done

# the rest are source files
if [ -z "$1" ]; then
    SRCS=-
else
    SRCS="$*"
fi
# stdin is a special case, since we scan the input files more than once
SAFE_SRCS=
STDIN=
for SRC in $SRCS; do
    if [ "$SRC" == "-" ]; then
	if [ -z "$STDIN" ]; then
	    STDIN=`mktemp`; trap 'rm -rf $STDIN' EXIT
	    cat - > $STDIN
	fi
	SAFE_SRCS+="$STDIN "
    else
	SAFE_SRCS+="$SRC "
    fi
done

grep -h "^#"    $SAFE_SRCS;                    RC=$?; if [ $RC != 0 ]; then exit $RC; fi
grep -h -v "^#" $SAFE_SRCS | sort $SORT_FLAGS; RC=$?; if [ $RC != 0 ]; then exit $RC; fi
exit 0
