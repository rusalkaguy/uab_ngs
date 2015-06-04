#!/bin/bash
#
# CORRECTLY sort a file that starts with col1=chr*, col2=pos
# to match GATK canonical sort order
# http://www.broadinstitute.org/gatk/guide/article?id=1204
#

# get leading flags 
SORT_FLAGS=
while [[ -n "$1" && "$1" != "-" && ("$1" == -* || ! -e "$1") ]]; do 
    SORT_FLAGS+="$1 "; shift; echo "SORT_FLAGS=$SORT_FLAGS"; 
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

# do the actual work
grep -h "^#"            	$SAFE_SRCS  # headers
grep -Ewh "^chr[M]"       	$SAFE_SRCS | sort -k1.4,1.5    -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then echo ERROR; exit $RC; fi
grep -Ewh "^chr[0-9]"		$SAFE_SRCS | sort -k1.4,1.4n   -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then echo ERROR; exit $RC; fi
grep -Ewh "^chr[0-9][0-9]"	$SAFE_SRCS | sort -k1.4,1.5n   -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then echo ERROR; exit $RC; fi
grep -Ewh "^chr[XY]"      	$SAFE_SRCS | sort -k1.4,1.5    -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then echo ERROR; exit $RC; fi
#grep -Ewh "^chr[0-9]+[^0-9]+"  $SAFE_SRCS | sort -k1.4,1.100  -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then echo ERROR; exit $RC; fi
grep -Ewh "^chr[0-9]_.*"	$SAFE_SRCS | sort -k1.4,1.4n   -k1.5,1.100  -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then echo ERROR; exit $RC; fi
grep -Ewh "^chr[0-9][0-9]_.*"	$SAFE_SRCS | sort -k1.4,1.5n   -k1.6,1.100  -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then echo ERROR; exit $RC; fi
grep -Ewh "^chrUn_g.*"  	$SAFE_SRCS | sort -k1.8,1.100  -k2,2n $SORT_FLAGS;   RC=$?; if [ $RC != 0 ]; then echo ERROR; exit $RC; fi
exit 0
