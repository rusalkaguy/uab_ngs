#!/bin/bash

if [ "$1" == "-q" ]; then QUIET=TRUE; shift 1; fi

if [ -z "$1" ]; then 
    awk 'NR%4==2{bases+=length($0)}END{print bases}'
    RC=$?; exit $RC
else
    for FQ in $*; do
	if [ -z "$QUIET" ]; then echo -n $FQ"	"; fi
	if [[ "$FQ" == *.gz ]]; then
	    awk 'NR%4==2{bases+=length($0)}END{print bases}' <(zcat $FQ)
	    RC=$?; if [ $RC != 0 ]; then exit $RC; fi
	else
	    awk 'NR%4==2{bases+=length($0)}END{print bases}' $FQ
	    RC=$?; if [ $RC != 0 ]; then exit $RC; fi
	fi
    done
fi
