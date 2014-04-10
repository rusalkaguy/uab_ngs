#!/bin/bash
#
# shared subroutine used by the qsub_* scripts in these subdirectories
#

# parse -xxx flags and declared parameters
parse_params () {
    parse_params_shift_count=0
    
    # capture the original CMD_LINE
    if [ -z "$CMD_LINE" ]; then  export CMD_LINE="$0 $*"; fi

    # check for debug or no qsub - so it will work outside the cluster
    QSUB=`which qsub 2>/dev/null`
    if [[ -z "$QSUB" ]]; then 
	export JOB_ID=run_now
	echo "**** NO QSUB FOUND ****" 
    fi

    # parse FLAGS on cmd-line
    while [[ "$1" == -* ]]; do 
	echo "PARSE FLAG: $1"
	if [[ "-debug" == "$1" || "-inline" == "$1" ]]; then
	    shift 1; parse_params_shift_count=$(( $parse_params_shift_count + 1 ))
	    export JOB_ID=run_now
	    echo "**** NO QSUB [NSLOTS=$NSLOTS] ****" 
	    continue
	fi
	# unknown flag
	echo "ERROR: unknown option: $1"
	qsub_exit 1
    done
    # parse PARAMS on cmd-line 
    for myvar in $CMD_LINE_PARAM_LIST ; do
	eval $myvar=$1
	export $myvar
	echo -n "Z: $myvar	:"; eval echo \$$myvar
	if [ -z "$1" ] ; then
	    echo "$myvar	: MISSING"
	    SCRIPT_NAME=`basename $0`
	    echo ""
	    echo "ERROR: $SCRIPT_NAME [-debug|-inline] ${CMD_LINE_PARAM_LIST}"
	    echo ""
	    echo ""
	    qsub_exit 1
	fi
	shift; parse_params_shift_count=$(( $parse_params_shift_count + 1 ))

    done
}

# put positional params in vars specified in CMD_LINE_PARAM_LIST
cmd_line_params () {
    for myvar in $CMD_LINE_PARAM_LIST ; do
	eval $myvar=$1
	export $myvar
	echo -n "${myvar}="; eval echo \$$myvar
	shift
    done
    return 0
}

print_cmd_line_params () {
    for myvar in $CMD_LINE_PARAM_LIST; do
	echo -n "$myvar="; eval echo \$$myvar
    done
    return 0
}

print_derived_params () {
    for myvar in $DERIVED_VAR_LIST; do
	echo -n "$myvar="; eval echo \$$myvar
    done
    return 0
}

# UTIL function

# exit from a qsub script adds unprintable chars to output
# send those chars to /dev/null
qsub_exit () { # ARGS: rc
    RC=$1
    exit $RC 2>&1 > /dev/null
    return 0
}

run_cmd () { # ARGS: stdout_dest_or_- cmd [args]
    TMPSTD=$1; shift  # redirect for stdout
    TMPERR=`mktemp`
    if [[ -z "$TMPSTD" || "-" == "$TMPSTD" ]]; then  
	# stdout to stdout
	echo "# CMD: $*" 
	eval $* 2>$TMPERR
    else
	# stdout to file
	echo "# CMD: $* 1>$TMPSTD" 
	eval $* 2>$TMPERR 1>$TMPSTD
    fi
    RC=$?
    if [ $RC != 0 ]; then 
	echo "ERROR: $*"
	echo "ERROR: "`cat $TMPERR`
	qsub_exit 1
    fi
    return $RC
}

#====================================================================== 
# UTIL FUNCTIONS
#====================================================================== 
run_step () { # ARGS: sample_name target_file step_name cmd_out cmd [cmd_args]
    # args
    RS_SAMPLE_NAME=$1; shift 1
    RS_TARGET=$1; shift 1
    RS_NAME=$1; shift 1
    RS_STDOUT=$1; shift 1       # redirect for stdout

    RS_TARGET_DONE=`dirname $RS_TARGET`/.done.`basename $RS_TARGET`

    # log start/skip
    if [[ ( -e "${RS_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${RS_TARGET_DONE}" ) ]]; then
	echo;echo `date`"	TS	SKIP ${RS_NAME}	${RS_SAMPLE_NAME}"; 
    else
	echo;echo `date`"	TS	START ${RS_NAME}	${RS_SAMPLE_NAME}"; 
	
	# run the cmd
	#RS_TMPERR=`mktemp`
	if [[ -z "$RS_STDOUT" || "-" == "$RS_STDOUT" ]]; then  
	    # stdout to stdout
	    echo "# CMD: $*" 
	    eval $* #2>$RS_TMPERR
	else
    	    # stdout to file
	    echo "# CMD: $* 1>$RS_STDOUT" 
	    eval $* 1>$RS_STDOUT #2>$RS_TMPERR 
	fi

	#  handle errors
	RC=$?
	if [ $RC != 0 ]; then 
	    echo "ERROR: RC=${RC}: $*"
	    #echo "ERROR: "`cat $RS_TMPERR`
	    echo;echo `date`"	TS	ERROR ${RS_NAME}	${RS_SAMPLE_NAME}"
	    qsub_exit $RC
	fi

	# handle success
	touch ${RS_TARGET_DONE}
	echo;echo `date`"	TS	DONE ${RS_NAME}	${RS_SAMPLE_NAME}"; 

    fi
    
    return $RC
}

