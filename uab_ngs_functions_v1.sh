#!/bin/bash
#
# shared subroutine used by the qsub_* scripts in these subdirectories
#

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

    # log start/skip
    if [[ ( -e "${RS_TARGET}" || -n "$DONE_ONLY" ) && ( -e "${RS_TARGET}.done" ) ]]; then
	echo;echo `date`"	TS	SKIP ${RS_NAME}	${RS_SAMPLE_NAME}"; 
    else
	echo;echo `date`"	TS	START ${RS_NAME}	${RS_SAMPLE_NAME}"; 
	
	# run the cmd
	RS_TMPERR=`mktemp`
	if [[ -z "$RS_STDOUT" || "-" == "$RS_STDOUT" ]]; then  
	    # stdout to stdout
	    echo "# CMD: $*" 
	    $* 2>$RS_TMPERR
	else
    	    # stdout to file
	    echo "# CMD: $* 1>$RS_STDOUT" 
	    $* 2>$RS_TMPERR 1>$RS_STDOUT
	fi
    fi
    #  handle errors
    RC=$?
    if [ $RC != 0 ]; then 
	echo "ERROR: $*"
	echo "ERROR: "`cat $RS_TMPERR`
	echo;echo `date`"	TS	ERROR ${RS_NAME}	${RS_SAMPLE_NAME}"
	qsub_exit $RC
    fi
    # handle success
    touch ${RS_TARGET}.done
    echo;echo `date`"	TS	DONE ${RS_NAME}	${RS_SAMPLE_NAME}"; 
    
    return $RC
}

