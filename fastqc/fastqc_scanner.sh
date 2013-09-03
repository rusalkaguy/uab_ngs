#/bin/bash
#
# scan for and fastqc *.fastq
#
if [ -z "$1" ]; then
    echo SYNTAX: `basename $0`" [-find dir [ext] | filename [filenames...]"
    exit 1
fi

if [[ "$1" == "-find" || "$1" == "-f" ]]; then
    shift
    FASTQ_DIR=..
    if [ -n "$1" ]; then FASTQ_DIR=$1; fi
    FASTQ_EXT=.fastq
    if [ -n "$2" ]; then FASTQ_EXT=$2; fi
    
    FILES=`find $FASTQ_DIR -name "*$FASTQ_EXT"`
else
    FILES="$*"
fi

LOGS_DIR=logs
# output headers
echo "STATUS	FASTQ	DATE	ELAPSED_SEC"
mkdir -p $LOGS_DIR
# iterate over availabel files
for fastq in $FILES ; do
    DONE_FILE=$LOGS_DIR/out.`basename $fastq`.done
    START_SEC=`date +%s`

    # skip already processed files
    FASTQC_OUT=`basename $fastq $FASTQ_EXT`_fastqc
    if [[ ( ! -e "$DONE_FILE" || ! -e $FASTQC_OUT ) && ! -e ${FASTQC_OUT}.zip ]]; then  
	echo "START	   $fastq   "`date`
	OUT=$LOGS_DIR/out.`basename $fastq`.txt
	fastqc -t 10 --extract --outdir=. $fastq  > $OUT 2>&1
	RC=$?

	if [ $RC != 0 ]; then 
            STATUS="DONE-ERROR"
	else
	    STATUS=DONE
	    touch $DONE_FILE
	fi
    else  # skip
	STATUS=DONE-SKIP
    fi

    # CSV summary
    CSV_OUT=`basename $fastq $FASTQ_EXT`.csv 
    if [[ ! -e $CSV_OUT && -e $FASTQC_OUT.zip ]]; then  
	unzip -p $FASTQC_OUT.zip "*/fastqc_data.txt" | ~/ics/fastqc_to_csv/fastqc_to_csv.sh > $CSV_OUT
    fi
    if [[ ! -e $CSV_OUT && -e $FASTQC_OUT/fastqc_data.txt ]]; then  
	~/ics/fastqc_to_csv/fastqc_to_csv.sh $FASTQC_OUT/fastqc_data.txt > $CSV_OUT
    fi

    # log the completion status
    DONE_SEC=`date +%s`
    echo "$STATUS	   $fastq   "`date`" "`expr $DONE_SEC - $START_SEC`
done 
