#!/bin/bash
module load ngs-ccts/ucsc_kent/2014-03-05

# check syntax
FA=$1
if [[ -z "$FA" || "$FA" != *.fasta ]]; then 
    echo "SYNTAX: $0 src.fasta"
    echo "Creates src/ and splits fasta into individual sequences"
    exit 1
fi

# do work
FDIR=`dirname $FA`/`basename $FA .fasta`
mkdir -p $FDIR
faSplit byname $FA $FDIR/
# fix names not to have |.s - replace with _
cd $FDIR
for x in *.fa; do 
    y=`echo $x | sed 's/|\././;s/|/_/g'`; 
    CMD="mv $x $y"; 
    if [ ! -e "$y" ]; then 
	echo $CMD;
	$CMD; 
    else 
	echo "skip $y"; 
    fi; 
done

