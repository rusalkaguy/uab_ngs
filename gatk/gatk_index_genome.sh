#!/bin/bash
FASTA=$1
if [ -z "$FASTA" ]; then echo "ERROR: $0 genome.fa"; exit 1; fi
if [ ! -e "$FASTA" ]; then echo "ERROR: file doesn't exist: $FASTA"; exit 1; fi
BASE=`dirname $FASTA`/`basename $FASTA .fa`
BASE=`dirname $FASTA`/`basename $FASTA .fas`
BASE=`dirname $FASTA`/`basename $FASTA .fasta`
if [ "$BASE" == `dirname $FASTA`/`basename $FASTA` ]; then
	echo "ERROR: un-recognized ending on $FASTA";
	exit 1
fi
module load ngs-ccts/samtools/0.1.19
samtools faidx $FASTA 
java -jar /share/apps/ngs-ccts/picard-tools/picard-tools-1.110/CreateSequenceDictionary.jar R= $FASTA O= $BASE.dict
