#!/bin/bash
# I've run this script passing parameters while qsubing 
#     qsub qsub_kraken_filter.sh $OUT_DIR $SAMPLE_NAME $IN_FASTA $TAXON $KRAKEN_DB_NAME
#**** QSUB ****
#$ -N kraken_filter
#$ -cwd # remember what dir I'm launched in 
#$ -l h_rt=5:00:00,s_rt=5:00:00
#$ -pe smp 4 -l vf=5.75G,h_vmem=6G
#$ -M curtish@uab.edu -m beas #email at:  Begining, End, Abort, Suspend
#$ -j y # merge stderr into stdout
#$ -o jobs/$JOB_NAME.$JOB_ID.out
#**** QSUB ****
module load ngs-ccts/kraken/0.10.3-beta

#Variables inputted from parameters
OUT_DIR=$1 #Output directory
SAMPLE_NAME=$2
IN_FASTA=$3 #Path to reference genome
TAXON_FILE=$4
KRAKEN_DB_NAME=$5

KRAKEN_DB=/scratch/share/public_datasets/ngs/databases/kraken/$KRAKEN_DB_NAME
if [[ -z "$OUT_DIR" || -z "$SAMPLE_NAME" || -z "$IN_FASTA" || -z "$TAXON_FILE" ]]; then
	echo "SYNTAX: $0 OUT_DIR SAMPLE_NAME IN_FASTA NCBI_TAXON_LIST [KRAKEN_DB]"
	exit 1
fi
if [ ! -e "$IN_FASTA" ]; then echo "FILE MISSING: IN_FASTA=$IN_FASTA"; exit 1; fi
if [ ! -e "$KRAKEN_DB" ]; then echo "FILE MISSING: KRAKEN_DB=$KRAKEN_DB"; exit 1; fi
if [[ "$TAXON_FILE" == "none" || "$TAXON_FILE" == "classified" || "$TAXON_FILE" == "unclassified" ]]; then 
    TAXON=$TAXON_FILE
else
    if [ ! -e "$TAXON_FILE" ]; then echo "FILE MISSING: TAXON_FILE=$TAXON_FILE"; exit 1; fi
    TAXON=`grep "^##ROOT_ID" $TAXON_FILE | cut -f 2` 
    if [ -z "$TAXON" ]; then echo "##ROOT_ID missing from TAXON_FILE=$TAXON_FILE"; exit 1; fi
fi
mkdir -p $OUT_DIR
if [ -z "$NSLOTS" ]; then 
    echo "WARNING: setting NSLOTS=1"
    NSLOTS=1
fi

#
#use Kraken for (FAST) metagenomic analysis
#

# 
# guess format of input file
# 
OUT_SUFFIX=fasta
INPUT_FLAGS=
# .dat  is a hack for dealing with Galaxy files
if [[ "$IN_FASTA" == *.fastq* || "$IN_FASTA" == *.fq* || "$IN_FASTA" == *.dat ]]; then INPUT_FLAGS="$INPUT_FLAGS --fastq-input"; OUT_SUFFIX="fastq"; fi
if [[ "$IN_FASTA" == *.gz ]]; then INPUT_FLAGS="$INPUT_FLAGS --gzip-compressed"; fi


# 
# output
#
OUTPUT=$OUT_DIR/$SAMPLE_NAME.kraken.txt
STATS=$OUT_DIR/$SAMPLE_NAME.kraken.stats
CLASSIFIED=$OUT_DIR/$SAMPLE_NAME.classified.$OUT_SUFFIX
UNCLASSIFIED=$OUT_DIR/$SAMPLE_NAME.unclassified.$OUT_SUFFIX
CONTIGS=$OUT_DIR/$SAMPLE_NAME.kraken.$TAXON.contigs.txt

# 
# format command
#
OUTPUT_FLAGS=
if [[ "$TAXON" == "classified" ]]; then OUTPUT_FLAGS="$OUTPUT_FLAGS --classified-out $CLASSIFIED "; fi
if [[ "$TAXON" == "unclassified" ]]; then OUTPUT_FLAGS="$OUTPUT_FLAGS --unclassified-out $UNCLASSIFIED "; fi
if [[ "$TAXON" == "none" ]]; then OUTPUT_FLAGS="$OUTPUT_FLAGS"; fi
CMD="kraken \
   -db $KRAKEN_DB \
   --threads $NSLOTS  \
   --preload \
   $INPUT_FLAGS $IN_FASTA \
   --output $OUTPUT \
   $OUTPUT_FLAGS \
   2> $STATS \
"
#CMD="echo skip kraken"
echo "#CMD=$CMD"
eval $CMD; RC=$?
if [ "$RC" == 0 ]; then 
	date >  $OUT_DIR/.done.${SAMPLE_NAME}.done
	echo "KRAKEN complete."
else
	echo "ERROR: KRAKEN FAILED: $RC"; 
	cat $STATS
	exit $RC
fi

if [[ -e "$TAXON_FILE" ]]; then 

    #
    # extract HCMV (species by NCBI TAXID) contigs
    #
    AWK_FILT=$OUT_DIR/taxon_filter.awk
    cat > $AWK_FILT <<EOF
# filter kraken output list based on taxid list
BEGIN { 
    FS="\t"; OFS="\t";
    #debug=0
}
# load taxa-of-interest table: c1 = taxid
(tax_filename == "") {
  if(1==debug){print "Loading " FILENAME}
  tax_filename=FILENAME
}
(FILENAME==tax_filename && ! /^#/ ) {
    taxarr[\$1]=1;
    if(1==debug){print "\tloaded["\$1"] => " \$0}
    next;
}
(tax_filename != FILENAME && tax_filename != "" && filter_start == "") {
  if(1==debug){print "Filtering " FILENAME}
  filter_start = FILENAME
}
# get seqIds from kreken.txt for matching taxa
(FILENAME!=tax_filename && "C"==\$1 && 1==taxarr[\$3]) {
  if(1==debug){
    print \$2, \$3, taxarr[\$3] # seqID, taxid, flag
  } else {
    print \$2  # print fasta seq ID
  }
}
END {
  if(1==debug){print "done"}
}
 
EOF
    echo -n "" > $CONTIGS 
    for x in `awk -f $AWK_FILT $TAXON_FILE $OUTPUT`; do \
	#echo "    extract seq $x from $CLASSIFIED"
	grep -w -A 1 "^>$x" $CLASSIFIED >> $CONTIGS \
    ; done
    echo "TAXON $TAXON CONTIGS: " `grep -c "^>" $CONTIGS`
fi
exit $RC
