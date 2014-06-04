#!/bin/bash
# I've run this script passing parameters while qsubing 
#     qsub blat.sh $Filename $Samplename
# There are more blat running options and output options, 
# but the .pslx file is the one I've been using for visualization in IGV.
# The pslx file can be viewed in Excel and includes alot of other information 
# that isn't readily available in IGV.
#**** QSUB ****
#$ -N kraken_mini_filter
#$ -cwd # remember what dir I'm launched in 
#$ -l h_rt=5:00:00,s_rt=5:00:00
#$ -l vf=9.5G,h_vmem=10G
#$ -M curtish@uab.edu -m beas #email at:  Begining, End, Abort, Suspend
#$ -j y # merge stderr into stdout
#$ -o jobs/$JOB_NAME.$JOB_ID.out
#**** QSUB ****
module load ngs-ccts/kraken/0.10.3-beta

#Variables inputted from parameters
OUT_DIR=$1 #Output directory
SAMPLE_NAME=$2
IN_FASTA=$3 #Path to reference genome
TAXON=$4

KRAKEN_MINIDB=/scratch/share/public_datasets/ngs/databases/kraken/minikraken_20140104

if [[ -z "$OUT_DIR" || -z "$SAMPLE_NAME" || -z "$IN_FASTA" || -z "$TAXON" ]]; then
	echo "SYNTAX: $0 OUT_DIR SAMPLE_NAME IN_FASTA NCBI_TAXON_ID"
	exit 1
fi
mkdir -p $OUT_DIR
if [ -z "$NSLOTS" ]; then 
    echo "WARNING: setting NSLOTS=1"
    NSLOTS=1
fi
 
#
#use Kraken for (FAST) metagenomic analysis
#
OUTPUT=$OUT_DIR/$SAMPLE_NAME.kraken.txt
STATS=$OUT_DIR/$SAMPLE_NAME.kraken.stats
CLASSIFIED=$OUT_DIR/$SAMPLE_NAME.classified.txt
UNCLASSIFIED=$OUT_DIR/$SAMPLE_NAME.unclassified.txt
CONTIGS=$OUT_DIR/$SAMPLE_NAME.kraken.$TAXON.contigs.txt

CMD="kraken \
   -db $KRAKEN_MINIDB \
   --threads $NSLOTS  \
   --preload \
   $IN_FASTA \
   --output $OUTPUT \
   --unclassified-out $UNCLASSIFIED  \
   --classified-out $CLASSIFIED \
   2> $STATS \
"
echo "#CMD=$CMD"
eval $CMD; RC=$?
if [ $RC == 0 ]; then 
	date >  $OUT_DIR/.done.${SAMPLE_NAME}.done
	echo "KRAKEN complete."
else
	echo "ERROR: KRAKEN FAILED: $RC"; 
	cat $STATS
	exit $RC
fi

#
# extract HCMV (species by NCBI TAXID) contigs
#
echo > $CONTIGS
for x in `awk '("'$TAXON'"==$3){print "^>" $2 "$";}' $OUTPUT`; do \
    grep -A 1 "$x" $CLASSIFIED >> $CONTIGS \
; done
echo "TAXON $TAXON CONTIGS: " `grep -c "^>" $CONTIGS`
exit $RC
