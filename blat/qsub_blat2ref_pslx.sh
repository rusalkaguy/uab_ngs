#! /bin/bash
# I've run this script passing parameters while qsubing 
#     qsub blat.sh $Filename $Samplename
# There are more blat running options and output options, 
# but the .pslx file is the one I've been using for visualization in IGV.
# The pslx file can be viewed in Excel and includes alot of other information 
# that isn't readily available in IGV.
#**** QSUB ****
#$ -N blat2ref_pslx
#$ -cwd
#$ -l h_rt=24:00:00,s_rt=24:00:00
#$ -l vf=31.5G,h_vmem=32G
#$ -M curtish@uab.edu -m beas #email at:  Begining, End, Abort, Suspend
#$ -j y # merge stderr into stdout
#$ -o jobs/$JOB_NAME.$JOB_ID.out
#**** QSUB ****
module load ngs-ccts/ucsc_kent/2014-03-05

#Variables inputted from parameters
OUT_DIR=$1 #Output directory
SAMPLE_NAME=$2
IN_FASTA=$3 #Path to reference genome
REF=$4


if [[ -z "$OUT_DIR" || -z "$SAMPLE_NAME" || -z "$IN_FASTA" || -z "$REF" ]]; then
	echo "SYNTAX: $0 OUT_DIR SAMPLE_NAME IN_FASTA REF"
	exit 1
fi
if [[ "$REF" == *.fa ]]; then
    if [  -e "${REF}.2bit" ]; then 
	REF=${REF}.2bit
    else
	echo "WARNING: no .2bit index for $REF "
    fi
fi

#use blat for long read alignment
CMD="blat -out=pslx $REF $IN_FASTA $OUT_DIR/${SAMPLE_NAME}.pslx"
echo "#CMD=$CMD"
$CMD; RC=$?
if [ $RC == 0 ]; then 
	touch $OUT_DIR/.done.${SAMPLE_NAME}.plsx
	echo "BLAT complete."
else
	echo "ERROR: BLAT FAILED: $RC"; 
fi
exit $RC
