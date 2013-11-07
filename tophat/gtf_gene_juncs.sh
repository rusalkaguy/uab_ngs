# get PKHD1 canonical junction lists for mm9 and mm10
module load galaxy/galaxy-command-line
GENE=pkhd1
echo "GENE=$GENE"
for MM in mm9 mm10; do
    GTF=/scratch/share/public_datasets/ngs/gtfs/${MM}_patched.gtf
    echo "GTF=$GTF"
    OUT=${MM}.${GENE}.juncs
    echo "extracting junctions for $GENE in $MM"
    /share/apps/ngs-ccts/tophat-2.0.9.Linux_x86_64/gtf_juncs \
	<(grep -i '"'$GENE'"' $GTF) \
	> $OUT
    wc -l $OUT
done
