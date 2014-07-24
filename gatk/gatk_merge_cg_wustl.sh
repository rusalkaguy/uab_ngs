# merge CGI data
# output is MERGED/ASW.CEU.Lupus.RA.YRI.vcfBeta-ALL-ASM.$TARGGENE.vcf

TARG_GENE_FILE="/home/vlaufer/uab_ngs/TARGET_GENE_LIST.txt"


#if [ -z "$NSLOTS" ]; then
#	NSLOTS=1
#fi


while read TARGGENE CHR LOCUS_BEGIN LOCUS_END; do
	cd /rstore/user/$USER/default/VCFs/GenesOfInterest/$TARGGENE  #vlaufer still hardcoded
	~/uab_ngs/gatk/gatk_cohort_merge_vcf.sh
	ALL_CG="ASW.CEU.Lupus.RA.YRI"
	awk -v CHR=$CHR 'BEGIN{OFS="\t";IFS="\t";}("#CHROM"==$1){print "##AWKCommandLine=<ID=Cheap_b37_to_hg19>"}(CHR==$1){$1="chr"CHR}{print $0}' \
MERGED/$ALL_CG.vcfBeta-ALL-ASM.$TARGGENE.vcf > MERGED/$ALL_CG.vcfBeta-ALL-ASM.$TARGGENE.hg19.vcf
	echo "Created MERGED/$ALL_CG.vcfBeta-ALL-ASM.$TARGGENE.hg19.vcf"
# output is MERGED/ASW.CEU.Lupus.RA.YRI.vcfBeta-ALL-ASM.$TARGGENE.hg19.vcf


# cut WUSTL600 down to TARGGENE
# WARNING: should put filtered VCF in a subdirectory like we do for CG
	WUSTL_SRC_VCF=/home/curtish/scratch/kimberly/annovar/src_vcfs/bwamem_pad200_hg19/bwamem_pad200_hg19.fsok.vcf
	WUSTL_TARGET_VCF=MERGED/bwamem_pad200_hg19.fsok.$TARGGENE.vcf
	awk -v CHR=$CHR -v LB=$LOCUS_BEGIN -v UB=$LOCUS_END 'BEGIN{OFS="\t";IFS="\t";}(/^#/){print $0; next;}("chr"CHR==$1 && $2 >= LB && $2 <= UB){print $0;}' \
$WUSTL_SRC_VCF > $WUSTL_TARGET_VCF
	echo "Created $WUSTL_TARGET_VCF"

# merge in WUSTL600
	module load java/jre1.7.0_51
	GATK=/share/apps/ngs-ccts/GenomeAnalysisTK/GenomeAnalysisTK-3.0-0/GenomeAnalysisTK.jar 
# WARNING: should pull REF from the WUSTL VCF!!!
	REF=/scratch/share/public_datasets/ngs/databases/gatk_bundle/2.5/hg19/ucsc.hg19.fasta
	MERGE_OPT=PRIORITIZE # UNIQUIFY

	ALL_OUT=WUSTL600.$ALL_CG
	OUT_VCF=MERGED/$ALL_OUT.$TARGGENE.hg19.vcf
	java -Xmx20g -Xms20g -jar $GATK \
		-nt 1 \
		-R $REF \
		-T CombineVariants \
		-o $OUT_VCF \
		-genotypeMergeOptions $MERGE_OPT \
		-priority cg,wustl \
		--variant:cg MERGED/$ALL_CG.vcfBeta-ALL-ASM.$TARGGENE.hg19.vcf \
		--variant:wustl $WUSTL_TARGET_VCF \
		> JOBS/$ALL_OUT.merge.so \
		2> JOBS/$ALL_OUT.merge.se 
	
	jobs
	wait

# bgzip 
	echo "bgzip/tabix -p vcf $OUT_VCF"
	~/uab_ngs/tabix/tabix_bgzip.sh  "$OUT_VCF"
done < "$TARG_GENE_FILE"
