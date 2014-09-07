#!/bin/bash

module load ngs-ccts/tabix/0.2.6

# bed file of high AAF dbSNP138 variants
TP=~/scratch/kimberly/ucsc/snp138/snp138.no_exceptions.weight1.pad0.txt.aaf20.gz  # bed-like file
# VCFs (bgzip;tabix -p vcf) prepared
WU=wustl_pad200_b37/wustl_pad200_b37.hg19.vcf.gz;
WUQC=filtered/wustl_pad200_b37.hg19.dp13.gq30.miss.vcf.gz;
UAB=bwamem_pad200_hg19/bwamem_pad200_hg19.raw.vcf.gz;
UABPASS=bwamem_pad200_hg19/bwamem_pad200_hg19.vcf.bgz
UABQC=filtered/bwamem_pad200_hg19.fs500.adp.dp13.gq30.miss.pass.left.dbsnp.vcf.gz;
(echo "chr dbSNP WU UAB UAB.PASS UAB.QC WU.QC"; 
for c in {1..22} X; do 
    echo \
	chr$c \
	`zcat $TP | grep -w "^chr$c" | wc -l` \
        `tabix $WU  -B  <(zcat $TP | grep -w "^chr$c") | wc -l` \
        `tabix $UAB -B  <(zcat $TP | grep -w "^chr$c") | wc -l` \
        `tabix $UABPASS -B  <(zcat $TP | grep -w "^chr$c") | wc -l` \
        `tabix $UABQC -B  <(zcat $TP | grep -w "^chr$c") | wc -l` \
        `tabix $WUQC -B  <(zcat $TP | grep -w "^chr$c") | wc -l` 
done) \
| column -t 
