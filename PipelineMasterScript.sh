#!/bin/bash


#Initialize Folders
MASTER_DIR=`pwd`
SCRATCH="/scratch/user/$USER/"
RSTORE="/rstore/user/$USER/default"
CGI_SUBSECT_DIR="$MASTER_DIR/vcf"
COHORT_SRC_FOLDER="$RSTORE/CohortInfoFiles"
MERGE_DIR="$MASTER_DIR/gatk"
ANNOVAR_DIR="/home/$USER/uab_ngs/annovar"


#Initialize Variable for Gene List
GENE_LIST="$MASTER_DIR/TARGET_GENE_LIST.txt"



#1
#This first script uses tabix to extract gene regions the CGI VCFs
$CGI_SUBSECT_DIR/WUSTL.CGI.VCF.PREPARER.sh

#2
#The second script merges the CGI VCFs, then the WUSTL and CGI VCFs
$MERGE_DIR/gatk_merge_cg_wustl.sh


#3
#Make folders for all of the Gene Reports and copy the necessary files there
while read TARGGENE CHR LB UB; do
	MERGED_VCF_DIR="/rstore/user/$USER/default/VCFs/GenesOfInterest/$TARGGENE/MERGED"
	if [[ ! -d  $SCRATCH/GENE_REPORTS/$TARGGENE ]]; then
		mkdir -p $SCRATCH/GENE_REPORTS/$TARGGENE	
	fi
	cp $COHORT_SRC_FOLDER/* $SCRATCH/GENE_REPORTS/$TARGGENE
	cd $SCRATCH/GENE_REPORTS/$TARGGENE
	bash $ANNOVAR_DIR/doit.annovar.ics223.unfilt.cg.sh $MERGED_VCF_DIR/WUSTL600.ASW.CEU.Lupus.RA.YRI.ITGAM.hg19.vcf HC200memFS.CG \
cohort_pheno_compare_coding.wustl.cg.txt > C200memFS.CG.log

done < "$GENE_LIST"

#4 Still need to recursively make Excel Files
