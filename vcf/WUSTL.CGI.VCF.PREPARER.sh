#!/bin/bash

module load ngs-ccts/tabix/0.2.6

SCRIPT_DIR="/home/$USER/uab_ngs"
GENE_LIST="$SCRIPT_DIR/TARGET_GENE_LIST.txt"
BASE_DIR="/rstore/user/$USER/default/VCFs"
FOLDERS="$BASE_DIR/Folders.txt"

while read FOLDER; do

	echo "$FOLDER"
	VCF_LIST="$BASE_DIR/$FOLDER/VCF_Uniq_IDs.txt"
	
	while read GENE CHR LOCUS_BEGIN LOCUS_END; do

		#Make target directories
		mkdir -p $BASE_DIR/GenesOfInterest/$GENE

		while read VCF; do

                        #Save Header, because tabix does not.
                        zcat $BASE_DIR/$FOLDER/$VCF.vcf.gz | head -200 | grep ^# > $BASE_DIR/GenesOfInterest/$GENE/$FOLDER.$VCF.$GENE.vcf
                        #Extract All Variants for the Gene(s) of Interest and remove those that failed CGIs pipeline and add back to header
                        tabix $BASE_DIR/$FOLDER/$VCF.vcf.gz $CHR:$LOCUS_BEGIN-$LOCUS_END | grep $GENE | grep -v VQLOW | grep -v AMBIGUOUS \
                        >> $BASE_DIR/GenesOfInterest/$GENE/$FOLDER.$VCF.$GENE.vcf

		done < "$VCF_LIST"
	done < "$GENE_LIST"
done < "$FOLDERS"
