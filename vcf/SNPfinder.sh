#!/bin/bash

BP="128578301"

echo $BP
for i in `seq 1 900`; do

	PERSON=$(grep "#CHROM" WUSTL600.ASW.CEU.Lupus.RA.YRI.IRF5.hg19.cgesrdVcgra.vcf | cut -f $i) 

	PLACE=$(grep $BP WUSTL600.ASW.CEU.Lupus.RA.YRI.IRF5.hg19.cgesrdVcgra.vcf | cut -f $i | grep 0/1)
	if [ $PLACE ]; then
		echo "$PERSON - heterozygote"
	fi 
        
	PLACE1=$(grep $BP WUSTL600.ASW.CEU.Lupus.RA.YRI.IRF5.hg19.cgesrdVcgra.vcf | cut -f $i | grep 1/0)
	if [ $PLACE1 ]; then
                echo "$PERSON - heterozygote"
        fi

        PLACE2=$(grep $BP WUSTL600.ASW.CEU.Lupus.RA.YRI.IRF5.hg19.cgesrdVcgra.vcf | cut -f $i | grep 1/1)
	if [ $PLACE2 ]; then
                echo "$PERSON - homozygote"
        fi
done
