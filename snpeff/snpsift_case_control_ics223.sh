#!/bin/bash
echo -n "START "; date
IN=$PWD/uab_pad200_hg19.vcf; NAME=HC200
if [ -n "$1" ]; then IN="$1"; fi
if [ -n "$2" ]; then NAME="$2"; fi
IN_BASE=`basename $IN`
OUT_PHE12=`basename $IN .vcf`.caseControl.vcf
OUT_X=`basename $IN .vcf`.caseControl.txt

#***************************************************************
echo "generate .tfam files from VCF $IN"
if [ ! -e $IN_BASE.sleVnorm.tfam ]; then 
    ~/uab_ngs/snpeff/vcf2tped_ics223.sh $IN
else
    echo SKIP
fi

PHE=sleVnorm #**************************************************
echo "caseControl first phenotype: $PHE"
OUT_PHE1=.$OUT.$PHE.vcf
if [ ! -e $OUT_PHE1 ]; then 
java -jar /share/apps/ngs-ccts/snpEff_3_3/SnpSift.jar caseControl \
    -tfam $IN_BASE.$PHE.tfam \
    -name _${NAME}_$PHE \
    $IN > $OUT_PHE1
else
    echo SKIP
fi
PHE=esrdVnorm #**************************************************
echo "caseControl first phenotype: $PHE"
if [ ! -e $OUT_PHE12 ]; then
    java -jar /share/apps/ngs-ccts/snpEff_3_3/SnpSift.jar caseControl \
	-tfam $IN_BASE.$PHE.tfam \
	-name _${NAME}_$PHE \
	$OUT_PHE1 > $OUT_PHE12

    echo "Created $OUT"
else
    echo SKIP
fi

echo "Excel Extract" #**************************************************
echo "Creating $OUT_X"
if [ ! -e $OUT_X ]; then 

    # extract field list
    echo "Extract field list"
    COL_LIST=`grep -v "^#" $OUT_PHE12 | grep Controls | head -n 1 \
              | cut -f 8 \
              | tr ";" "\n" \
              | egrep "^(Cases|Control|CC_)" \
              | cut -d = -f 1 \
              | sort \
              | uniq`
    echo "    "`echo $COL_LIST | wc -w`" fields ("`echo $COL_LIST | tr "	" "\n"| cut -d _ -f 1 | sort | uniq | perl -pe 's/\n//;'`")"
    echo "    $COL_LIST"
    echo "Subscript field list"
    for f in $COL_LIST; do
	if [[ "$f" == Cases* || "$f" == Controls* ]]; then
	    COL_LIST_SUB="$COL_LIST_SUB ${f}[0] ${f}[1] ${f}[2]"
	    MOD_LIST="$MOD_LIST $f"
	else
	    COL_LIST_SUB="$COL_LIST_SUB $f"
	fi
    done
    echo "run SnpSift.jar/extractFields"
    java -jar /share/apps/ngs-ccts/snpEff_3_3/SnpSift.jar extractFields \
	$OUT_PHE12 \
	CHROM POS REF ALT QUAL FILTER AD DP GQ GT PL AC AF AN BaseQRankSum ClippingRankSum DB DP DS FS HaplotypeScore InbreedingCoeff MLEAC MLEAF MQ MQ0 MQRankSum QD ReadPosRankSum $COL_LIST_SUB \
	| perl -pe 's/norm\[0\]/norm_HOM/g;s/norm\[1\]/norm_HET/g;s/norm\[2\]/norm_TOT/g;' \
	> $OUT_X
else
    echo "SKIP"
fi

echo -n "COMPLETE "; date
