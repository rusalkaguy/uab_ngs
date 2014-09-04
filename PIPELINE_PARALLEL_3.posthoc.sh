#!/bin/bash

#Set Dirs
BASE_DIR=`pwd`
SCRATCH="/scratch/user/$USER"
RSTORE="/rstore/user/$USER/default"
CGI_SUBSECT_DIR="$MASTER_DIR/vcf"
COHORT_SRC_FOLDER="$RSTORE/CohortInfoFiles"
ANNOVAR_DIR="/home/$USER/uab_ngs/annovar"

#Initialize Variable for Gene List
GENE_LIST="$BASE_DIR/uab_ngs/TARGET_GENE_LIST.txt"

#Set Amount of RAM
RAM="4"
RAM_FREE=$(($RAM - 1))

#Make Dirs for scripts
if [[ ! -d $BASE_DIR/JOBS ]]; then
	mkdir -p $BASE_DIR/JOBS
	mkdir -p $BASE_DIR/JOBS/OUT
fi

#Make Dirs where slave scripts will reside
JOBS_DIR=$BASE_DIR/JOBS
ERROR_DIR=$JOBS_DIR/OUT

#Start Outer Loops
#head -2 "$GENE_LIST" |
cat "$GENE_LIST" |
while read TARGGENE CHR LB UB; do

	QSCRIPT=$JOBS_DIR/$TARGGENE.RV.POSTHOC1.sh
	DONE_FILE=$QSCRIPT.done
	QSUB_NAME=$TARGGENE.RV.POSTHOC1.qsub
	cat <<EOT > $QSCRIPT
#!/bin/sh
#$ -cwd
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -j y # merge stderr into stdout
#$ -o OUT/\$JOB_NAME.\$JOB_ID
# email at:  Beginning, End, Abort, Suspend
#$ -m beas -M $USER@uab.edu,curtish@uab.edu
#$ -N $QSUB_NAME
# *** DRMAA ****
#$ -l vf=${RAM_FREE}.9G -l h_vmem=${RAM}G  # observed max of 1.87G for chr1
#$ -l h_rt=10:00:00          # observed max of 30 min for chr1
date

MERGED_VCF_DIR="/rstore/user/vlaufer/default/VCFs/GenesOfInterest/$TARGGENE/MERGED"
#Run the Script
cd $SCRATCH/GENE_REPORTS/$TARGGENE
pwd
# convert VCF, with groups (uab_ngs: f6ceac76a485117bf009301d0b5e0dd86ebc25fb)
VCF=\`\ls *.caseControl.vcf\`
OUT=\${VCF}.genotypes.txt
$BASE_DIR/uab_ngs/vcf/vcf_genotype_export.pl \
        --groups $COHORT_SRC_FOLDER/cohort_sample_list.wustl.cg.txt \
        < \$VCF \
        > \$OUT


RC=\$?
echo `date` "RC = \$RC"
if [ "\$RC" == 0 ]; then
	date > $DONE_FILE
fi
exit \$RC 2>&1 > /dev/null
EOT
	chmod +x $QSCRIPT
        #Queue up for submission to compute node
	if [ -e "$DONE_FILE" ]; then
		echo "SKIPPED: Out file may already exist; job is $QSCRIPT"
	else
		(cd $JOBS_DIR; qsub $QSCRIPT) 
	fi
	
        #Close the loops opened before the slave script is written
done 
