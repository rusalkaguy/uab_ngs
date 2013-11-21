#!/bin/bash
#*** QSUB FLAGS ***
## avoid SSG nodes, which cause "mem map file errors" (SNPeff annotate & annovar/convert2annovar.pl)
## -q all.q,sipsey.q
# *** std qsub ****
#$ -S /bin/bash
#$ -cwd # remember what dir I'm launched in
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -m beas  #email at:  Begining, End, Abort, Suspend
# *** output logs ***
#$ -j y # merge stderr into stdout
#$ -o $JOB_NAME.$JOB_ID.out
# run time
#$ -l h_rt=119:00:00 -l s_rt=120:55:00 
#$ -l vf=6G -l h_vmem=6" # 1 CPUs, 6G RAM
#*** END QSUB ****
# libraries to load
. /etc/profile.d/modules.sh          # enable module loading
module load ngs-ccts/annovar/2013.09.11
module list

if [ -n "$1" ]; then VCF=$1; fi
if [ -z "$VCF" ]; then VCF=uab_pad200_hg19.vcf; fi
if [[ "$VCF" == *hc_* ]]; then 
    # UAB generated HaploType caller - includes sites that failed QC, FILTER column is set appropriately
    FILTER="-filter PASS"
else
    # WUSTL - these VCFs already have the QC failures removed, but do NOT have PASS in the filter column, only "."
    FILTER=""
fi

echo "VCF=$VCF"
echo "FILTER=$FILTER"

DEST=vcf4old.inc.com/`basename $VCF .vcf`.avinput
mkdir -p `dirname $DEST`
echo -n "START "; date;
echo -n "HOST  $HOSTNAME aka "; hostname
CMD="convert2annovar.pl -format vcf4old $VCF -outfile $DEST $FILTER -include -comment "
echo $CMD
$CMD; RC=$?
echo -n "DONE RC=$RC "; date

DEST=vcf4old.allallele.inc.com/`basename $VCF .vcf`.avinput
mkdir -p `dirname $DEST`
echo -n "START "; date;
echo -n "HOST  $HOSTNAME aka "; hostname
CMD="convert2annovar.pl -format vcf4old $VCF -outfile $DEST $FILTER -include -comment -allallele"
echo $CMD
$CMD; RC=$?
echo -n "DONE RC=$RC "; date
