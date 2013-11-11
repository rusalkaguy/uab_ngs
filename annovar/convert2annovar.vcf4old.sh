#!/bin/bash
#*** QSUB FLAGS ***
# avoid SSG nodes, which cause "mem map file errors" (SNPeff annotate & annovar/convert2annovar.pl)
#$ -q all.q,sipsey.q
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
VCF=gatk2.5-2-var_filt-snps_indels_filtered-all.vcf


DEST=vcf4old.inc.com/uab_pad200_hg19 
mkdir -p `dirname $DEST`
echo -n "START "; date;
echo -n "HOST  $HOSTNAME aka "; hostname
cd /scratch/user/curtish/kimberly_lupus/annovar/src_vcfs/uab_pad200_hg19
convert2annovar.pl -format vcf4old $VCF -outfile $DEST -include -comment  
echo -n "DONE  "; date

DEST=vcf4old.allallele.inc.com/uab_pad200_hg19 
mkdir -p `dirname $DEST`
echo -n "START "; date;
echo -n "HOST  $HOSTNAME aka "; hostname
cd /scratch/user/curtish/kimberly_lupus/annovar/src_vcfs/uab_pad200_hg19
convert2annovar.pl -format vcf4old $VCF -outfile $DEST -include -comment  -allallele
echo -n "DONE  "; date
