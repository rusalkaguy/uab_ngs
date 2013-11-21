#!/bin/bash
########################################################################
# 
# Run an ANNOVAR pipeline on a SAMPLE
# Filter variants by a list of conditions
#
# http://www.openbioinformatics.org/annovar/annovar_accessary.html#var
#
# -operation argument instruct what operation are used: gene-based (g), reverse region-based (rr), region-based (r), filter-based (f), model-based (m), respectively.
#
########################################################################
#$ -S /bin/bash
#$ -cwd # remember what dir I'm launched in
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -m beas  #email at:  Begining, End, Abort, Suspend
#$ -M $USER@uab.edu
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

if [ -n "$1" ]; then IN=$1; fi
if [ -z "$IN" ]; then IN=wustl_pad200_b37.hg19.ESRD_SLE_1.avinput; fi
if [ -z "$OUT" ]; then OUT=`dirname $IN`/`basename $IN .avinput`.reduce; fi

echo "IN=$IN"
echo "OUT=$OUT"
FLAGS=

# gene annotations
OPERATOR+=g; PROTOCOL+=nonsyn_splicing; OP_ARGS+=  # keep only nonsyn and splice-site variants
# region annotations

#OPERATOR+=,r; PROTOCOL+=,phastConsElements46way; OP_ARGS+=,
#OPERATOR+=,r; PROTOCOL+=,tfbs; OP_ARGS+=,
#OPERATOR+=,r; PROTOCOL+=,dgv; OP_ARGS+=,
#OPERATOR+=,r; PROTOCOL+=,gwascatalog; OP_ARGS+=,
#OPERATOR+=,r; PROTOCOL+=,wgEncodeBroadHmmGm12878HMM; OP_ARGS+=, # blood cell line promoters/repressor

# filter annotations

OPERATOR+=,f; PROTOCOL+=,1000g2012apr_eur; OP_ARGS+=, #OP_ARGS+=",-maf 0.05", # 1000 Genomes Project (2012 April) EUR allele freq
  FLAGS+=" --aaf_threshold 0.05"

  # -maf 0.05 # drop all variants with MAF > 0.05 in 100G_euro
  # -maf 0.05 -reverse # drop all variants with MAF < 0.05 in 1000G_euro

#OPERATOR+=,f; PROTOCOL+=,1000g2012apr_all; OP_ARGS+=, # 1000 Genomes Project (2012 April) ALL-samples allele freq
#OPERATOR+=,f; PROTOCOL+=,popfreq_max; OP_ARGS+=, # Max pop freq across 
#OPERATOR+=,f; PROTOCOL+=,popfreq_all; OP_ARGS+=, # DOWNLOAD MISSING - CSV allele/population freq values: PopFreqMax 1000G2012APR_ALL 1000G2012APR_AFR 1000G2012APR_AMR 1000G2012APR_ASN 1000G2012APR_EUR ESP6500si_AA ESP6500si_EA CG46 NCI60 SNP137 COSMIC65 DISEASE
#OPERATOR+=,f; PROTOCOL+=,esp6500si_all; OP_ARGS+=, # NHLBI exome over 6000 healthy + disease http://www.openbioinformatics.org/annovar/annovar_filter.html#esp
#OPERATOR+=,f; PROTOCOL+=,esp6500si_ea; OP_ARGS+=, # European Ancestry
#OPERATOR+=,f; PROTOCOL+=,cg46; OP_ARGS+=, # Complete Genomics 46 unrelated healthy people (WGS)
#OPERATOR+=,f; PROTOCOL+=,snp129; OP_ARGS+=, # last pre 1000g dbSNP
#OPERATOR+=,f; PROTOCOL+=,snp137; OP_ARGS+=, # latest
#OPERATOR+=,f; PROTOCOL+=,snp137NonFlagged; OP_ARGS+=, # Flagged SNPs include SNPs < 1% minor allele frequency (MAF) (or unknown), mapping only once to reference assembly, flagged in dbSnp as "clinically associated". 
OPERATOR+=,f; PROTOCOL+=,ljb2_sift; OP_ARGS+=, # remove "benign" snps

OPERATOR+=,f; PROTOCOL+=,ljb2_pp2hdiv; OP_ARGS+=, # remove "benign" snps
  # ljb2_pp2hvar should be used for diagnostics of Mendelian diseases, which requires distinguishing mutations with drastic effects from all the remaining human variation, including abundant mildly deleterious alleles.The authors recommend calling "probably damaging" if the score is between 0.909 and 1, and "possibly damaging" if the score is between 0.447 and 0.908, and "benign" is the score is between 0 and 0.446
  # ljb2_pp2hdiv should be used when evaluating rare alleles at loci potentially involved in complex phenotypes, dense mapping of regions identified by genome-wide association studies, and analysis of natural selection from sequence data. The authors recommend calling "probably damaging" if the score is between 0.957 and 1, and "possibly damaging" if the score is between 0.453 and 0.956, and "benign" is the score is between 0 and 0.452.

# models

# OPERATOR+=,m; PROTOCOL+=,dominant; OP_ARGS+=, # 
echo -n "START "; date

CMD="variants_reduction.pl $IN $ANNOVAR_HG19_DB --buildver hg19 \
    --outfile $OUT \
    --protocol  $PROTOCOL \
    --operation $OPERATOR \
    --argument  $OP_ARGS \
    $FLAGS
"
echo $CMD
$CMD
RC=$?
echo "RC=$RC"

echo -n "COMPLETE "; date
exit $RC 2>&1 > /dev/null


