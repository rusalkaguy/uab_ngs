#!/bin/bash
########################################################################
# 
# Run an ANNOVAR pipeline on a SAMPLE
# Summarize into EXCEL format (TSV)
#
# -operation argument instruct what operation are used: gene-based (g), reverse region-based (rr), region-based (r), filter-based (f), model-based (m), respectively.
#
########################################################################
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

if [ -z "$IN" ]; then IN=wustl_pad200_b37.hg19.ESRD_SLE_1.avinput; fi
if [ -z "$OUT" ]; then OUT=$IN; fi

OUT=$IN

# gene annotations
OPERATOR+=g; PROTOCOL+=refGene; OP_ARGS+=

# region annotations

OPERATOR+=,r; PROTOCOL+=,phastConsElements46way; OP_ARGS+=,
OPERATOR+=,r; PROTOCOL+=,tfbs; OP_ARGS+=,
OPERATOR+=,r; PROTOCOL+=,dgv; OP_ARGS+=,
OPERATOR+=,r; PROTOCOL+=,gwascatalog; OP_ARGS+=,
OPERATOR+=,r; PROTOCOL+=,wgEncodeBroadHmmGm12878HMM; OP_ARGS+=, # blood cell line promoters/repressor

# filter annotations

OPERATOR+=,f; PROTOCOL+=,cg46; OP_ARGS+=, # blood cell line promoters/repressor
OPERATOR+=,f; PROTOCOL+=,esp6500si_all; OP_ARGS+=, # NHLBI exome over 6000 healthy + disease http://www.openbioinformatics.org/annovar/annovar_filter.html#esp
OPERATOR+=,f; PROTOCOL+=,esp6500si_ea; OP_ARGS+=, # European Ancestry
OPERATOR+=,f; PROTOCOL+=,cg46; OP_ARGS+=, # Complete Genomics 46 unrelated healthy people (WGS)
OPERATOR+=,f; PROTOCOL+=,snp129; OP_ARGS+=, # last pre 1000g dbSNP
OPERATOR+=,f; PROTOCOL+=,snp137; OP_ARGS+=, # latest
OPERATOR+=,f; PROTOCOL+=,snp137NonFlagged; OP_ARGS+=, # Flagged SNPs include SNPs < 1% minor allele frequency (MAF) (or unknown), mapping only once to reference assembly, flagged in dbSnp as "clinically associated". 
OPERATOR+=,f; PROTOCOL+=,ljb2_all; OP_ARGS+=, # 

# DOWNLOAD?: popfreq_all - CSV allele/population freq values: PopFreqMax 1000G2012APR_ALL 1000G2012APR_AFR 1000G2012APR_AMR 1000G2012APR_ASN 1000G2012APR_EUR ESP6500si_AA ESP6500si_EA CG46 NCI60 SNP137 COSMIC65 DISEASE
# DOWNLOAD 1000g2012apr (1000g2012apr,1000g2012apr_EA/EUR) - pop freq  - can apply -maf 0.05 argument
echo -n "START "; date

table_annovar.pl $IN $ANNOVAR_HG19_DB --buildver hg19 \
    -protocol  $PROTOCOL \
    -operation $OPERATOR \
    -arg       $OP_ARGS \
    -nastring NA \
    -outfile $OUT
RC=$?
echo "RC=$RC"

echo -n "COMPLETE "; date
exit $RC 2>&1 > /dev/null


