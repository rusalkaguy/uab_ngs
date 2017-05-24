# ======================================================================
# run DESeq2 basic analysis and setup for QC
#
# MUST CUSTOMIZE: DESeqDataSet: design = ~ ??? 
# inputs: se.gene
# outputs: dds.gene, rld.gene
# ======================================================================

# pre-load libraries
library(DESeq2)

#
# read counts are loaded
# 
load("se.gene.union.RData", verbose=TRUE)
se <- se.gene

#
# create the DESeq2 dataset
#
# Note: default settings of the package, you should put the
#   * variable of interest at the END of the formula
#   * and make sure the control level is the 1ST LEVEL.

library(DESeq2)
#dds <- DESeqDataSet(se, design= ~ cell_line + fraction)
# basic version always works, but less powerful
dds <- DESeqDataSet(se, design= ~ sample)  
# DESeq2 1.3.6 Pre-Filtering; remove rows with low expression
dds<- dds[rowSums(counts(dds)) > 1, ]
# 2-3 minutes to build model
dds <- DESeq(dds)  

#
# create RLOG version for QC
# blind = don't use design
#
rld <- rlog(dds, blind=TRUE)

#
# save
#
dds.gene <- dds
rld.gene <- rld
save(file="dds.gene.union.RData", dds.gene, rld.gene)

