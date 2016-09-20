# ======================================================================
# QC: heatmap of most variable genes
#
# CUSTOMIZE:
#    intgroups  = what sampleTable columns define experimental groups
#    ntop       = how many genes to use for qc_heatmap
# input: rld.gene
# ======================================================================
library(genefilter)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# CUSTOMIZE
ntop <- 50 
intgroups <- c("cell_line", "fraction")
ref_genes <- c()

# which version
load("dds.gene.union.ref_norm.RData", verbose=TRUE)
rld <- rld.gene

#
# find most variable genes
#
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),ntop)

# add reference genes, if you have a list in ref_genes
refGenes <- which(rownames(assay(rld)) %in% ref_genes)
topGenes <- c(refGenes,topVarGenes)

#
# build heatmap
#
mat <- assay(rld)[ topGenes, ]
mat <- mat - rowMeans(mat)

# label the reference genes, if any
rownames(mat)[rownames(mat) %in% ref_genes] = paste(
	rownames(mat)[rownames(mat) %in% ref_genes], 
	rep(" (ref)", length(ref_genes))
)

df <- as.data.frame(colData(rld)[,intgroups])
rownames(df) <- rld$sample
colnames(mat) <- rownames(df)


# 
# PLOT QC heatmap (basic)
#
library("pheatmap")

fig_title <- unlist(strsplit(x=c("QC",ntop,"most variable genes and ref_genes"),split=" "))

# OPEN PDF ----------------------------------------
pdf(paste(paste(fig_title,collapse="_"),".pdf",sep=""), height=9.5, width=7.5)

pheatmap(mat, annotation_col=df,cex=0.9, main=paste(fig_title, collapse=" "))

dev.off()
# CLOSE PDF  --------------------------------------


#
# PCA (basic)
#
fig_title <- unlist(strsplit(x=c("QC basic PCA",ntop,"most variable genes"),split=" "))

# OPEN PDF ----------------------------------------
pdf(paste(paste(fig_title,collapse="_"),".pdf",sep=""), height=9.5, width=7.5)
#
plotPCA(rld, intgroup = intgroups, ntop=ntop) #,main=paste(fig_title, collapse=" ") )
#
dev.off()
# CLOSE PDF  --------------------------------------

#
# PCA (advanced)
library(ggplot2); library(ggrepel)

fig_title <- unlist(strsplit(x=c("QC pretty PCA",ntop,"most variable genes"),split=" "))

(pcaData <- plotPCA(rld, intgroup = intgroups, ntop=ntop, returnData=TRUE))
pcaPercentVar <- round(100 * attr(pcaData , "percentVar"), digits=1)

# how to we get color/shape assignments dynamically!?
#aes(PC1, PC2, color=intgroups[1], shape=intgroups[2])) +
# see r0.*.R for defining color maps
#scale_color_manual(values=displayColors$cell_subset ) + 

# OPEN PDF ----------------------------------------
pdf(paste(paste(fig_title,collapse="_"),".pdf",sep=""), height=9.5, width=7.5)
#
ggplot(pcaData , aes(PC1, PC2, color=fraction, shape=cell_line)) +
	geom_point(size=3) +
	xlab(paste0("PC1: ",pcaPercentVar[1],"% variance")) +
	ylab(paste0("PC2: ",pcaPercentVar[2],"% variance")) +
	ggtitle(paste(fig_title, collapse=" ")) 
#
dev.off()
# CLOSE PDF ----------------------------------------

