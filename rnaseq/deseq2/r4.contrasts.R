#
#
dds <- dds.gene
# collapse replicates (none in this case!)
ddsc <- collapseReplicates(dds, groupby=paste(dds$cell_line,dds$fraction,sep="."), run=dds$sample)
rldc <- rlog(ddsc,blind=FALSE) # non-blinded for display actual results (not QC)


alpha <- 0.10
contrast_list <- list(
	c("cell_line", "Ctrl", "shS25"),
	c("fraction","Mono","Poly")
)

contrast <- contrast[[1]]

# get the results (fold change)
res <- results(ddsc, alpha=alpha, contrast)

#----------------------------------------------------------------------
# volcano plot
#----------------------------------------------------------------------

#save to PDF
pdf(paste(figure_title,".pdf",sep=""), height=8, width=8)

# plot expr vs fold-change, color by pvalue
(lfcLim=c(min(res$log2FoldChange),max(res$log2FoldChange)))
(logpLim=c( min(-log10(res$pvalue),na.rm=T),max(-log10(res$pvalue),na.rm=T)))
(logpadjLim=c( min(-log10(res$padj),na.rm=T),max(-log10(res$padj),na.rm=T)))

cex<-0.3
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj>=.05 & abs(log2FoldChange)<=1), 
	plot(log2FoldChange, -log10(pvalue), pch=20, col="black",cex=cex,
		main=paste("Volcano plot ",paste0(contrast,collapse="."), ", alpha=",alpha,sep=""),
		xlim=lfcLim, ylim=logpLim
	)
)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20,cex=cex, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20,cex=cex, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20,cex=1.5*cex, col="green"))

dev.off()

#----------------------------------------------------------------------
# scatter plot
#----------------------------------------------------------------------


# get the expression levels
geneData <- data.frame(x=assay(rldc)[,contrast[2]], y=assay(rldc)[,contrast[3]])

