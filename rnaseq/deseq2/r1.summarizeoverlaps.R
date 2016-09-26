# ======================================================================
# setup sample metadata, BAMs and count reads (multi-core)
#
# sample_list.txt = sample metadata TSV (sample, condition, bam, ...)
# env: NSLOTS = available cores
# env: REF_GTF = path to transcript GTF
# ======================================================================

# pre-load needed libraries
library(GenomicAlignments)
library(BiocParallel)
library(Rsamtools)
library(GenomicFeatures)

#
# load sample list
#
# required columns:
#   sample: unique sample name
#   bam:    path to bam file
#
(sampleTable <- read.table("sample_list.txt",header=T))
sampleTable$bam <- as.character(sampleTable$bam) # un-factor file names

# if you want to order factors, you can do that now
#
#sampleTable$cell <- factor(sampleTable$cell <- subset,c("naive","unstim","minus", "plus")) # not syntax correct


# check if files exist
if( length(sampleTable$bam) == sum(file.exists(sampleTable$bam)) ) {
	print(paste("all bams, pos   sort: ", sum(file.exists(sampleTable$bam)),    " present"))
} else {
	print(paste("some bams, pos   sort: ", sum(FALSE==(file.exists(sampleTable$bam))),    " missing"))
	sampleTable$bam[!file.exists(sampleTable$bam)]
	stop("missing bam files")
}

	
#
# make "bam" file reader
#
library("Rsamtools")
# create BAM list
# asMates:  (Bioconductor > 2.12) use asMates=T for paired-end, even if sorted by position
# bieldSize: limits how man reads are loaded into memory at a one time (batch size), limit RAM usage.
Rbamfiles  <- BamFileList(sampleTable$bam, yieldSize=2000000, asMates=TRUE)
names(bamfiles) <- sampleTable$sample_name

# print out genome for manual validation
seqinfo(bamfiles[1])

#--------------------------------------------------
# read gene annotations
#--------------------------------------------------
library(GenomicFeatures)
(ref_gtf <- Sys.getenv("REF_GTF",unset="../../ucsc.hg19.iGenomes.genes.plusChrM.plusFIREFLYluc_txid.gtf"))

(txdb <- makeTxDbFromGFF(ref_gtf, format="gtf", circ_seqs=character()))
#(txdb <- makeTxDbFromGFF(ref_gtf, format="gtf", circ_seqs=character(), organism="Mus musculus (mm10.igenomes.ucsc.2015)", taxonomyId="10090"))

# build gene/tx->exon map
ebt <- exonsBy(txdb, by="tx", use.names=TRUE)
(ebg <- exonsBy(txdb, by="gene"))

#----------------------------------------------------------------------
# compute read counts
#
# http://www.bioconductor.org/help/workflows/rnaseqGene/#count
#
# cheaha: 30 minutes for all 22 Zajac samples on 15 CPU @ 2 G RAM each
#---------------------------------------------------------------------- 

# several cores on one machine, if available
nCores <-strtoi(Sys.getenv("NSLOTS",unset=1))
nCores <-strtoi(Sys.getenv("SLURM_CORES",unset=nCores))
register(MulticoreParam(workers = nCores), default = TRUE)
registered()[1]

print(paste("start:", date()))
# compute GENE values
se.gene <- summarizeOverlaps(features=ebg, 
			reads=bamfiles,
                        mode="Union",  # perhaps "IntersectionStrict" is better - esp for transcripts/isoforms!
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
# add metadata
colData(se.gene) <- DataFrame(sampleTable)
save(file="se.gene.union.RData",se.gene)
print(paste("done:", date()))

print(paste("start:", date()))
# compute TRANSCRIPT values
se.tx <- summarizeOverlaps(features=ebt, 
			reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
# add metadata
colData(se.tx) <- DataFrame(sampleTable)
save(file="se.tx.union.RData",se.tx)
print(paste("done:", date()))

#----------------------------------------------------------------------
# save and exit
#---------------------------------------------------------------------- 
save.image()
