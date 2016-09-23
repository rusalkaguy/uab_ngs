#!/bin/bash 
# cheaha1
module load ngs-ccts/bedtools
module load ngs-ccts/ucsc_kent/2014-03-05
genome=hg19
chromSizes=${genome}.chrom.sizes
if [ ! -e "${chromSizes} " ]; then echo fetching ${chromSizes}; fetchChromSizes  ${genome} > ${chromSizes}; fi
bamlist="$*"
if [ -z "${bamlist}" ]; then bamlist=$(ls -1 *.bam); fi
echo BAMLIST=${bamlist}
for bam in ${bamlist}; do
	# takes 10 minutes to run? 
	# bedGraph format: https://genome.ucsc.edu/goldenPath/help/bedgraph.html
	# RNA-Seq use case http://www.cureffi.org/2013/11/18/an-mrna-seq-pipeline-using-gsnap-samtools-cufflinks-and-bedtools/ 
	# genomeCoverageBed docs: http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html 
	# -ibam = input bam file
	# -bg = Report depth in BedGraph format
	# -split =  parse BAM CIGAR strings for alignment gaps and don.t count the gaps (introns)
	# -pc = for paired-end reads, coverage is calculated as the number of fragments covering each base pair
	echo computing ${bam}.bedgraph
	genomeCoverageBed -ibam ${bam} -bg -split  > ${bam}.bedgraph

	# covert to big wig
	# https://genome.ucsc.edu/goldenpath/help/bigWig.html#Ex3
	echo converting to ${bam}.bw
	bedGraphToBigWig ${bam}.bedgraph ${chromSizes} ${bam}.bw 

	# cleanup
	rm ${bam}.bedgraph
done
