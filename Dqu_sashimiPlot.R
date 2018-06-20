

## now draw an example - Loci 10014
# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")

library("GenomicRanges")
library("Gviz")

options(ucscChromosomeNames=FALSE)
setwd("/Users/afarre/Unimelb/Doublesex/Dquadriceps/necklace/")


# make granges for each of the genes
gene_lengths=read.delim("SuperDuper.fasta.fai",stringsAsFactors=F,header=F)[,1:2]
n=dim(gene_lengths)[1]
gene_gr=GRanges(gene_lengths[,1],IRanges(rep(1,n),gene_lengths[,2]),rep("*",n))


gene="gene5793:MSTRG.5916"

#genomic gap
novel_pos=gene_gr[seqnames(gene_gr)==gene]

#alignment track
alTrack2 <- AlignmentsTrack("dsx2ST.bam",type = c("sashimi"),sashimiHeight=100,name=" ",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=1)

#to index bam files 
library(Rsamtools)
indexBam("dsx2ST_reproductive.bam")
indexBam("dsx2ST_non-reproductive.bam")
indexBam("dsx2ST.bam")


alTrackQ <- AlignmentsTrack("dsx2ST_reproductive.bam",type = c("coverage"),name="Reproductive",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)
alTrackW <- AlignmentsTrack("dsx2ST_non-reproductive.bam",type = c("coverage"),name="Non-reproductive",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)


#add in the transcripts
library(rtracklayer)
transGranges<-as(import.gff3("Dqu_dsx_only_exons.gff"), "GRanges")
elementMetadata(transGranges)$transcript=gsub("_mid1","",mcols(transGranges)$Parent)
elementMetadata(transGranges)$transcript=gsub("MSTRG.15287.1","DsxF",mcols(transGranges)$transcript)
elementMetadata(transGranges)$transcript=gsub("rna24360","DsxM",mcols(transGranges)$transcript)

col<-c("blue", "red", "red")
transTrack <- GeneRegionTrack(transGranges,fill=col,name="Transcripts",background.title ="#619CFF", size=3)

#queens, males and workers separately
ht2 <- HighlightTrack(trackList = list(alTrackQ,alTrackW,alTrack2,transTrack),#alTrack2,transTrack),
                      novel_pos,alpha=0,inBackground=F)

pdf(file="Dqu_sashimiQW.pdf")
plotTracks(c(GenomeAxisTrack(),ht2), chromosome=gene,
           transcriptAnnotation = "transcript")
dev.off()





