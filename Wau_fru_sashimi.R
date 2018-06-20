## now draw an example - Loci 10014
# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")

library("GenomicRanges")
library("Gviz")

options(ucscChromosomeNames=FALSE)
setwd("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/necklace_ame/sashimiPlot/")


# make granges for each of the genes
gene_lengths=read.delim("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/SuperDuper.fasta.fai",stringsAsFactors=F,header=F)[,1:2]
n=dim(gene_lengths)[1]
gene_gr=GRanges(gene_lengths[,1],IRanges(rep(1,n),gene_lengths[,2]),rep("*",n))

#fruitless
gene="gene839:MSTRG.967"



# afrom <- -300 #1
# ato <- gene_lengths$V2[gene_lengths$V1==gene]

#genomic gap
novel_pos=gene_gr[seqnames(gene_gr)==gene]

#alignment track
alTrack2 <- AlignmentsTrack("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST.bam",type = c("sashimi"),sashimiHeight=100,name=" ",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=1)

#to index bam files 
library(Rsamtools)
indexBam("fru2ST_males.bam")
indexBam("fru2ST_queens.bam")
indexBam("fru2ST_workers.bam")
indexBam("fru2ST.bam")


alTrackM <- AlignmentsTrack("fru2ST_males.bam",type = c("coverage"),name="Males",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)
alTrackW <- AlignmentsTrack("fru2ST_workers.bam",type = c("coverage"),name="Workers",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)
alTrackQ <- AlignmentsTrack("fru2ST_queens.bam",type = c("coverage"),name="Queens",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)


#add in the transcripts
library(rtracklayer)
transGranges<-as(import.gff3("Wau_fru_exons.gff"), "GRanges")
elementMetadata(transGranges)$transcript=gsub("_mid1","",mcols(transGranges)$Parent)
elementMetadata(transGranges)$transcript=gsub("MSTRG.15287.1","DsxF",mcols(transGranges)$transcript)
elementMetadata(transGranges)$transcript=gsub("rna24360","DsxM",mcols(transGranges)$transcript)

col<-c("red")
transTrack <- GeneRegionTrack(transGranges,fill=col,name="Transcripts",background.title ="#619CFF", size=3)


#queens, males and workers separately
ht2 <- HighlightTrack(trackList = list(alTrackM,alTrackQ,alTrackW,alTrack2,transTrack),#alTrack2,transTrack),
                      novel_pos,alpha=0,inBackground=F)

pdf(file="/Users/afarre/Unimelb/Doublesex/Wauropunctatat/necklace_ame/sashimiPlot/Wau_fru_sashimiQWM.pdf")
plotTracks(c(GenomeAxisTrack(),ht2), chromosome=gene,
           transcriptAnnotation = "transcript")
dev.off()


