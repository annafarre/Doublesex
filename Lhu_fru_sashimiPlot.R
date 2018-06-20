## now draw an example - Loci 10014
# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")

library("GenomicRanges")
library("Gviz")

options(ucscChromosomeNames=FALSE)
setwd("/Users/afarre/Unimelb/Doublesex/Lhumile/")


# make granges for each of the genes
gene_lengths=read.delim("SuperDuper.fasta.fai",stringsAsFactors=F,header=F)[,1:2]
n=dim(gene_lengths)[1]
gene_gr=GRanges(gene_lengths[,1],IRanges(rep(1,n),gene_lengths[,2]),rep("*",n))

#fruitless
gene="gene3721:MSTRG.4113"



# afrom <- -300 #1
# ato <- gene_lengths$V2[gene_lengths$V1==gene]

#genomic gap
novel_pos=gene_gr[seqnames(gene_gr)==gene]

#alignment track
alTrack2 <- AlignmentsTrack("fru2ST.bam",type = c("sashimi"),sashimiHeight=100,name=" ",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=1)

#to index bam files 
library(Rsamtools)
indexBam("fru2ST.bam")
indexBam("fru2ST_queensPRJDB4088.bam")
indexBam("fru2ST_workersPRJDB4088.bam")
indexBam("fru2ST_workersPRJNA279338.bam")


alTrackW <- AlignmentsTrack("fru2ST_workersPRJDB4088.bam",type = c("coverage"),name="Workers",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)
alTrackW2 <- AlignmentsTrack("fru2ST_workersPRJNA279338.bam",type = c("coverage"),name="Queens Other",
                             background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)
alTrackQ <- AlignmentsTrack("fru2ST_queensPRJDB4088.bam",type = c("coverage"),name="Queens",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)


#add in the transcripts
library(rtracklayer)
transGranges<-as(import.gff3("Lhu_fru_exons.gff"), "GRanges")
elementMetadata(transGranges)$transcript=gsub("_mid1","",mcols(transGranges)$Parent)
# elementMetadata(transGranges)$transcript=gsub("MSTRG.15287.1","DsxF",mcols(transGranges)$transcript)
# elementMetadata(transGranges)$transcript=gsub("rna24360","DsxM",mcols(transGranges)$transcript)

col<-c("red")
transTrack <- GeneRegionTrack(transGranges,fill=col,name="Transcripts",background.title ="#619CFF", size=3)


#queens, males and workers separately
ht2 <- HighlightTrack(trackList = list(alTrackW,alTrackQ,alTrackW2,alTrack2,transTrack),#alTrack2,transTrack),
                      novel_pos,alpha=0,inBackground=F)
ht2 <- HighlightTrack(trackList = list(alTrackW,alTrackQ),#alTrack2,transTrack),
                      novel_pos,alpha=0,inBackground=F)
pdf(file="Lhu_fru_sashimiQW.pdf")
plotTracks(c(GenomeAxisTrack(),ht2), chromosome=gene,
           transcriptAnnotation = "transcript")
dev.off()


