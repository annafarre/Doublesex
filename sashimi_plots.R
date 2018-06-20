

## now draw an example - Loci 10014
# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")

library("GenomicRanges")
library("Gviz")

options(ucscChromosomeNames=FALSE)



# make granges for each of the genes
gene_lengths=read.delim("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/SuperDuper.fasta.fai",stringsAsFactors=F,header=F)[,1:2]
n=dim(gene_lengths)[1]
gene_gr=GRanges(gene_lengths[,1],IRanges(rep(1,n),gene_lengths[,2]),rep("*",n))


gene="gene13749:MSTRG.15287"



# afrom <- -300 #1
# ato <- gene_lengths$V2[gene_lengths$V1==gene]

#genomic gap
novel_pos=gene_gr[seqnames(gene_gr)==gene]

#alignment track
alTrack1 <- AlignmentsTrack("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST.bam",type = c("coverage"),name="Alignments",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A")
alTrack2 <- AlignmentsTrack("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST.bam",type = c("sashimi"),sashimiHeight=100,name=" ",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=1)

#to index bam files 
library(Rsamtools)
indexBam("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST_male.bam")
indexBam("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST_worker.bam")
indexBam("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST_queen.bam")


alTrackM <- AlignmentsTrack("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST_male.bam",type = c("coverage"),name="Males",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)
alTrackW <- AlignmentsTrack("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST_worker.bam",type = c("coverage"),name="Workers",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)
alTrackQ <- AlignmentsTrack("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST_queen.bam",type = c("coverage"),name="Queens",
                            background.title ="#042E8A",col="#042E8A",fill="#042E8A", size=2)

alTrackMsh <- AlignmentsTrack("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST_male.bam",type = c("sashimi"),sashimiHeight=100,name=" ",
                              background.title ="#042E8A",col="#042E8A",fill="#042E8A")
alTrackWsh <- AlignmentsTrack("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST_worker.bam",type = c("sashimi"),sashimiHeight=100,name=" ",
                              background.title ="#042E8A",col="#042E8A",fill="#042E8A")
alTrackQsh <- AlignmentsTrack("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/dsx2ST_queen.bam",type = c("sashimi"),sashimiHeight=100,name=" ",
                              background.title ="#042E8A",col="#042E8A",fill="#042E8A")

# #human match
# human_pos=human_gr[seqnames(human_gr)==gene]
# humTrack <- AnnotationTrack(human_pos,fill="black",name="CDS",
#                             group="Conserved CDS",
#                             groupAnnotation = "group")

#add in the transcripts
library(rtracklayer)
transGranges<-as(import.gff3("/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/Wau_dsx_exons.gff"), "GRanges")
elementMetadata(transGranges)$transcript=gsub("_mid1","",mcols(transGranges)$Parent)
elementMetadata(transGranges)$transcript=gsub("MSTRG.15287.1","DsxF",mcols(transGranges)$transcript)
elementMetadata(transGranges)$transcript=gsub("rna24360","DsxM",mcols(transGranges)$transcript)

col<-c("blue", "red", "red")
transTrack <- GeneRegionTrack(transGranges,fill=col,name="Transcripts",background.title ="#619CFF", size=3)


#AnnotationTrack(starts=c(100,400),width=100,chromosome=gene)

#genomic gap highlight
#not_in_genome_pos=not_in_genome[seqnames(not_in_genome)==gene]
ht <- HighlightTrack(trackList = list(alTrack1,alTrack2,transTrack),#alTrack2,transTrack),
                     novel_pos,alpha=0.3,inBackground=F)

plotTracks(c(GenomeAxisTrack(),ht), chromosome=gene,
           transcriptAnnotation = "transcript")

#queens, males and workers separately
ht2 <- HighlightTrack(trackList = list(alTrackM,alTrackQ,alTrackW,alTrack2,transTrack),#alTrack2,transTrack),
                     novel_pos,alpha=0,inBackground=F)

pdf(file="/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/Wau_sashimiQWM.pdf")
plotTracks(c(GenomeAxisTrack(),ht2), chromosome=gene,
           transcriptAnnotation = "transcript")
dev.off()





#plot queens only 
#queens, males and workers separately
htQ <- HighlightTrack(trackList = list(alTrackQ,alTrackQsh,transTrack),#alTrack2,transTrack),
                      novel_pos,alpha=0.3,inBackground=F)

#pdf(file="/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/Wau_sashimiQWM.pdf")
plotTracks(c(GenomeAxisTrack(),htQ), chromosome=gene,
           transcriptAnnotation = "transcript")
#dev.off()

#plot males only 
#queens, males and workers separately
htM <- HighlightTrack(trackList = list(alTrackM,alTrackMsh,transTrack),#alTrack2,transTrack),
                      novel_pos,alpha=0.3,inBackground=F)

#pdf(file="/Users/afarre/Unimelb/Doublesex/Wauropunctatat/sashimiPlots/Wau_sashimiQWM.pdf")
plotTracks(c(GenomeAxisTrack(),htM), chromosome=gene,
           transcriptAnnotation = "transcript")
#dev.off()




# plotTracks(c(GenomeAxisTrack(),ht), from = afrom, 
#            to = ato,chromosome=gene,
#            transcriptAnnotation = "transcript",
#            sizes=c(0.25,0.8,0.4,0.1,1))
