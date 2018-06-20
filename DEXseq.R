########################################################################
################ Author: Anthony Hawkins, Nadia Davidson ###############
# Description: A script to extract differential exon usage from 
# a counts matrix using DEXseq. 
########################################################################

library(dplyr)
library(DEXSeq)


#DexSeq requires a table with information about the samples as well as the raw counts for each exon in each sample
#And grouping information for which exons belong to which genes/clusters
sample_id<-read.table("/Users/afarre/Unimelb/Doublesex/Mpharaonis/data_info/Mph_pe_SRR_Acc_List.txt", col.names = "sampleID", stringsAsFactors = F )
sample_id 

info_samples<-read.table("/Users/afarre/Unimelb/Doublesex/Mpharaonis/data_info/Mph_pe_SraRunTable.txt", sep = "\t", header = T,  stringsAsFactors = F )

info_samples <- info_samples[info_samples$BioProject_s=="PRJDB4088",]
info_samples

#Get the sample names per caste 
Queens <- info_samples[grep(".*Q\\d$", info_samples$sample_name_s, perl = T),]$Run_s
Workers <- c(info_samples[grep(".*W\\d$", info_samples$sample_name_s, perl = T),]$Run_s, info_samples[grep("^PHARAO\\..*COL\\d$", info_samples$sample_name_s, perl = T),]$Run_s)

DRP002877 <- info_samples[grep(".*[Q|W]\\d$", info_samples$sample_name_s, perl = T),]$Run_s

#add a population and caste column
info_samples$casteName <- ""
info_samples$caste <- ""

info_samples[info_samples$Run_s %in% Queens,]$casteName <- "Queen"
info_samples[info_samples$Run_s %in% Workers,]$casteName <- "Worker"

info_samples[info_samples$Run_s %in% Queens,]$caste <- "1"
info_samples[info_samples$Run_s %in% Workers,]$caste <- "2"


data_info <- dplyr::select(info_samples, Run_s, casteName)
row.names(data_info) <- data_info$Run_s
data_info$casteName <- as.factor(data_info$casteName)
colnames(data_info) <- c("sample", "condition")

#First get counts table
counts <- read.table("/Users/afarre/Unimelb/Doublesex/Mpharaonis/necklace_ame/block.counts",sep="\t",header=TRUE,stringsAsFactors=F)
names(counts) <- gsub(".*(DRR.*)\\.bam","\\1", names(counts))




#each exon needs a unique name
#create a column with unique exon ids 
#exonid = geneid + exon number
counts <- do.call(rbind, lapply(unique(counts$Geneid), function(x){ 
  mutate(dplyr::filter(counts,Geneid==x), exon_ids=
           paste(dplyr::filter(counts,Geneid==x)$Geneid, 
                 rownames(dplyr::filter(counts,Geneid==x)), sep = "."))
} ))

row.names(counts) <- counts$exon_ids


#create the matrix with only the counts
count_matrix <- data.matrix(select(counts, -Geneid, -Chr, -Start, -End, -Strand, -Length, -exon_ids))
count_matrix <- round(count_matrix)
gene_ids <- counts$Geneid
exon_ids <- counts$exon_ids



#design

dxd <- DEXSeqDataSet(count_matrix,design=~sample + exon + condition:exon, featureID=as.factor(exon_ids),
                     groupID=as.factor(gene_ids),sampleData=data_info)
##Estimate size factors and dispersions (DEXseq does this based on a negative bionmial distribution
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)


#Test for DEU
dxd <- testForDEU(dxd)

#Extract the results
res <- DEXSeqResults(dxd)

#Get p-value per gene based on the p-values for all exons in the gene/cluster
pgq <- perGeneQValue(res, p = "pvalue")

## Save results to a text file and R object
save(dxd, res, pgq, file = "Mph_DEXseq.Rdata")
tmp <- cbind(gene = names(pgq), "adjP" = pgq)
write.table(tmp, file = "Mph_DEXseq_byGene.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(DEXSeqResults(dxd), file="Mph_DEXseq_byExon.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t" )


