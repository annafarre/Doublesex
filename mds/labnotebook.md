[TOC]

####17/04/2018

Using Wau_necklace_limma.rmd script
Filter set to cpm>0 as cpm>1 was too strict

Use the script for Apis mellifera

####18/04/2018 Find Dsx and Fem in _Apis mellifera_ ST

Find Dsx and Fem in _Apis mellifera_ SuperTranscripts (ST) produced by necklace. 

Since we used an annotation that already contained the genes we can just grep the names with 

```
grep "Dsx" genome_superTranscriptome/genome_merged.gft
grep "Fem" genome_superTranscriptome/genome_merged.gft 
```

Result: 

```
$ grep "Fem" genome_superTranscriptome/genome_merged.gft 
NC_007072.3	StringTie	transcript	11113282	11123371	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11113282	11113303	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "1"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11118364	11118687	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "2"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11119008	11119192	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "3"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11120905	11121015	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "4"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11121090	11121186	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "5"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11122040	11122087	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "6"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11122178	11122321	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "7"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11122741	11122888	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "8"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11122952	11123099	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "9"; gene_name "Fem"; ref_gene_id "gene3145"; 
NC_007072.3	StringTie	exon	11123167	11123371	1000	+	.	gene_id "MSTRG.3247"; transcript_id "rna7001"; exon_number "10"; gene_name "Fem"; ref_gene_id "gene3145"; 
```

Then we can grep the StringTie gene name in the SuperDuper file so we get the complete ST name of the gene:

```
$ grep "MSTRG.3247" superTranscriptome/SuperDuper.gff 
gene3145:MSTRG.3247	SuperTranscript	exon	1	2766	.	.	0	.
```

####23/04/2018

Created and md to report on the *Vollenhovia emeryi* issue: excess of splice sites reported by HISAT2 which generates an excessive number of blocks for ST. 

~/Unimelb/Doublesex/mds/Vem_issue_excessBlocks.md 

####24/04/2018

R analysis for Ame, Aec and Wau. Instersections of genes showing alternative splicing. 

Rmd scripts. 

Analysis Monomorium pharaonis data.

Intersections done using UpsetR in: Intersections_AS.Rmd

![Intersections](/Users/afarre/Unimelb/Doublesex/figures/Intersection_QW_ame_aec_wau_mph_ame4day.pdf)

####25/04/2018

Merged bam files of last necklace mapping (to ST) and used StringTie on it. Extracted splice junctions and made blocks. Reduced greatly the number of blocks:

**STATS**

Species |#StringTie exons|#ST blocks|#Stringtie exons|#Stringtie blocks
:--|:--|:--|:--|:--
V.emeryi|277,505|999,514|80,807|71,446
W.auropunctata|254,989|172,634|62,330|49,058
A.mellifera|291,739|247,922|NA|NA
M.pharaonis|260,650|159,248|52,250|NA
Aechinatior|233,719|177,313|44,438|NA

In wasmannia auropunctata only one exon of dsx was preserved after the extra steps. Probably the extra step are filtering too much. Droped as option. 

####02/05/2018
Got list of proteins involved in sex determination in drosophila. Trying to find the ortholgs in Wasmannia auropunctata. 

List og genes involved in Sex det in Drosophila melanogaster: 

`Sex_det_prot_Dme_list.txt`

```
transformer 2
transformer
stand still
sisterless A
female lethal d
daughterless
extra macrochaetae
SRm160
doublesex
YTH domain containing 1
runt
ovo
ovarian tumor
hopscotch
PHD finger protein 7
fruitless
intersex
sans fille
hermaphrodite
groucho
Signal-transducer and activator of transcription protein at 92E
Methyltransferase like 14
scute
Sex lethal
Chronologically inappropriate morphogenesis
Methyltransferase like 3
virilizer
degringolade
unpaired 1
deadpan
Nuclear receptor binding SET domain protein
dissatisfaction
```


[This paper](https://www.cell.com/action/showImagesData?pii=S1534-5807%2814%2900738-2) identified doublesex target genes. I downloaded the S1 table which contains the results of their occupancy testing. Filtered the genes to get only the ones belonging to Peaksum kmeans cluster number 5. These are the ones with highest occupancy and the ones that show tissue nonspecific expression. 




