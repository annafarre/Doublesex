# Finding gene orthologs in necklace output

Runing necklace with a non-model organism (e.g. _Wasmannia auropunctata_) and a model related species (e.g. _Drosophila melanogaster_)

[toc] 

##Necklace output 
```
$ pwd
/home/afarre/Wau_necklace
$ head counts/block.counts
# Program:featureCounts v1.5.3; Command:"/home/afarre/.local/necklace-necklace_v0.9/tools/bin/featureCounts" "-T" "32" "--primary" "-p" "-t" "exon" "-g" "gene_id" "--fraction" "-O" "-f" "-a" "counts/blocks.gtf" "-o" "counts/block.counts" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029066.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029067.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029068.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029069.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029070.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029071.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029072.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029073.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029074.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029075.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029076.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029077.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029078.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029079.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029080.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029081.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029082.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029083.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029084.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029085.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029086.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029087.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029088.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029089.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029090.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029091.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029092.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029093.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029094.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029095.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029096.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029097.bam" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029098.bam" 
Geneid	Chr	Start	End	Strand	Length	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029066.bam/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029067.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029068.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029069.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029070.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029071.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029072.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029073.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029074.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029075.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029076.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029077.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029078.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029079.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029080.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029081.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029082.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029083.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029084.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029085.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029086.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029087.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029088.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029089.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029090.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029091.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029092.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029093.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029094.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029095.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029096.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029097.bam	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029098.bam
MSTRG.1	MSTRG.1	1	315	+	315	0	0	0	0	0	0	0	0	2.80	0	0	0	0	1	0	0	0.20	0	0	0	0	0	48.27	50.32	119.25	237.67	254.82	16.95	43.68	67.35	40.63
MSTRG.10	MSTRG.10	1	303	+	303	0	0	0	0	0	0	1.50	39.50	24.50	14.50	0	0.25	17
MSTRG.100	MSTRG.100	1	316	+	316	1.75	7.08	4.75	4.83	3.75	9.25	2.92	9.25	3.75	7.75	3	2.75	2.25	4.33	4.42	8.50	4.25	2.33	1	10.17	6.33	8.92	13.83	8.75	27.92	32.83	41.00	27.25	17.25	6	17.33	21.50	11.83
MSTRG.10000	MSTRG.10000	1	177	+	177	0	1.75	1.75	1.67	2	1.25	0.75	3.75	1.75	1.50	2.67	0.50	2	5.08	3.75	1.67	8.17	0.83	5	4.75	2.50	0	2.92	3.49	7.50	8.50	16.75	6	6	2.25	13.50	10.25	4.92
MSTRG.10000	MSTRG.10000	178	735	+	558	0	1.25	0.50	0.17	0.50	0.25	0.75	1.25	2.75	0	0.67	0	2.50	0.75	2.25	0.67	3.17	0	0.25	0.25	0.50	0	0.92	4.42	34.83	59.83	15.25	5	9	3	19.25	14.75	6.42
MSTRG.10000	MSTRG.10000	736	1455	+	720	0	0	0	0	0	0	2.50	0	0	0	0	0.50	0	1	0	0	0	0	1	0	2.50	643.58	944.17	31	0	4	0.50	1	1	2
MSTRG.10001	MSTRG.10001	1	384	+	384	0	0	0	0	0	0	0.50	0	0	0	0	0	0	0	0	0	0	0	0	0	100.25	194	7.50	0	0	0	0	1	0
MSTRG.10008	MSTRG.10008	1	186	+	186	131.67	18	76	16.92	11.83	27.50	17.75	68.42	76.25	28	15.25	8	65.25	39.75	63.58	54.17	60.08	53.83	64.25	79.17	48.58	54.33	36.75	50.58	137.17	171.33	369.33	234.17	218.83	32.33	290.83	332.83	162.92
```

Gene names come from StringTie, which gives no information about the gene. 

Necklace uses blat to find orthologs of the target species in a given close relative. Output of the blat is found in `cluster_files`: 



```
$ pwd
/home/afarre/Wau_necklace/cluster_files
$ head relST_genomeST.psl
psLayout version 3

match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
---------------------------------------------------------------------------------------------------------------------------------------------------------------
289	20	0	0	1	35	1	36	-+	MSTRG.4	410	11	355	gene:FBgn0267501:gene:FBgn0267502:gene:FBgn0267503:gene:FBgn0267504:gene:FBgn0267507:MSTRG.16881	8118	1821	2166	172,137,	55,262,	1821,2029,
289	20	0	0	1	35	1	36	-+	MSTRG.4	410	11	355	gene:FBgn0267498:gene:FBgn0267499:gene:FBgn0267500:gene:FBgn0267506:MSTRG.16880	4191	1821	2166	2	172,137,55,262,	1821,2029,
243	17	0	0	1	35	1	36	-+	MSTRG.4	410	11	306	gene:FBgn0250731:gene:FBgn0267516:gene:FBgn0267521:gene:FBgn0267522:MSTRG.16768	5725	1769	2065	2	123,137,104,262,	1769,1928,
290	25	0	0	1	26	1	27	-+	MSTRG.4	410	11	352	gene:FBgn0085817:MSTRG.33	1266	367	709	2	178,137,	58,262,	367,572,
266	16	0	0	2	61	2	48	-+	MSTRG.4	410	11	354	gene:FBgn0085813:MSTRG.16761	1975	955	1285	3	54,91,137,	56,136,262,	955,1021,1148,
```
`relST_genomeST.psl` gives gene IDs found in the gtf file of the related species (_Drosophila melanogaster_ in this case). Use gtf to find gene names for each gene ID. 

```
$ head /home/afarre/original_genomes/Drosophila_melanogaster.BDGP6.37.gtf
211000022278158	FlyBase	exon	592	1036	.	+	.	transcript_id "transcript:FBtr0114187"; gene_id "gene:FBgn0085737"; gene_name "CR40502";
211000022278279	FlyBase	exon	11634	12709	.	-	.	transcript_id "transcript:FBtr0347033"; gene_id "gene:FBgn0267594"; gene_name "CR45932";
211000022278282	FlyBase	exon	518	2511	.	+	.	transcript_id "transcript:FBtr0114261"; gene_id "gene:FBgn0085807"; gene_name "CR41590";
211000022278298	FlyBase	exon	648	1614	.	+	.	transcript_id "transcript:FBtr0114196"; gene_id "gene:FBgn0085746"; gene_name "CR40571";
211000022278307	FlyBase	exon	1	824	.	-	.	transcript_id "transcript:FBtr0114201"; gene_id "gene:FBgn0085751"; gene_name "CR40582";
211000022278309	FlyBase	exon	1	1110	.	+	.	transcript_id "transcript:FBtr0114198"; gene_id "gene:FBgn0085748"; gene_name "CR40573";
211000022278436	FlyBase	exon	1551	2815	.	-	.	transcript_id "transcript:FBtr0300146"; gene_id "gene:FBgn0259849"; gene_name "Su(Ste):CR42418";
211000022278449	FlyBase	exon	528	1187	.	-	.	transcript_id "transcript:FBtr0308947"; gene_id "gene:FBgn0085494"; gene_name "Mst77Y-16Psi";
211000022278449	FlyBase	exon	1245	1635	.	-	.	transcript_id "transcript:FBtr0308947"; gene_id "gene:FBgn0085494"; gene_name "Mst77Y-16Psi";
211000022278498	FlyBase	exon	13189	13940	.	-	.	transcript_id "transcript:FBtr0114286"; gene_id "gene:FBgn0085828"; gene_name "CR42195";
```