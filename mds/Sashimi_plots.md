#How to make sashimi plots
[TOC]
##Identify the gene of interest in the genome of the target species
Use Blat or BLAST to determine the orthologue of the gene of interest (e.g. doublesex). 

Fruitless: 

```
grep "gene7455:MSTRG.7359" cluster_files/relST_genomeST.psl
```



##Extract list SRA ID numbers for queens and workers
Will use the lists of queens, males and workers to merge the bam files later.

```
$ grep "worker" ~/Wauropunctata/data_info/Wau_females_SraRunTable.txt | grep -o "DRR[[:digit:]]*" > ~/Wauropunctata/data_info/Wau_workers_SRR_Acc_List.txt
$ grep "queen" ~/Wauropunctata/data_info/Wau_females_SraRunTable.txt | grep -o "DRR[[:digit:]]*" > ~/Wauropunctata/data_info/Wau_queens_SRR_Acc_List.txt
```
##Extract the gene of interest from the bam files

For each bam file, subsets it to only the gene of interest.

```
module load SAMtools 
for i in `ls DRR*.bam`; do  samtools view -h $i gene13749:MSTRG.15287 >> dsx_$i ; done
```
##Merge bam files 

Use the list of queens, males and workers to merge bam files by caste.

Merge also all bams together.

Can run it in command line: 

```
samtools merge dsx2ST.bam sample1.bam sample2.bam [...]
```
Easy way to list all samples per caste: 

```
for i in `cat ~/Wauropunctata/data_info/Wau_queens_SRR_Acc_List.txt`; do printf "fru_$i.bam " ; done 
```

Or as a script:

```
#!/bin/bash


#SBATCH -p physical-cx4
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=mergebam
#SBATCH --mem=20GB
#SBATCH --out=slurm_%j.out  
#SBATCH --err=slurm_%j.err  
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --mail-type=ALL

module load SAMtools

samtools merge dsx2ST.bam $@
```

Output: 

1. All bams merged 
2. Queen merged bams
3. Worker merged bams 
4. Male merged bams  

##Map the transcripts to the SuperTranscript

The SuperTranscript does not give us information about isoforms.

To map the transcripts to the ST, first subset the genome_merged.gtf file (Necklace output) to get only the gene of interest. 

```
grep "LOC106744688" genome_superTranscriptome/genome_merged.gft > Dqu_dsx_exons.gtf
```


Then, using the gtf of the gene, extract the fasta sequence: 

```
module load Cufflinks
gffread Wau_dsx_exons.gtf -g ~/Wauropunctata/genomes/GCF_000956235.1_wasmannia.A_1.0_genomic.fna -w Wau_dsx_exons.fasta
```

Use the fasta sequence of the gene of interest in Blat againt the SuperTranscriptome.

BLAT v. 35

```
module load BLAT
blat superTranscriptome/SuperDuper.fasta Wau_dsx_exons.fasta -minScore=200 -minIdentity=98 Wau_dsx_exons.psl
```
Use the `blat2gff.pl` script to convert the psl output from blat to a gff file.

```
~/scripts/blat2gff.pl < Wau_dsx_exons.psl > Wau_dsx_exons.gff
```

##Make the sashimi plot