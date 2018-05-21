### Necklace _Wasmania auropunctata_ 

After running Necklace with the code: 

```
#!/bin/bash

#SBATCH -p physical
#SBATCH --time=3-12
#SBATCH --job-name=necklace
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_%j.out  
#SBATCH --err=slurm_%j.err 

module load Java

cd /home/afarre/Wau_necklace/

/home/afarre/.local/necklace-necklace_v0.9/tools/bin/bpipe run /home/afarre/.local/necklace-necklace_v0.9/necklace.groovy /home/afarre/Wau_necklace/Wau_necklace_trimmed.config 
```

and Wau\_necklace_config: 

```
// sequencing data
reads_R1="/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029066_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029067_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029068_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029069_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029070_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029071_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029072_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029073_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029074_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029075_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029076_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029077_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029078_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029079_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029080_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029081_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029082_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029083_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029084_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029085_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029086_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029087_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029088_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029089_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029090_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029091_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029092_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029093_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029094_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029095_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029096_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029097_pass_1_val_1.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029098_pass_1_val_1.fastq.gz" 
reads_R2="/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029066_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029067_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029068_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029069_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029070_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029071_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029072_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029073_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029074_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029075_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029076_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029077_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029078_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029079_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029080_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029081_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029082_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029083_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029084_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029085_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029086_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029087_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029088_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029089_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029090_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029091_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029092_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029093_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029094_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029095_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029096_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029097_pass_2_val_2.fastq.gz,/data/projects/punim0304/Wauropunctata/SRA/trimGalore/necklaceInput/DRR029098_pass_2_val_2.fastq.gz"

//The reference genome and its annotation
annotation="/data/projects/punim0304/Wauropunctata/genomes/GCF_000956235.1_wasmannia.A_1.0_genomic.gtf"
genome="/data/projects/punim0304/Wauropunctata/genomes/GCF_000956235.1_wasmannia.A_1.0_genomic.fna"

//The genome and annotation of a related species
annotation_related_species="/home/afarre/original_genomes/Drosophila_melanogaster.BDGP6.37.gtf"
genome_related_species="/home/afarre/original_genomes/Drosophila_melanogaster.BDGP6.dna.toplevel.fa"
```
Created a BLAST database with the output file SuperDuper.fasta 

And looked for doublesex in the fasta file. The result was saved in Wau_dsx_superTranscripts_pos.txt

```
XM_011706995.1	gene13749:MSTRG.15287	100.000	1351	0	0	1	1351	660	2010	0.0	2495
XM_011706995.1	gene13749:MSTRG.15287	100.000	601	0	0	1351	1951	8469	9069	0.0	1110

```

Looked for the blocks created by Necklace in blocks.gtf

Saved it in Wau\_doublesex_blocks.txt

```
gene13749:MSTRG.15287	superTranscript_blocks	exon	1	2008	.	+	.	gene_id "gene13749:MSTRG.15287"; transcript_id "gene13749:MSTRG.15287"
gene13749:MSTRG.15287	superTranscript_blocks	exon	2009	2009	.	+	.	gene_id "gene13749:MSTRG.15287"; transcript_id "gene13749:MSTRG.15287"
gene13749:MSTRG.15287	superTranscript_blocks	exon	2010	5366	.	+	.	gene_id "gene13749:MSTRG.15287"; transcript_id "gene13749:MSTRG.15287"
gene13749:MSTRG.15287	superTranscript_blocks	exon	5367	5634	.	+	.	gene_id "gene13749:MSTRG.15287"; transcript_id "gene13749:MSTRG.15287"
gene13749:MSTRG.15287	superTranscript_blocks	exon	5635	8468	.	+	.	gene_id "gene13749:MSTRG.15287"; transcript_id "gene13749:MSTRG.15287"
gene13749:MSTRG.15287	superTranscript_blocks	exon	8469	8826	.	+	.	gene_id "gene13749:MSTRG.15287"; transcript_id "gene13749:MSTRG.15287"
gene13749:MSTRG.15287	superTranscript_blocks	exon	8827	9069	.	+	.	gene_id "gene13749:MSTRG.15287"; transcript_id "gene13749:MSTRG.15287"

``` 

