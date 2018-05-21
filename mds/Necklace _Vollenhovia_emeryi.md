### Necklace _Vollenhovia emeryi_

Running Necklace with run_necklace.sh :

```
#!/bin/bash

#SBATCH -p physical
#SBATCH --time=3-12
#SBATCH --job-name=necklace
#SBATCH --cpus-per-task=10
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_%j.out  
#SBATCH --err=slurm_%j.err 

module load Java

cd /home/afarre/Vem_necklace/

MAX_JAVA_MEM=2g /home/afarre/.local/necklace-necklace_v0.9/tools/bin/bpipe run /home/afarre/.local/necklace-necklace_v0.9/necklace.groovy /home/afarre/Vem_necklace/Vem_necklace.config 
```

And Vem_necklace.config

```
// sequencing data
reads_R1="/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030152_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030153_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030154_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030155_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030156_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030157_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030158_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030159_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030160_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030161_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030162_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030163_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030164_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030165_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030166_pass_1_val_1.fastq.gz" 
reads_R2="/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030152_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030153_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030154_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030155_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030156_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030157_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030158_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030159_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030160_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030161_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030162_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030163_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030164_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030165_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030166_pass_2_val_2.fastq.gz"

//The reference genome and its annotation
annotation="/data/projects/punim0304/Vemeryi/genomes/GCA_000949405.1_V.emery_V1.0_genomic.gff"
genome="/data/projects/punim0304/Vemeryi/genomes/GCA_000949405.1_V.emery_V1.0_genomic.fna"

//The genome and annotation of a related species
annotation_related_species="/home/afarre/original_genomes/Drosophila_melanogaster.BDGP6.37.gtf"
genome_related_species="/home/afarre/original_genomes/Drosophila_melanogaster.BDGP6.dna.toplevel.fa"
```

Got the error message

```
$ bpipe errors

============================== Found 1 failed commands from run 2840 ===============================

============================================ Command 59 ============================================

Command    : cat /data/projects/punim0304/Vemeryi/genomes/GCA_000949405.1_V.emery_V1.0_genomic.gff > genome_superTranscriptome/ref_annotations_combined.gtf ; /home/afarre/.local/necklace-necklace_v0.9/tools/bin/stringtie --merge -G genome_superTranscriptome/ref_annotations_combined.gtf -o genome_superTranscriptome/genome_merged.gft genome_guided_assembly/genome_assembly.gtf /data/projects/punim0304/Vemeryi/genomes/GCA_000949405.1_V.emery_V1.0_genomic.gff
Exit Code  : 1
Output     : 

No command output found in most recent log file
```

###SOLUTION: 

Genome and annotation used were GenBank
 
Use RefSeq instead <- contains predicted exons, transcripts, etc. 

New Vem_necklace.config:

```
// sequencing data
reads_R1="/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030152_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030153_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030154_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030155_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030156_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030157_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030158_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030159_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030160_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030161_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030162_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030163_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030164_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030165_pass_1_val_1.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030166_pass_1_val_1.fastq.gz" 
reads_R2="/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030152_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030153_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030154_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030155_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030156_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030157_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030158_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030159_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030160_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030161_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030162_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030163_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030164_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030165_pass_2_val_2.fastq.gz,/data/projects/punim0304/Vemeryi/SRA/trimGalore/necklaceInput/DRR030166_pass_2_val_2.fastq.gz"

//The reference genome and its annotation
annotation="/data/projects/punim0304/Vemeryi/genomes/GCF_000949405.1_V.emery_V1.0_genomic.gff"
genome="/data/projects/punim0304/Vemeryi/genomes/GCF_000949405.1_V.emery_V1.0_genomic.fna"

//The genome and annotation of a related species
annotation_related_species="/home/afarre/original_genomes/Drosophila_melanogaster.BDGP6.37.gtf"
genome_related_species="/home/afarre/original_genomes/Drosophila_melanogaster.BDGP6.dna.toplevel.fa"
```
