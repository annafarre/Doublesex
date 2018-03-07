#!/bin/bash

#SBATCH -p physical
#SBATCH --time=09:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=featureCounts_${filename}
#SBATCH --mem=40000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_%j.out
#SBATCH --err=slurm_%j.err



featureCounts -J -p -t exon -g exon -T 8 -G /data/projects/punim0329/genome/GCF_000002195.4_Amel_4.5_genomic.fna \
-a /data/projects/punim0329/genome/GCF_000002195.4_Amel_4.5_genomic.gfff  \
 -o /data/projects/punim0329/Amel_counts.txt $@ \

