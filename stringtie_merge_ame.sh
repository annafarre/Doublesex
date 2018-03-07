#!/bin/bash

#merge all gtfs create with stringtie

#SBATCH -p physical
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=stringtie_merge
#SBATCH --mem=10000
#SBATCH --out=slurm_%j.out  
#SBATCH --err=slurm_%j.err  
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --mail-type=ALL

module load StringTie

stringtie --merge -G /data/projects/punim0329/genome/GCF_000002195.4_Amel_4.5_genomic.gff -o merged_stringtie.gtf $@ 