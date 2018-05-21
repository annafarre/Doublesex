#!/bin/bash


#SBATCH -p cloud
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

samtools merge STAR_adultandlarvae.sortedByCoord.out.bam $@