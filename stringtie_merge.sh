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

stringtie --merge -o merged_stringtie.gtf $@ 