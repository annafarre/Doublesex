#!/bin/bash

#create an index for kallisto
#input: $1 index name $2 fasta file (transcripts)


#SBATCH -p physical
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=kallisto_idx
#SBATCH --mem=20GB
#SBATCH --mail-type=ALL
#SBATCH --out=slurm_%j.out
#SBATCH --err=slurm_%j.err
#SBATCH --mail-user=afarre@student.unimelb.edu.au


kallisto index --index=$1 $2