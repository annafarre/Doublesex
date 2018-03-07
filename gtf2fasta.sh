#!/bin/bash

#take a gtf and genome ref and extract transcript sequences 
#input: gtf and ref genome fasta
#output: fasta


#SBATCH -p physical
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=gtf2fasta
#SBATCH --mem=20GB
#SBATCH --mail-type=ALL
#SBATCH --out=slurm_%j.out
#SBATCH --err=slurm_%j.err
#SBATCH --mail-user=afarre@student.unimelb.edu.au

module load TopHat

gtf_to_fasta $1 $2 $3