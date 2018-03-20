#!/bin/bash


#SBATCH -p cloud
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=align
#SBATCH --mem=5GB
#SBATCH --out=slurm_%j.out  
#SBATCH --err=slurm_%j.err  
#SBATCH --mem-per-cpu=5000
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --mail-type=ALL


#### READ ME
#align proteins using MUSCLE
#pass 2 arguments through STDIN 
# $1 filein with prots 
# $2 output, alignment 

module load MUSCLE
muscle  -in "$1" -phyiout "$2"
