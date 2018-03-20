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
#create a tree with FastTree from (prot or others) alignment
#tree made with WAG gamma
# $1 filein with alignment
# $2 output, tree 


module load FastTree
FastTree -wag -gamma < "$1"  > "$2"
