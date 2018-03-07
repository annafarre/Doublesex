#!/bin/bash

# $1 input phylip file
# $2 output tree (only name, no directory path!)

#SBATCH -p physical
#SBATCH --time=1-06
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=raxml
#SBATCH --mem=20GB
#SBATCH --out=slurm_%j.out  
#SBATCH --err=slurm_%j.err  
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --mail-type=ALL
 
module load RAxML/7.7.5-intel-2016.u3-mt-sse3
 
raxmlHPC-PTHREADS-SSE3 -s "$1" -n "$2" -m PROTGAMMAIWAG -T 8 -f a  -p 12345 -x 12345 -N 1000 
