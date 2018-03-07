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

/home/afarre/.local/necklace-necklace_v0.9/tools/bin/bpipe run /home/afarre/.local/necklace-necklace_v0.9/necklace.groovy /home/afarre/Vem_necklace/Vem_necklace.config 