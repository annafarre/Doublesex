#!/bin/bash
 
#SBATCH -p physical
#SBATCH --time=09:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --job-name=StringTie_${filename}
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=stringtie_%j.out
#SBATCH --err=stringtie_%j.err
 
module load StringTie
 
stringtie -v -o stringtie_merged_out_bam.gtf merged_out.bam -p 3