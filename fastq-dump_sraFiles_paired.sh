#!/bin/bash

#for PAIRED samples
#create scripts that will fastq-dump sra files
#give through STDIN a list of sra files (give the path to the list or run in same directory)
#either:
#	the scripts have to be in the same directory as the sra files 
#	path to the sra samples is specified through STDIN $2

#STDIN


while read line
do 
filename=${pwd}${line:0:9}_fastq-dump.sh
cat <<-EOF > ${filename}
	#!/bin/bash 
 	#SBATCH -p cloud 
 	#SBATCH --time=24 
	#SBATCH --nodes=1 
	#SBATCH --ntasks=1 
	#SBATCH --cpus-per-task=1 
	#SBATCH --job-name=fastq-dump 
	#SBATCH --mem-per-cpu=50000 
	#SBATCH --mail-type=ALL 
	#SBATCH --mail-user=afarre@student.unimelb.edu.au 
	#SBATCH --out=slurm_%j.out 
	#SBATCH --err=slurm_%j.err 

	fastq-dump --split3 --readids --skip-technical --clip --read-filter pass --dumpbase $2${line}.sra 
EOF
done <$1
