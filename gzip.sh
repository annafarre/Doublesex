#!/bin/bash

# Gunzip all files listed
# Input through STDIN 
# Use wildcards to zip multiple files
# Creates a script per file

for filename in `ls $@`
do

pathScript=$(pwd)/gzip_${filename}.sh

cat <<-EOF > ${pathScript}
#!/bin/bash
 
#SBATCH -p physical
#SBATCH --time=9-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=gzip_${filename}
#SBATCH --mem=25GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_%j.out
#SBATCH --err=slurm_%j.err
 
gzip ${filename}

EOF
done 
