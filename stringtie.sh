#!/bin/bash

#Create scripts to run StringTie
#Input: list of bam files from STAR (passed as arguments in the command line, with wildcards)
#The scripts created separate folder for each sample


for longfilename in `ls $@`
do

filename=$(basename ${longfilename} Aligned.sortedByCoord.out.bam)
prefix="STAR_"
filename=${filename#$prefix}
mkdir -p $(pwd)/stringtie_${filename}
pathScript=$(pwd)/stringtie_${filename}/stringtie_${filename}.sh

cat <<-EOF > ${pathScript}
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
#SBATCH --out=slurm_${filename}_%j.out
#SBATCH --err=slurm_${filename}_%j.err
 
cd \"$(pwd)/stringtie_${filename}\"
 
input=$(pwd)/${longfilename}
 
module load StringTie
 
stringtie -v -o stringtie_${filename}.gtf \${input} -p 3
EOF

done 