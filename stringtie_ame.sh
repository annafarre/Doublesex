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

echo "#!/bin/bash" >${pathScript}
echo " " >>${pathScript}
echo "#SBATCH -p physical" >>${pathScript}
echo "#SBATCH --time=09:00:00" >>${pathScript}
echo "#SBATCH --nodes=1" >>${pathScript}
echo "#SBATCH --ntasks=1" >>${pathScript}
echo "#SBATCH --cpus-per-task=1" >>${pathScript}
echo "#SBATCH --job-name=StringTie_${filename}" >>${pathScript}
echo "#SBATCH --mem-per-cpu=4000" >>${pathScript}
echo "#SBATCH --mail-type=ALL" >>${pathScript}
echo "#SBATCH --mail-user=afarre@student.unimelb.edu.au" >>${pathScript}
echo "#SBATCH --out=slurm_${filename}_%j.out" >>${pathScript}
echo "#SBATCH --err=slurm_${filename}_%j.err" >>${pathScript}
echo " " >>${pathScript}
echo "cd \"$(pwd)/stringtie_${filename}\"" >>${pathScript}
echo " " >>${pathScript}
echo "input=$(pwd)/${longfilename}" >>${pathScript}
echo " " >>${pathScript}
echo "module load StringTie" >>${pathScript}
echo " " >>${pathScript}
echo "stringtie -v -G /data/projects/punim0329/genome/GCF_000002195.4_Amel_4.5_genomic.gff  -o stringtie_${filename}.gtf \${input} -p 3" >>${pathScript}

done 