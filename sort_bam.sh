#!/bin/bash

for longfilename in `ls $@`
do

filename=$(basename ${longfilename} .out.bam)

pathScript=$(pwd)/sort_bam_${filename}.sh

echo "#!/bin/bash" >${pathScript}
echo " " >>${pathScript}
echo "#SBATCH -p cloud" >>${pathScript}
echo "#SBATCH --time=09:00:00" >>${pathScript}
echo "#SBATCH --nodes=1" >>${pathScript}
echo "#SBATCH --ntasks=1" >>${pathScript}
echo "#SBATCH --cpus-per-task=3" >>${pathScript}
echo "#SBATCH --job-name=tsort_${filename}" >>${pathScript}
echo "#SBATCH --mem=10000" >>${pathScript}
echo "#SBATCH --mail-type=ALL" >>${pathScript}
echo "#SBATCH --mail-user=afarre@student.unimelb.edu.au" >>${pathScript}
echo "#SBATCH --out=slurm_%j.out" >>${pathScript}
echo "#SBATCH --err=slurm_%j.err" >>${pathScript}
echo " " >>${pathScript}
echo "input=${longfilename}" >>${pathScript}
echo "sorted=$(pwd)/${filename}_sorted.bam" >>${pathScript}
echo " " >>${pathScript}
echo "module load SAMtools" >>${pathScript}
echo " " >>${pathScript}
echo "samtools sort -@ 3 \${input} > \${sorted}" >>${pathScript}

done 