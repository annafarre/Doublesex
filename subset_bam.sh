#!/bin/bash

for longfilename in `ls $@`
do

filename=$(basename ${longfilename} .out.bam)

pathScript=$(pwd)/subset_bam_${filename}.sh
#dependon=$(seq -s : 921301 921424)

echo "#!/bin/bash" >${pathScript}

echo "#SBATCH -p cloud" >>${pathScript}
echo "#SBATCH --time=09:00:00" >>${pathScript}
echo "#SBATCH --nodes=1" >>${pathScript}
echo "#SBATCH --ntasks=1" >>${pathScript}
echo "#SBATCH --cpus-per-task=1" >>${pathScript}
echo "#SBATCH --job-name=subset_bam_${filename}" >>${pathScript}
#echo "#SBATCH --dependency=afterok:${dependon}" >>${pathScript}
echo "#SBATCH --mem-per-cpu=10000" >>${pathScript}
echo "#SBATCH --mail-type=ALL" >>${pathScript}
echo "#SBATCH --mail-user=afarre@student.unimelb.edu.au" >>${pathScript}
echo "#SBATCH --out=slurm_%j.out" >>${pathScript}
echo "#SBATCH --err=slurm_%j.err" >>${pathScript}


echo "chrom=DF955771.1" >>${pathScript}
echo "input=${longfilename}" >>${pathScript}
echo "output=$(pwd)/${filename}_\${chrom}.bam" >>${pathScript}

echo "module load SAMtools" >>${pathScript}
echo "samtools view -h \${input} \${chrom} >> \${output}" >>${pathScript}

done 