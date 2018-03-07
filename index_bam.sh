#!/bin/bash

for longfilename in `ls $@`
do

filename=$(basename ${longfilename} .bam)

pathScript=$(pwd)/index_bam_${filename}.sh

cat <<-EOF > ${pathScript}
#!/bin/bash

#SBATCH -p physical
#SBATCH --time=09:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=index_bam_${filename}
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_%j.out
#SBATCH --err=slurm_%j.err


input=${longfilename}

module load SAMtools
samtools index \${input} 

EOF
done 