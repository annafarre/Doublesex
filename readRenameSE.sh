#!/bin/bash

#Rename reads in fastq.gz files coming from SRA 
#to be used in Necklace 

#Example: 
	#Old name: DRRXXX.YYYY.Z 
	#New name: DRRXXX_YYYY/Z 
	
	#DRRXXX SRA sample name/number
	#YYYY read number
	#Z read pair number - 1 or 2


for filename in `ls $@`
do
cat <<-EOF > ${filename}_readRename.sh
#!/bin/bash

#SBATCH -p physical
#SBATCH --time=01:00:00
#SBATCH --job-name=readconvert
#SBATCH --mem=10GB
#SBATCH --out=slurm_%j.out  
#SBATCH --err=slurm_%j.err  
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --mail-type=ALL

mkdir -p necklaceInput/

gunzip -c ${filename} | sed 's/\(^@.\)*\.\(.*\)/\1_\2\/1/g' | gzip > necklaceInput/${filename}
gunzip -c ${filename} | sed 's/\(^+.\)*\.\(.*\)/\1_\2\/1/g' | gzip > necklaceInput/${filename}
EOF
done