#!/bin/bash
# ------------------------------------------------------------------
# [Author] 	Anna Farre Orteu, 2017
#			afarre@student.unimelb.edu.au
#          	Create scripts to run kallisto
# ------------------------------------------------------------------

bold=$(tput bold)
normal=$(tput sgr0)

subject=kallisto
usage="
${bold}DESCRIPTION: 
${normal}Run kallisto for paired end data
use v0.43.0 when multithreading, v0.43.1 crashes

kallisto.sh [options] reads_1.fastq.gz reads_2.fastq.gz ...

where: 
	-o output path
	-i kallisto index file
	-f input path (where fastq files are)

"

# --- Option processing --------------------------------------------
if [ $# == 0 ] ; then
    echo "$usage"
    exit 1;
fi
while getopts ":o:i:f:h" optname
  do
    case "$optname" in
      "i")
		kallisto_index="$(pwd)/$OPTARG"
        ;;
      "o") 
      	mkdir -p $OPTARG
      	
 		if [[ $OPTARG == /* ]] ; then
 			if [ ${OPTARG: -1} == "/" ] ; then
        		outputPath="$OPTARG"
        	else 
        		outputPath="$OPTARG/"
        	fi
       
        else 
        	if [ ${OPTARG: -1} == "/" ] ; then
        		outputPath="$(pwd)/$OPTARG"
        	else 
        		outputPath="$(pwd)/$OPTARG/"
        	fi
        fi
        ;;      
      "f")       	
 		if [[ $OPTARG == /* ]] ; then
 			if [ ${OPTARG: -1} == "/" ] ; then
        		inputPath="$OPTARG"
        	else 
        		inputPath="$OPTARG/"
        	fi
       
        else 
        	if [ ${OPTARG: -1} == "/" ] ; then
        		inputPath="$(pwd)/$OPTARG"
        	else 
        		inputPath="$(pwd)/$OPTARG/"
        	fi
        fi
        ;;      
      "h")
        echo "$usage"
        exit 0;
        ;;
      "?")
        echo "Unknown option $OPTARG"
        exit 0;
        ;;
      ":")
        echo "No argument value for option $OPTARG"
        exit 0;
        ;;
      *)
        echo "Unknown error while processing options"
        exit 0;
        ;;
    esac
  done
  
if [ $OPTIND -eq 1 ]; then 
	printf "\nNo options were passed
	$usage" 
	exit 1;
fi

shift $(($OPTIND - 1))

# -----------------------------------------------------------------
#  SCRIPT LOGIC GOES HERE
# -----------------------------------------------------------------


list="$@"
while read filename 
do

mkdir -p ${outputPath}${subject}_${filename}/
pathScript=${outputPath}${subject}_${filename}/${subject}_${filename}.sh

cat <<-EOF > ${pathScript}
#!/bin/bash

#SBATCH -p physical
#SBATCH --time=03:00:00
#SBATCH --job-name=${subject}_${filename}
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_${subject}_${filename}_%j.out
#SBATCH --err=slurm_${subject}_${filename}_%j.err

/home/afarre/kallisto_linux-v0.43.0/kallisto quant --bootstrap-samples 100 --index=${kallisto_index} --output-dir=${outputPath}${subject}_${filename}/ --threads=8 --single ${inputPath}${filename}_pass_trimmed.fq.gz 

EOF
done<"${list}"

