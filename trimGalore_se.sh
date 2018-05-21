#!/bin/bash
# ------------------------------------------------------------------
# [Author] 	Anna Farre Orteu, 2017
#			afarre@student.unimelb.edu.au
#          	Create scripts to run TrimGalore!
# ------------------------------------------------------------------

bold=$(tput bold)
normal=$(tput sgr0)

subject=TrimGalore!
usage="
${bold}DESCRIPTION: 
${normal}Create scripts that will run TrimGalore!
Trims adapters, low quality read ends and runs FastQC
Single end data ONLY
gz compressed data
extention must be: 
	*_pass.fastq.gz 
	
trimGalore.sh [options] list_SRA.txt

where: 
	-o output path
	-i input path

${bold}REQUIEREMENTS
${normal}cutadapt
fastqc
"

# --- Option processing --------------------------------------------
if [ $# == 0 ] ; then
    echo "$usage"
    exit 1;
fi
while getopts ":o:i:h" optname
  do
    case "$optname" in
      "i")
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


list=$@

while read filename 
do

pathScript=${outputPath}trimGalore_${filename}.sh

cat <<-EOF > ${pathScript}
#!/bin/bash

#run TrimGalore on raw data 

#SBATCH -p physical
#SBATCH --time=03:00:00
#SBATCH --job-name=trimGalore_${filename}
#SBATCH --mem=25GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_trimGalore_${filename}_%j.out
#SBATCH --err=slurm_trimGalore_${filename}_%j.err


#path to list of files to use

module load Python
module load fastqc

/home/afarre/.local/TrimGalore-0.4.5/trim_galore --fastqc --gzip --output_dir ${outputPath} ${inputPath}${filename}_pass.fastq.gz
EOF
done <"${list}"