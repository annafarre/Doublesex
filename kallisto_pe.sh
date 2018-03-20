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

Only for ${bold}paired-end ${normal}samples

kallisto.sh [options] list.txt

where: 
	-o output path
	-i kallisto index file
	-f input path (where fastq files are)
	-r extension for R1 (read 1) e.g _pass_1_val_1.fq.gz
	-l extension for R2 (read 2) e.g _pass_1_val_1.fq.gz
	
	list.txt contains all SRA sample IDs

"

# --- Option processing --------------------------------------------
if [ $# == 0 ] ; then
    echo "$usage"
    exit 1;
fi
while getopts ":o:i:f:r:l:h" optname
  do
    case "$optname" in
      "i")
		if [[ $OPTARG == /* ]] ; then
        	kallisto_index="$OPTARG"
        else 
       		kallisto_index="$(pwd)/$OPTARG"
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
      "r")
      	read1Name="$OPTARG"
      	;;
      "l")
        read2Name="$OPTARG"
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


/home/afarre/.local/kallisto_linux-v0.43.0/kallisto quant --bootstrap-samples 100 --index=${kallisto_index} --output-dir=${outputPath}${subject}_${filename}/ --threads=8 ${inputPath}${filename}${read1Name} ${inputPath}${filename}${read2Name}

EOF
done<"${list}"

