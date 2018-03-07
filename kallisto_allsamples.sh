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


array=( "$@" )
array=( "${array[@]/#/$(pwd)/}" )
printf -v list '%s ' "${array[@]}"
printf '%s ' "${array[@]}"


pathScript=${outputPath}${subject}.sh

echo "#!/bin/bash" >${pathScript}


echo "#SBATCH -p physical" >>${pathScript}
echo "#SBATCH --time=09:00:00" >>${pathScript}
echo "#SBATCH --job-name=${subject}" >>${pathScript}
echo "#SBATCH --cpus-per-task=8" >>${pathScript}
echo "#SBATCH --mem=40GB" >>${pathScript}
echo "#SBATCH --mail-type=ALL" >>${pathScript}
echo "#SBATCH --mail-user=afarre@student.unimelb.edu.au" >>${pathScript}
echo "#SBATCH --out=slurm_${subject}_%j.out" >>${pathScript}
echo "#SBATCH --err=slurm_${subject}_%j.err" >>${pathScript}


echo "/home/afarre/kallisto_linux-v0.43.0/kallisto quant --index=${kallisto_index} --output-dir=${outputPath} --plaintext --threads=8 ${list}" >>${pathScript}
