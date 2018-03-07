#!/bin/bash
# ------------------------------------------------------------------
# [Author] 	Anna Farre Orteu
#			afarre@student.unimelb.edu.au
#          	Fastqc for paired and single end fastq files
# ------------------------------------------------------------------

bold=$(tput bold)
normal=$(tput sgr0)

subject=fatsq-dump
version=0.1.0
usage="
${bold}DESCRIPTION: 
${normal}Run fastqc in paired and single end fastq files
    
${bold}USAGE: 
${normal}fastqc.sh [OPTION] [list of fastq files] 

where:
    -h  show this help text
    -o	path where output has to be saved
    
"
# --- Option processing --------------------------------------------
if [ $# == 0 ] ; then
    echo "$usage"
    exit 1;
fi

while getopts ":o:h" optname
  do
    case "$optname" in
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

param1=$@

# -----------------------------------------------------------------
#  SCRIPT LOGIC GOES HERE
# ------------------------------------------------------------

for filename in `ls ${param1}`
do

pathScript=fastqc_${filename}.sh
cat <<-EOF > ${pathScript}
#!/bin/bash

#SBATCH -p cloud
#SBATCH --time=09:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=fastqc_${filename}
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_%j.out
#SBATCH --err=slurm_%j.err 

module load fastqc

fastqc -o ${outputPath} -t 8 ${filename}
EOF
done

