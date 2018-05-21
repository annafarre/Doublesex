#!/bin/bash
# ------------------------------------------------------------------
# [Author] 	Anna Farre Orteu
#			afarre@student.unimelb.edu.au
#          	Fastq-dump for paired and single end sra files
# ------------------------------------------------------------------

bold=$(tput bold)
normal=$(tput sgr0)

subject=fatsq-dump
version=0.1.0
usage="
${bold}DESCRIPTION: 
${normal}Create scripts that will fastq-dump sra files for paired and single end libraries
Paired end and single end have to be run ${bold}separately
    
${bold}USAGE: 
${normal}dumpthemall [OPTION] [list of SRA files to fastq-dump] 

where:
    -h  show this help text
    -p	path to SRA files
    -r	library type: ${bold}p ${normal}(paired end) and ${bold}s ${normal}(single end)

"
# --- Option processing --------------------------------------------
if [ $# == 0 ] ; then
    echo "$usage"
    exit 1;
fi

while getopts ":p:r:h" optname
  do
    case "$optname" in
      "p")
        if [ ${OPTARG: -1} == "/" ] ; then
        	path="$OPTARG"
        else 
        	path="$OPTARG/"
        fi
        ;;
      "r")
        if [ $OPTARG == "p" ] ; then
        	options="--defline-seq '@$sn[_$rn]/$ri' --split-3 --readids --skip-technical --clip --read-filter pass --dumpbase --gzip"       
       		echo $options
        elif [ $OPTARG == "s" ] ; then
        	options="--defline-seq '@$sn[_$rn]/$ri' --readids --skip-technical --clip --read-filter pass --dumpbase  --gzip"          
        	echo $options
        else 
        	echo "Invalid argument. Specify p (paired) or s (single) end library"        
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

param1=$1

# -----------------------------------------------------------------
#  SCRIPT LOGIC GOES HERE
# ------------------------------------------------------------
while read line
do 
	pathScript=${path}${line}_fastq-dump.sh
	cat <<-EOF > ${pathScript}
	#!/bin/bash 
 	#SBATCH -p physical 
 	#SBATCH --time=03:00:00 
	#SBATCH --nodes=1 
	#SBATCH --ntasks=1 
	#SBATCH --cpus-per-task=1 
	#SBATCH --job-name=fastq-dump 
	#SBATCH --mem-per-cpu=50000 
	#SBATCH --mail-type=ALL 
	#SBATCH --mail-user=afarre@student.unimelb.edu.au 
	#SBATCH --out=slurm_%j.out 
	#SBATCH --err=slurm_%j.err 

	fastq-dump $options ${path}${line}.sra 
	EOF
	
done <"$1"

