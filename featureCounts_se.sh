#!/bin/bash
# ------------------------------------------------------------------
# [Author] 	Anna Farre Orteu, 2017
#			afarre@student.unimelb.edu.au
#          	Create scripts to run featureCounts
# ------------------------------------------------------------------

bold=$(tput bold)
normal=$(tput sgr0)

subject=featureCounts
usage="
${bold}DESCRIPTION: 
${normal}Create scripts that will run featureCounts
Count number of reads per exon and transcript in RNAseq data.
	
featureCounts.sh [options] list_SRA.txt

where: 
	-a gtf file with annotation of chromosome and trasncript coordinates (with path)
	-g reference genome (with path)
	-o name of output (with path)

"

# --- Option processing --------------------------------------------
if [ $# == 0 ] ; then
    echo "$usage"
    exit 1;
fi
while getopts ":a:o:g:h" optname
  do
    case "$optname" in
      "a")
 		gtf=$OPTARG
        ;;
      "g") 
      	genome=$OPTARG
        ;;  
      "o") 
      	mkdir -p $OPTARG
      	output=$OPTARG
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

pathScript=${outputPath}featureCounts.sh

cat <<-EOF > ${pathScript}
#!/bin/bash

#SBATCH -p physical
#SBATCH --time=02:00:00
#SBATCH --job-name=featureCounts
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_featureCounts_%j.out
#SBATCH --err=slurm_featureCounts_%j.err

featureCounts -J -t exon -g transcript_id -f -O -T 3 -s 0 -G ${genome} \\
-a ${gtf} \\
-o ${output} ${list} \\

EOF
done<"${list}"
