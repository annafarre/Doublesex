#!/bin/bash
# ------------------------------------------------------------------
# [Author] 	Anna Farre Orteu, 2017
#			afarre@student.unimelb.edu.au
#          	Create scripts to run STAR
# ------------------------------------------------------------------

bold=$(tput bold)
normal=$(tput sgr0)

subject=STAR
version=0.1.0
usage="
${bold}DESCRIPTION: 
${normal}Create scripts that will run STAR. 
Align paired and single end libraries to a genome.
${bold}Paired ${normal}end and ${bold}single ${normal}end have to be run ${bold}separately
    
${bold}USAGE: 
${normal}starthemall [OPTION] [list of junction files (*optional)] [list of fastq files to align] 

where:
    -h  show this help text
    -i	path to input, the fastq files
    -g	specifes path to the genome directory where genome indices where generated
    -o	path where output should be saved (STAR scripts and output)
    -r	library type: ${bold}p ${normal}(paired end) and ${bold}s ${normal}(single end)
    -c	Are files compressed (*.gz)? Yes/No 
    -p	Partition to be use in Spartan: physical/cloud
    -d	dependencies. List of job ids (e.g. 123:456:789) Job can begin after the specified jobs have run to completion with an exit code of zero
    
    *list of junction files can be passed as wildcard expressions
    *list of fastq is a txt file - list of SRA accessions -> *SRR_Acc_List.txt
"
# --- Option processing --------------------------------------------
if [ $# == 0 ] ; then
    echo "$usage"
    exit 1;
fi

while getopts ":i:g:r:o:d:c:p:h" optname
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
      "g")
 		if [[ $OPTARG == /* ]] ; then
        	sampleGenome="$OPTARG"
        else 
        	sampleGenome="$(pwd)/$OPTARG"
        fi
        ;;
      "r")
        if [ $OPTARG == "p" ] ; then
        	options="--runThreadN 8 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No --outFileNamePrefix \${output}STAR_\${sample}"       
        	samples="--genomeDir \${genome} --readFilesIn \${samplepath}_pass_1_val_1\${file_extension} \${samplepath}_pass_2_val_2\${file_extension}"
        elif [ $OPTARG == "s" ] ; then
        	options="--runThreadN 8 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No --outFileNamePrefix \${output}STAR_\${sample}"
        	samples="--genomeDir \${genome} --readFilesIn \${samplepath}_pass_trimmed\${file_extension}"          
        else 
        	echo "Invalid argument. Specify p (paired) or s (single) end library"        
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
      "b")
 		if [[ $OPTARG == /* ]] ; then
        	sampleGff="$OPTARG"
        else 
        	sampleGff="$(pwd)/$OPTARG"
        fi
        ;;       
      "c")
      	if [ $OPTARG == "Yes" ] ; then
       		compressed="--readFilesCommand zcat"
        	file_extension=".fq.gz"
      	elif [ $OPTARG == "No" ] ; then
      		file_extension=".fq"
      	else  
        	echo "Invalid argument. Specify file format: Yes (compressed) or No (not compressed)" 
      	fi
        ;;  
      "p")
        partition="$OPTARG"
        ;;
      "d")
        dependencies="#SBATCH --dependency=afterok:$OPTARG"
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
 
if [[ $compressed ]]; then 
 	options="${options} ${compressed}"
fi

if [ $OPTIND -eq 1 ]; then 
	printf "\nNo options were passed
	$usage" 
	exit 1;
fi

shift $(($OPTIND - 1))

list=$@

junctions_array=`printf "$(pwd)/%s\n" ${list} | grep ".*SJ.out.tab"`
printf -v junction_files '%s ' ${junctions_array}

if [[ ${junction_files} ]]; then 
 	options="${options} --sjdbFileChrStartEnd \${junctions}"
fi

list=`printf '%s\n' $list | grep -v ".*SJ.out.tab"`



# -----------------------------------------------------------------
#  SCRIPT LOGIC GOES HERE
# -----------------------------------------------------------------


while read filename 
do

pathScript=${outputPath}STAR_${filename}.sh

cat <<-EOF > ${pathScript}
#!/bin/bash
 
#SBATCH -p ${partition}
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
${dependencies}
#SBATCH --job-name=STAR_${filename}
##SBATCH --mem=25GB
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=afarre@student.unimelb.edu.au
#SBATCH --out=slurm_${filename}_%j.out
#SBATCH --err=slurm_${filename}_%j.err
 
samplepath=${inputPath}${filename}
sample=${filename}
genome=${sampleGenome}
output=${outputPath}
junctions=\"${junction_files}\"
file_extension=${file_extension}
 
module load STAR
 
STAR ${options} ${samples}
EOF 

done <"${list}"

