#!/bin/bash

bold=$(tput bold)
normal=$(tput sgr0)

usage="
${bold}NAME: 
${normal}runthemall
$0

${bold}DESCRIPTION: 
${normal}A script to run them all
Program to run a list of bash scripts in a Slurm cluster as independent jobs
Wildcards can be used
    
${bold}USAGE: 
${normal}runthemall [-h] [list of scripts to run] 

where:
    -h  show this help text
    
${bold}EXAMPLES:
    ${normal}runthemall *.sh
    runthemall script1.sh script2.sh ...
"

if [ "$1" == "-h" ]; then
  echo "$usage"
  exit 0

elif [ $# -lt 1 ]; then
echo -e 1>&2 "$0: missing list of scripts to run

${bold}Usage: 
\t${normal}runthemall [-h] [list of scripts to run] 

use option -h for more information
	"
	exit 2

else 
	for script in `ls $@`
		do
			sbatch $script
		done
fi
