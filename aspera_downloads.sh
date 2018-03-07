#!/bin/bash
# Download a bunch of *.sra files from the NCBI SRA, using the aspera client


max_bandwidth_mbps=5000

# SRA files written line by line in a STDIN file

while read file
do 
 /Users/afarre/Applications/Aspera\ Connect.app/Contents/Resources/ascp \
 -i /Users/afarre/Applications/Aspera\ Connect.app/Contents/Resources/asperaweb_id_dsa.openssh \
 -k1 -QTr -l${max_bandwidth_mbps}m \
 anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/${file:0:3}/${file:0:6}/${file}/${file}.sra /Users/afarre/SRA/
done <"$1"