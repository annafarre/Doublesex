#!/bin/perl
# ----------------------------------------------------------------------------------------------------
# [Author] 	Anna Farre Orteu. 2017
#			afarre@student.unimelb.edu.au
#          	Change header of sample names in count output of ST
# ----------------------------------------------------------------------------------------------------

#Changes header of Necklace (SuperTranscripts) count output
#From:
#	# Program:featureCounts v1.5.3; Command:"/home/afarre/.local/necklace-necklace_v0.9/tools/bin/featureCounts" "-T" "32" "--primary" "-p" "-t" "exon" "-g" "gene_id" "--fraction" "-O" "-f" "-a" "counts/blocks.gtf" "-o" "counts/block.counts" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029066" "/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029067"
#	Geneid	Chr	Start	End	Strand	Length	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029066	/data/projects/punim0352/Wau_necklace/mapped_reads/DRR029067
#To: 
#	# Program:featureCounts v1.5.3; Command:"/home/afarre/.local/necklace-necklace_v0.9/tools/bin/featureCounts" "-T" "32" "--primary" "-p" "-t" "exon" "-g" "gene_id" "--fraction" "-O" "-f" "-a" "counts/blocks.gtf" "-o" "counts/block.counts" "DRR029066" "DRR029067"	
#	Geneid	Chr	Start	End	Strand	Length	DRR029066	DRR029067

#Script changes ONLY the first two lines, all else stays the same


use warnings;
use Cwd;
#use strict;

my $dir = getcwd;
chdir $dir;

my $file = "$dir/$ARGV[0]";
my $outfile = "$dir/$ARGV[1]";

open my $info, $file  or die "Could not open $file: $!";
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";

my $count = 0;
while( my $line = <$info>){  
	chomp $line;
	if($count <2){
		$line =~ s/(\/\S*\/)*(DRR\d*)\.bam/\2/g;
		print $fh "$line\n";
		$count++;
		}
	else{ print $fh "$line\n"
	}
}

close $info;
close $file;
close $outfile;
