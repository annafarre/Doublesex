#!/bin/perl
# ----------------------------------------------------------------------------
# [Author] 	Anna Farre Orteu. 2017
#			afarre@student.unimelb.edu.au
#          	Change names of protein/DRA/RNA records downloaded from NCBI
# ----------------------------------------------------------------------------

# Use:
# perl prots_name_change_4tree.pl filein.fasta fileout.fasta
#
# Ex: 
#	filein: 	>XP_014616309.1 PREDICTED: doublesex- and mab-3-related transcription factor A2 [Polistes canadensis]
#	fileout: 	>Pca_Dsx-Mab_A2_XP_014616309.1
#
#	filein:		>XP_014612935.1 PREDICTED: programmed cell death protein 2 [Polistes canadensis]
#	fileout:	>Pca_programmed_cell_death_protein_2_XP_014612935.1



use warnings;
use Cwd;
use strict;

my $dir = getcwd;
chdir $dir;

my $file = "$dir/$ARGV[0]";
my $outfile = "$dir/$ARGV[1]";

open my $info, $file  or die "Could not open $file: $!";
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";


while( my $line = <$info>){  
	chomp $line;
	if(grep {/mab/i} $line){
		$line =~ m/>(.*?)\s.*?(doublesex[-?]\sand\smab.3.related|isoform\s\S*?|\S*?\sisoform|\S*?\sdoublesex|doublesex.*\sisoform|factor.*truncated-like|\S*?)\s\[(\w)\w*?\s(\w\w)\w*?\].*/i;
		$line = ">$3$4_Dsx-Mab_$2_$1";
		$line =~ s/\s/_/g;
		$line =~ s/,/_/g;
		$line =~ s/__/_/g;
		print $fh "$line\n";
		}
	elsif(grep {/doublesex/i} $line){
		$line =~ m/>(.*?)\s.*?(isoform\s\S*?|\S*?\sisoform|factor.*truncated-like|\S*?)\s\[(\w)\w*?\s(\w\w)\w*?\].*/i;
		$line = ">$3$4_Doublesex_$2_$1";
		$line =~ s/\s/_/g;
		$line =~ s/,/_/g;
		$line =~ s/Doublesex_doublesex/Doublesex/ig;
		$line =~ s/__/_/g;
		print $fh "$line\n";
		}
		
	elsif($line =~ m/>(.*?)\s.*?(LOC\d*\sisoform\s\S*?|LOC\d*|\S*?\sisoform|factor.*truncated-like)\s\[(\w)\w*?\s(\w\w)\w*?\].*/i){
		$line = ">$3$4_$2_$1";
		$line =~ s/\s/_/g;
		$line =~ s/,/_/g;
		$line =~ s/__/_/g;
		$line =~ s/PREDICTED://g;
		print $fh "$line\n";
		}
	elsif($line =~ m/>(.*?)\s.*?(.*)\s\[(\w)\w*?\s(\w\w)\w*?\].*/i){
		$line = ">$3$4_$2_$1";
		$line =~ s/\s/_/g;
		$line =~ s/,/_/g;
		$line =~ s/PREDICTED://g;
		$line =~ s/__/_/g;
		print $fh "$line\n";
	}
	else{ print $fh "$line\n"
	}
}

close $info;
close $file;
close $outfile;
