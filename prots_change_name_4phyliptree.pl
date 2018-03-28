#!/bin/perl
# ----------------------------------------------------------------------------------------------------
# [Author] 	Anna Farre Orteu. 2017
#			afarre@student.unimelb.edu.au
#          	Change names of protein/DRA/RNA records downloaded from NCBI to phylip compatible names
# ----------------------------------------------------------------------------------------------------

# Use:
# perl prots_change_name_4phyliptree.pl filein.fasta fileout.fasta
#
# Ex: 
#	filein: 	>XP_014616309.1 PREDICTED: doublesex- and mab-3-related transcription factor A2 [Polistes canadensis]
#	fileout: 	>014616309
#
#	filein:		>XP_014612935.1 PREDICTED: programmed cell death protein 2 [Polistes canadensis]
#	fileout:	>014612935


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
	if($line =~ m/>.*?([0-9a-fA-F]*)\..\s.*/){
	$line = ">$1";
	$line =~ s/\s/_/g;
	$line =~ s/,/_/g;
	print $fh "$line\n";
	}
	else{ print $fh "$line\n"
	}
}

close $info;
close $file;
close $outfile;
