#!/bin/perl

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
