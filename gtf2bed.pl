#!/bin/perl

#Convert gtf to bed 


use warnings;
use Cwd;
#use strict;

my $dir = getcwd;
chdir $dir;

my $filein = "$dir/$ARGV[0]";
my $fileout = "$dir/$ARGV[1]";
open my $info, $filein or die "Could not open $filein: $!";
open(my $output, '>', $fileout) or die "Could not open file '$fileout' $!";

my $count=1;

while( my $line = <$info>){  
	chomp $line;
	$line =~ m/(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(.*)/;
	my $chromosome= "$1";
	my $start=$4-1;
	my $end=$5;
	print $output "${chromosome}\t${start}\t${end}\n";
}

close $info;
