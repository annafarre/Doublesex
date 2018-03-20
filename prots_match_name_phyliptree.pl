#!/bin/perl

#use warnings;
use Cwd;
#use strict;
use Data::Dumper;

my $dir = getcwd;
chdir $dir;

my $fasta = "$dir/$ARGV[0]";
my $phylip = "$dir/$ARGV[1]";
my $outfile = "$dir/$ARGV[2]";

open my $longname, $fasta  or die "Could not open $fasta: $!";
open my $shortname, $phylip  or die "Could not open $phylip: $!";
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";

#my $count = 0;
while( my $line = <$longname>){  
	chomp $line;
	my $longline = $line;
	$longline =~ s/>//;
	if($line =~ m/>.*?([0-9a-fA-F]*)\..*/){
	$line = "$1";
	$line =~ s/\s/_/g;
	$line =~ s/,/_/g;
	$short2long{$1}=$longline;
	#$count++;
	#print "$longline\n";
	}
}
#print "$count\n";


while( my $line = <$shortname>){  
	chomp $line;
	@matches = ( $line =~ m/[\(|,]([\w|\d]+?):/g);
	#print $line;
	#print @matches;
	#my $count = 0;
	for $name (@matches){
		if (exists $short2long{$name}){
			$line =~ s/([\(|,])$name(:)/$1$short2long{$name}$2/g;
			#print "$name\n";
			#print "$short2long{$name}\n";
			#$count++;
		}
	}
	print $fh $line;
	#print "$count\n";
}


close $longname;
close $shortname;
close $outfile;
