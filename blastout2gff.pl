#!/bin/perl

#convert blast output to gff
#If strand exists in blastOut use commented lines

#code adds count # to exon id 
#to differentiate exons in those cases that gene names are the same (such as when using SuperTranscripts)


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
	$line =~ m/(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)/;
#	$line =~ m/(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)/;
	my $seqid=$2;
	my $source="CustomAnnotAOrteu";
	my $type="exon";
	my $start=$9;
	my $end=$10;
	my $score=".";
#	my $strand=$13;
	my $strand=".";
	my $phase="0";
	my $attributes= "exon_id=$1.$count";
#	${strand} =~ s/minus/-/g;
#	${strand} =~ s/plus/+/g;
	${seqid} =~ s/ref\|(.*)\|/$1/g;
#	print $output "${seqid}\t${source}\t${type}\t${start}\t${end}\t${score}\t${strand}\t${phase}\t${attributes}\n";
	print $output "${seqid}\t${source}\t${type}\t${start}\t${end}\t${score}\t${strand}\t${phase}\t${attributes}\n";
	$count++;
}

close $info;
