#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

#Max Lundberg 2019

#This is script to parse an interpro tsv file and store each information
#of a specific gene on a single row

my ($interpro_output_file,$help);

GetOptions ("interpro_output_file=s"	 =>\$interpro_output_file,
	    "help"	         	 =>\$help	
	    );


my $USAGE = <<"USAGE";

Usage: summarize_interpro_output.pl --interpro_output_file <file.tsv>

--interpro_output_file	output file (.tsv) from interpro
--help			display this message

USAGE

if(!$interpro_output_file){die "Must specify an interpro output file\n$USAGE\n"}


open(INTERPRO_DATA,$interpro_output_file) or die "cannot open $interpro_output_file\n$USAGE\n";


my %interpro_annotations;
my %analyses;


while(my $interpro_input=<INTERPRO_DATA>){

	my ($transcript,$md5,$length,$analysis,$sign_accession,$start,$stop,$score,$status,$date,$ipr_accession,$ipr_annotation)=split("\t",$interpro_input);

	my $annotation=join(":",$analysis,$sign_accession);
	
	

	if(exists($interpro_annotations{$transcript}{$annotation})){
		$interpro_annotations{$transcript}{$annotation}++
		}

	else{
		$interpro_annotations{$transcript}{$annotation}=1
		}


	}

close(INTERPRO_DATA);




print "gene\tannotations\n";


for my $transcript (sort keys %interpro_annotations){

	$transcript=~/([^.]+)/;
	
	print $1,"\t";

	my $index=1;

	for my $annotation (sort keys %{$interpro_annotations{$transcript}}){

		unless($index==1){
			print ","
			}
		print $annotation;

		$index++
		
		}


	print "\n"
		


	}


