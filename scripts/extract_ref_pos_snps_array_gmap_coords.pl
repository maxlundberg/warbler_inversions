#!/usr/bin/perl

use warnings;
use strict; 

#Max Lundberg 2019

#This script is used to extract position of array SNPs in the ww genome
#with gmap and with a "coords" output. It assumes that the focal SNP
#is marked with a low letter ambig nucleotide (r,y,k or m).It will provide
#a list with the position of the SNP and the reference base of this position


my $probe;
my %probe_data;
 

while(my $input=<STDIN>){

	chomp($input);

	if($input=~m/^>(.+)/){
		$probe=$1;
		}
	else{
		if($input=~m/(r|y|k|m)\s+([+-])([^:]+):(\d+)\s\d+\s([A-Z])/){

			my $allele=$1;
			my $orientation=$2;		
			my $scaffold=$3;
			my $pos=$4;
			my $ref=$5;

	
		if($orientation eq "-"){
			$ref=~tr/ACGT/TGCA/;		
			}



		$probe_data{$probe}=join("\t",$scaffold,$pos,$probe,$ref); 
	
		}


	}
	}

print "##Scaffold\tpos\tSNP_id\tref\n"; 

for(keys %probe_data){
	print $probe_data{$_},"\n"; 
	}
