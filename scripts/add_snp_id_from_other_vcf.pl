#!/usr/bin/perl

use warnings;
use strict; 
use Getopt::Long;


#Max Lundberg 2020

#This scripts transfers a SNPid from another vcf file based on overlapping positions



my ($ref_vcf_file,$help);

GetOptions ("ref_vcf_file=s"	 =>\$ref_vcf_file,
	    "help"	         =>\$help	
	    );

my $USAGE = <<"USAGE";

Usage: add_snp_id_from_other_vcf.pl --ref_vcf_file <file> --target_vcf_file <file>

--ref_vcf_file		reference vcf file from which SNPids should be transferred
--help			display this message

USAGE

if($help){
	die $USAGE;
	}

if(!$ref_vcf_file){
	die "\nNeed to specify ref vcf file\n$USAGE"
	}


#my $ref_vcf_file="snparray_data_filt.new.vcf";

open(REF_VCF,$ref_vcf_file) or die "cannot open ref vcf file\n";


my %snp_ids;

while(my $ref_input=<REF_VCF>){

	next if($ref_input=~/^#/);

	my @ref_input=split("\t",$ref_input);
	
	my $pos=join("-",$ref_input[0],$ref_input[1]);

	my $ID=$ref_input[2];
		
	$snp_ids{$pos}=$ID;

	}

close(REF_VCF); 


$,="\t";

#open(TARGET_VCF,$target_vcf_file) or die "cannot open target vcf file\n";

while(my $input=<STDIN>){

	if($input=~/^#/){print $input;next};

	my @input=split("\t",$input);
	
	my $pos=join("-",@input[0..1]);

	print @input[0..1],$snp_ids{$pos},@input[3..$#input];
	

	}

#close(TARGET_VCF); 



