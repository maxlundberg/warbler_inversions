#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

#Max Lundberg 2022

#This is script to polarize alleles in a vcf file based on outgroup
#species information. 

my ($ref_species_vcf,$outgroup_info,$help);

GetOptions ("ref_species_vcf=s"	 =>\$ref_species_vcf,
	    "outgroup_info=s"	 =>\$outgroup_info,
	    "help"	         =>\$help	
	    );


my $USAGE = <<"USAGE";

Usage: add_outgroup_geno_vcf.pl --ref_species_vcf <FILE.vcf> --outgroup_info <FILE> 
--ref_species_vcf	vcf file with data from the reference species
--outgroup_info		file containing data from outgroup
--help			display this message

USAGE

if(!$ref_species_vcf){die "\nMust specify a reference vcf file!\n$USAGE\n"}
if(!$outgroup_info){die "\nMust specify an outgroup info file!\n$USAGE\n"}



#Open outgroup info file 

my %ancestral_allele;

open(OUTINFO,$outgroup_info) or die "cannot open outgroup info file $outgroup_info\n";

while(my $outgroup_info=<OUTINFO>){

	chomp($outgroup_info);
	
	my @outgroup_info=split("\t",$outgroup_info);
	my $out_interval=join(":",@outgroup_info[0..1]);
	
	my $ancestral=$outgroup_info[10];

	$ancestral_allele{$out_interval}=$ancestral
	}



$,="\t";


if($ref_species_vcf=~/gz$/){
	open(VCF_REF,"zcat $ref_species_vcf |") or die "cannot open vcf ref file $ref_species_vcf !\n";	
	}
else{
	open(VCF_REF,$ref_species_vcf) or die "cannot open vcf ref file $ref_species_vcf ! \n";
	}




while(my $input=<VCF_REF>){

	if($input=~/^##/){
		print $input
		}
	elsif($input=~/#CHROM/){
		chomp($input);
		print "$input\tOUTGROUP\n"
		}
			
	else{

		chomp($input);	
			
		my @input=split("\t",$input); 

		my ($scaffold,$pos)=@input[0..1];

		my $interval=join(":",$scaffold,$pos);

		my ($ref,$alt)=@input[3..4];

		if(exists($ancestral_allele{$interval}) && $ancestral_allele{$interval} ne "."){

			print $input;

			if($ancestral_allele{$interval} eq "REF"){
				print "\t0|0\n"
				}
			else{
				print "\t1|1\n"			
				}
				
			}
	
		}



	}






