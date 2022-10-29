#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;


#Max Lundberg 2022

#This is script to find ancestral alleles based on data from outgroup species 

my ($ref_species_vcf,$outgroup1_file,$outgroup2_file,$outgroup1_cov,$outgroup2_cov,$help);

GetOptions ("ref_species_vcf=s"	 =>\$ref_species_vcf,
	    "outgroup1_file=s"	 =>\$outgroup1_file,
            "outgroup2_file=s"	 =>\$outgroup2_file,
	    "outgroup1_cov=s"	 =>\$outgroup1_cov,
	    "outgroup2_cov=s"	 =>\$outgroup2_cov,	
	    "help"	         =>\$help	
	    );

my $USAGE = <<"USAGE";

Usage: find_ancestral_allele_outgroup_species.2.pl --ref_species_vcf <FILE.vcf> --outgroup1_file <FILE> --outgroup2_file <FILE>
--ref_species_vcf	vcf file with data from the reference species
--outgroup1_file	file containing data from outgroup species 1
--outgroup2_file	file containing data from outgroup species 2
--outgroup1_cov		outgroup 1 mean coverage
--outgroup2_cov		outgroup 2 mean coverage
--help			display this message

USAGE

if(!$ref_species_vcf){die "\nMust specify a reference vcf file!\n$USAGE\n"}
if(!$outgroup1_file){die "\nMust specify an outgroup species 1 file!\n$USAGE\n"}
if(!$outgroup2_file){die "\nMust specify an outgroup species 2 file!\n$USAGE\n"}
if(!$outgroup1_cov){die "\nMust specify an outgroup mean cov!\n$USAGE\n"}
if(!$outgroup2_cov){die "\nMust specify an outgroup mean cov!\n$USAGE\n"}


my @nuc=("A","C","G","T");

#open the outgroup allele counts
open(OUTGROUP1,$outgroup1_file) or die "cannot open $outgroup1_file !\n";

my %allele_counts_out1;

while(my $out1_input=<OUTGROUP1>){

	chomp($out1_input);

	my @out1_input=split("\t",$out1_input); 
	
	my $interval=join(":",@out1_input[0..1]);

	for(0..3){
		if($out1_input[$_+2]>0){
			$allele_counts_out1{$interval}{$nuc[$_]}=$out1_input[$_+2]
			}

		}
	$allele_counts_out1{$interval}{"all"}=$out1_input[2]+$out1_input[3]+$out1_input[4]+$out1_input[5];
	}

close(OUTGROUP1);

my %allele_counts_out2;

open(OUTGROUP2,$outgroup2_file) or die "cannot open $outgroup2_file !\n";

while(my $out2_input=<OUTGROUP2>){

	chomp($out2_input);

	my @out2_input=split("\t",$out2_input); 

	my $interval=join(":",@out2_input[0..1]);
	
	for(0..3){
		if($out2_input[$_+2]>0){
			$allele_counts_out2{$interval}{$nuc[$_]}=$out2_input[$_+2]
			}

		}


	$allele_counts_out2{$interval}{"all"}=$out2_input[2]+$out2_input[3]+$out2_input[4]+$out2_input[5];
	}

close(OUTGROUP2);


my ($ref_count,$alt_count);


if($ref_species_vcf=~/gz$/){
	open(VCF_REF,"zcat $ref_species_vcf |") or die "cannot open vcf ref file $ref_species_vcf !\n";	
	}
else{
	open(VCF_REF,$ref_species_vcf) or die "cannot open vcf ref file $ref_species_vcf ! \n";
	}


while(my $input=<VCF_REF>){

	next if($input=~/^#/);

	my @input=split("\t",$input); 

	my ($scaffold,$pos)=@input[0..1];

	my $interval=join(":",$scaffold,$pos);

	my ($ref,$alt)=@input[3..4];

#	print $ref,"\t",$alt,"\n";

	my ($out1_ref,$out1_alt,$out2_ref,$out2_alt);

	my $low_cov1="no";
	my $low_cov2="no";

	if(exists($allele_counts_out1{$interval}{$ref})){
		$out1_ref=$allele_counts_out1{$interval}{$ref};
		}
	else{
		$out1_ref=0
		}
	if(exists($allele_counts_out1{$interval}{$alt})){
		$out1_alt=$allele_counts_out1{$interval}{$alt};	
		}
	else{
		$out1_alt=0;
		}
	if(exists($allele_counts_out1{$interval}{"all"})){
		if($allele_counts_out1{$interval}{"all"}<$outgroup1_cov/3){$low_cov1="yes"}
		}
	
	if(exists($allele_counts_out2{$interval}{$ref})){
		$out2_ref=$allele_counts_out2{$interval}{$ref};
		}
	else{
		$out2_ref=0
		}
	if(exists($allele_counts_out2{$interval}{$alt})){
		$out2_alt=$allele_counts_out2{$interval}{$alt};	
		}
	else{
		$out2_alt=0;
		}

	if(exists($allele_counts_out2{$interval}{"all"})){
		if($allele_counts_out2{$interval}{"all"}<$outgroup2_cov/3){$low_cov2="yes"}
		}
	

	#Evaluate
	
	my $out1_eval;
 
	if($low_cov1 eq "yes"){
		$out1_eval="low cov"
		}
	else{
		if($out1_ref > 0 && $out1_alt==0){$out1_eval="hom_ref"}
		elsif($out1_ref==0 && $out1_alt>0){$out1_eval="hom_alt"}
		elsif($out1_ref>0 && $out1_alt>0){$out1_eval="both"}
		elsif($out1_ref==0 && $out1_alt==0){$out1_eval="no_cov"}
		}

	my $out2_eval;

	if($low_cov2 eq "yes"){
		$out2_eval="low cov"
		}
	else{
		if($out2_ref > 0 && $out2_alt==0){$out2_eval="hom_ref"}
		elsif($out2_ref==0 && $out2_alt>0){$out2_eval="hom_alt"}
		elsif($out2_ref>0 && $out2_alt>0){$out2_eval="both"}
		elsif($out2_ref==0 && $out2_alt==0){$out2_eval="no_cov"}
		}

	my $consensus=".";

	if($out1_eval=~/hom/ && $out1_eval eq $out2_eval){
		if($out1_eval eq "hom_ref"){
			$consensus="REF"	
			}		

		else{
			$consensus="ALT"
			}
		}

	print $scaffold,"\t",$pos,"\t",$ref,"\t",$alt,"\t",$out1_ref,"\t",$out1_alt,"\t",$out2_ref,"\t",$out2_alt,"\t",$out1_eval,"\t",$out2_eval,"\t",$consensus,"\n";

}

