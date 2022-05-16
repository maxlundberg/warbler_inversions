#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;


#Max Lundberg 2021

#Usage: cat <vcf file> | change_name_sequences_bionano_sv_vcf.pl --key_file <key file> 

my $key_file;

GetOptions ("key_file=s"	 =>\$key_file);

if(!$key_file){
	die "\nNeed to specify a key file\n"
	}


#################

#Process the key file

#################


open(KEY_FILE,$key_file) or die "cannot open key file\n";


my %key_sequences;

while(my $key_input=<KEY_FILE>){

	if($key_input=~/^\d+/){
	
		my @key_input=split("\t",$key_input);
		
		my $new_id=join("","chr",$key_input[0]);
		my $old_id=$key_input[1];
		
		$key_sequences{$new_id}=$old_id
		
		}

	}

close(KEY_FILE);


#################

#Process the vcf file

#################

while(my $vcf_input=<STDIN>){
	
	if($vcf_input=~/(chr[0-9]+)/){
		
		my $new_id=$1;
	
		$vcf_input=~s/$new_id/$key_sequences{$new_id}/g;
	
	}

	print $vcf_input
	
	}










