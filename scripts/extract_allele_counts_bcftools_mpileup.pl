#!/usr/bin/perl

use warnings; 
use strict; 

while(my $input=<STDIN>){

my %alleles;

#Ignore comment lines and indels
next if($input=~/^#/ || $input=~/INDEL/);

my @input=split("\t",$input);

my ($scaffold,$pos)=@input[0..1];

my $ref=$input[3];

$input[7]=~/AD=([^;]+)/;

my $depths=$1;

my @depths=split(",",$depths);

$alleles{$ref}=$depths[0];

if($input[4] ne "<*>"){
	my @alts=split(",",$input[4]);
	for(0..$#alts){
		$alleles{$alts[$_]}=$depths[$_+1];
		}
	}

print $scaffold,"\t",$pos;

for("A","C","G","T"){
	if(exists($alleles{$_})){	
		print "\t$alleles{$_}";
		}
	else{print "\t0"}
	}

print "\n";

}
