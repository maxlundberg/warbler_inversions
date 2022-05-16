#!/usr/bin/perl

use warnings;
use strict;

#This script will identify overlapping contigs within
#bionano super-scaffolds


#Usage:  cat <file.agp> | get_supercontigs_from_bionano_scaffolding.pl

my $contig_index=1;
my $super_scaffold; 

my %super_contigs;
my %super_contig_count;



while(my $input=<STDIN>){

	my ($Obj_Name,$Obj_Start,$Obj_End,$PartNum,$Compnt_Type,$CompntId_GapLength,$CompntStart_GapType,$CompntEnd_Linkage,$Orientation_LinkageEvidence)=split("\t",$input);

	if($Obj_Name=~/^Super/){


		#For the first entry
		if(!$super_scaffold){
			$super_scaffold=$Obj_Name			
			}
		else{
			if($Obj_Name ne $super_scaffold){
				$super_scaffold=$Obj_Name;
				$contig_index++
				}
			}
			

		if($Compnt_Type eq "N" && $CompntId_GapLength!=13){
			$contig_index++
			}




		#If the line contains scaffold information 
		if($Compnt_Type eq "W"){

			$super_contigs{$contig_index}.="$contig_index\t$Obj_Name\t$Obj_Start\t$Obj_End\t$CompntId_GapLength\n";
			$super_contig_count{$contig_index}++;
			
			
	
			}
	
		}

	}

for (sort {$a <=> $b} keys %super_contigs){

	if($super_contig_count{$_}>1){

		print $super_contigs{$_}
		}
	}




