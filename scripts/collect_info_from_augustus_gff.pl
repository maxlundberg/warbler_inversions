#!/usr/bin/perl

use warnings;
use strict;

#Max Lundberg 2019

#This script collects information about genes and transcripts from an augustus gff3 file

#Usage: cat augustus.predictions.gff3 | collect_info_from_augustus_gff.pl > summary_augustus.tbl



my ($gene,$transcript);

my %gene_info;
my %transcript_support;
my %gene_cds_number;
my %gene_cds_length;

my $cds_info_flag=0;

my %transcript_cds_protein_support;
my %transcript_cds_rnaseq_support;


while(my $input=<STDIN>){


	chomp($input); 

	if($input=~/start gene (\w+)/){
	
		$gene=$1;
		

		}
	else{
		if($input!~/^#/){
			
			my ($scaffold,$source,$feature_type,$start,$end,$score,$strand,$phase,$attributes)=split("\t",$input);
			
			if($feature_type eq "gene"){
				$gene_info{$gene}.="$scaffold\t$start\t$end\t$strand";
				}
			if($feature_type eq "transcript"){
				
				$attributes=~/ID=([^;]+)/;
				$transcript=$1;
				}
			if($feature_type eq "CDS"){
				my $cds_id=join("-",$scaffold,$start,$end);
				if(exists($gene_cds_number{$gene}{$cds_id})){
					$gene_cds_number{$gene}{$cds_id}++;
					}
				else{
					$gene_cds_number{$gene}{$cds_id}=1;
					$gene_cds_length{$gene}+=$end-$start+1;				
					}
				}


			}
		
		else{
			if($input=~/any source\)\: (\d+)/){
				my $support=$1;
				$transcript_support{$gene}{$transcript}=$support;
				}
			elsif($input=~/CDS exons: (\d+)\/(\d+)/){
				my $supported_cds=$1;	
				my $total_cds=$2;
				$cds_info_flag=1
				}
			else{
				if($input!~/CDS introns/ && $cds_info_flag==1){
					
					$input=~/#\s+([A-Z])\:\s+(\d+)/;
					my ($type,$support)=($1,$2);
						
					if($type eq "P"){
						$transcript_cds_protein_support{$gene}{$transcript}=$support
						}
					else{
						$transcript_cds_rnaseq_support{$gene}{$transcript}=$support
						}	
					
					

					}	
				else{
					if($input=~/CDS introns/){
						$cds_info_flag=0
						}
  					}						
				}
					#print $1,"\t",$2,"\n"
			}

					
					
		}
	}


print "gene\tscaffold\tstart\tend\tstrand\tn_transcripts\thighest_support\ttranscripts\tcds_number\tcds_length\tcds_supported_by_rnaseq\tcds_supported_by_protein\n";

for my $gene (keys %transcript_support){

	print "$gene\t$gene_info{$gene}\t";


	my $n_transcripts=scalar(keys %{$transcript_support{$gene}});
		
 	my $highest_support=0;
	my $transcripts;


	my $highest_rnaseq_cds_support=0;
	my $highest_protein_cds_support=0;		
	
	for my $transcript (keys %{$transcript_support{$gene}}){

		my $support=$transcript_support{$gene}{$transcript};

		if($support>$highest_support){
			$highest_support=$support
			}
		$transcripts.="$transcript-($support);";

		if(exists($transcript_cds_rnaseq_support{$gene}{$transcript})){
			my $rnaseq_cds_support=$transcript_cds_rnaseq_support{$gene}{$transcript};	

			if($rnaseq_cds_support>$highest_rnaseq_cds_support){	
				$highest_rnaseq_cds_support=$rnaseq_cds_support
				}		
			}

		if(exists($transcript_cds_protein_support{$gene}{$transcript})){
			my $protein_cds_support=$transcript_cds_protein_support{$gene}{$transcript};	

			if($protein_cds_support>$highest_protein_cds_support){	
				$highest_protein_cds_support=$protein_cds_support
				}		
			}	
	

	

		}
		
	my $n_cds=scalar(keys %{$gene_cds_number{$gene}});

	print "$n_transcripts\t$highest_support\t$transcripts\t$n_cds\t$gene_cds_length{$gene}\t$highest_rnaseq_cds_support\t$highest_protein_cds_support\n";		
	
	}
