#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;


#This is script to process reciprocal blast data
#Max Lundberg 2021

my ($query_vs_target,$target_vs_query,$help);

GetOptions ("query_vs_target=s"	 =>\$query_vs_target,
	    "target_vs_query=s"	 =>\$target_vs_query,
	    "help"	         =>\$help	
	    );

my $USAGE = <<"USAGE";

Usage: parce_reciprocal_blast_data.new.pl --query_vs_target <file> --target_vs_query <file>

--query_vs_target	blast output file from gene models and ref proteins as query and target, respectively
--target_vs_query	blast output file from ref proteins and gene  models as query and target, respectively
--help			display this message

Note that both blast output files should be in a tabular format (fmt 6)
The name of the gene models should follow that provided by augustus and the
target file contain proteins with ensembl_id|ensemb_transcript_id|associated_gene_name

USAGE

if($help){
	die $USAGE;
	}

if(!$query_vs_target || !$target_vs_query){
	die "\nNeed to specify blast output files\n$USAGE"
	}




#######################################

#Process the query vs target blast file

#######################################


my %blast_info_query;
my %blast_target_query;
my %target_gene_transcript_overlap;


open(QUERY_BLAST,$query_vs_target) or die "\nCannot open file: $query_vs_target\n$USAGE";


while(my $blast_query=<QUERY_BLAST>){

	chomp($blast_query); 

	my ($query,$target,$perc_id,$al_length,$mismatches,$gap,$query_start,$query_end,$target_start,$target_end,$evalue,$bit_score)=split("\t",$blast_query); 

	$query=~/([^.]+)/;

	my $gene=$1;

	
	$target=~/([^:]+)\:(\S+)/;
	my $target_gene_id=$1; 
	my $target_transcript_id=$2;

		
	#store the overlap with each gene
	$target_gene_transcript_overlap{$target_gene_id}{$gene}=[($perc_id,$target_start,$target_end)];
	$blast_target_query{$gene}=$target_gene_id;
	$blast_info_query{$gene}.="$query\t$target_gene_id\t$target_transcript_id\t$target_start\t$target_end\t$perc_id\t$al_length\t$evalue\t$bit_score";


	}

close(QUERY_BLAST);


#######################################

#Process the target vs query blast file

#######################################

open(BLAST_TARGET,$target_vs_query) or die "\nCannot open file: $target_vs_query\n$USAGE";; 


my %target_query_hit;
my %target_query_hit_info;


while(my $blast_target=<BLAST_TARGET>){

	#print $blast_target;

	my ($query,$target,$perc_id,$al_length,$mismatches,$gap,$query_start,$query_end,$target_start,$target_end,$evalue,$bit_score)=split("\t",$blast_target);


	$query=~/([^:]+)/;
	my $query_id=$1; 


	$target=~/([^.]+)/;
	my $target_gene=$1;

	$target_query_hit{$query_id}=$target_gene;

	$target_query_hit_info{$query_id}="$perc_id\t$al_length";
	
	}

close(BLAST_TARGET); 




#######################################

#Summarize the data from both files

#######################################


print "gene\ttranscript_used\ttarget_name\ttarget_transcript\ttarget_start\ttarget_end\tperc_id\tal_length\tevalue\tscore\tn_hits_same_target\t";
print "reciprocity\treciprocal_gene\treciprocal_gene_perc\treciprocal_gene_aln\n";

for my $gene (keys %blast_info_query){

	print "$gene\t$blast_info_query{$gene}\t"; 


	my $target_gene_id=$blast_target_query{$gene};

	my $n_target_gene_hits=scalar(keys %{$target_gene_transcript_overlap{$target_gene_id}});

	print $n_target_gene_hits,"\t";


	

	if(exists($target_query_hit{$target_gene_id})){
		if($gene eq $target_query_hit{$target_gene_id}){
			print "RECIPROCAL\t"
			}
		else{
			print "NOT_RECIPROCAL\t";
			}
		print "$target_query_hit{$target_gene_id}\t$target_query_hit_info{$target_gene_id}\n";
		}
	else{
		print "NO HIT\tNA\tNA\tNA\n"
		}
	}


