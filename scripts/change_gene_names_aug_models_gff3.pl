#!/usr/bin/perl

use warnings; 
use strict;
use Getopt::Long;


#This is a script to change the names of augustus gene models according to a
#table with new and old names
 
#Max Lundberg 2017

my ($aug_gff,$name_list,$help);

GetOptions ("aug_gff=s"	 	 =>\$aug_gff,
	    "name_list=s"        =>\$name_list,    
            "help"	         =>\$help	
	    );

my $USAGE = <<"USAGE";

Usage: change_gene_names_aug_models_gff3.pl --aug_gff <aug_predictions.gff3> --name_list <old_new_names.tbl>

--aug_gff		augustus predictions in gff3 format
--name_list		two column table with old augustus name and new name for each gene
--help			display this message

USAGE

if($help){
	die $USAGE;
	}


if(!$aug_gff || !$name_list){die "\nNeed to specify input files\n$USAGE\n"}

 


#############################

#extract information from name list file

#############################


open(NAME_LIST,$name_list) or die "\ncannot open name list file: $name_list\n$USAGE\n"; 


my %names;

while(my $input=<NAME_LIST>){

	next if($input=~/^gene/);

	chomp($input); 
	my ($old_name,$new_name)=split("\t",$input);
	$names{$old_name}=$new_name;
	}


close(NAME_LIST);

############################



open(AUG_GFF,$aug_gff) or die "\ncannot open augustus gff3 file: $aug_gff\n$USAGE\n"; 



while(my $gff3_input=<AUG_GFF>){

	print "$gff3_input" if($.==1);

	next if($gff3_input=~/^#/);
	
	chomp($gff3_input);
	

	my $old_gene_name;
	my $feat_name;
	
	if($gff3_input=~/ID\=([^.]+)/){
		$old_gene_name=$1;
		}
	else{
		$gff3_input=~/Parent\=([^.]+)/;
		$old_gene_name=$1
		}

	
	if(exists($names{$old_gene_name})){
		$gff3_input=~s/(?<=\=)$old_gene_name/$names{$old_gene_name}/g;
		$gff3_input=~s/\,/\-/g;
		$gff3_input=~s/\btranscript\b/mRNA/;
		
		if($gff3_input=~/ID=([^;]+)/){
			$feat_name=$1
			}
		else{
			$gff3_input=~/Parent\=([^;]+)/;	
			$feat_name=$1
			}
		

		print $gff3_input,";Name=$feat_name","\n";
		}


	}

close(AUG_GFF);

