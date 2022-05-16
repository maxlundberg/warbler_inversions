#!/usr/bin/perl

use warnings;
use strict;
use List::Util qw(max);
use List::MoreUtils qw(uniq);
use Getopt::Long;

#Max Lundberg 2021

#This is script to combine annotation information

my ($gene_info_file,$blast_summary_file,$blast_target_info_file,$kraken_intersect_file,$swissprot_blast,$swissprot_info_file,$interpro_summary,$help);

GetOptions ("gene_info_file=s"		 =>\$gene_info_file,
	    "blast_summary_file=s"	 =>\$blast_summary_file,
	    "blast_target_info_file=s"	 =>\$blast_target_info_file, 	
	    "kraken_intersect_file=s"	 =>\$kraken_intersect_file,
	    "swissprot_blast=s"		 =>\$swissprot_blast,
	    "swissprot_info_file=s"	 =>\$swissprot_info_file,	 
	    "interpro_summary=s"	 =>\$interpro_summary,			 	
	    "help"	         	 =>\$help	
	    );


my $USAGE = <<"USAGE";

Usage: summarize_annotation_data.pl --gene_info_file <file> --blast_summary_file <file> --kraken_intersect_file <file> --swissprot_blast <file> --swissprot_info_file <file> 	--interpro output summary <file>


--gene_info_file	table with statistics of genes, one gene per row. The contents of this file
			will not be modified, but added to the first columns of the output
--blast_summary_file	summary file of blast reciprocal data
--blast_target_info_file	info about target genes in the reciprocal blast analysis
--kraken_intersect_file	summary file of intersection of gene model transcripts and kraken synteny-transferred models
--swissprot_blast	swissprot blast results file
--swissprot_info_file	file containing metadata about each swissprot entry 
--interpro_summary	interpro output summary


--help			display this message

USAGE

if($help){
	die $USAGE;
	}

if(!$gene_info_file || !$blast_summary_file || !$blast_target_info_file || !$kraken_intersect_file || !$swissprot_blast || !$swissprot_info_file || !$interpro_summary){
	die "\nNeed to specify input files\n$USAGE"
	}


$,="\t";


#############################################

#extract information from the reciprocal blast

#############################################

my %blast_summary;

open(BLAST_SUMMARY,$blast_summary_file) or die "\ncannot open blast summary file $blast_summary_file\n$USAGE\n"; 

while(my $blast_input=<BLAST_SUMMARY>){
	
	chomp($blast_input);

	my @blast_input=split("\t",$blast_input);

	my $gene_model=$blast_input[0];
	
	$blast_summary{$gene_model}=[@blast_input[1..$#blast_input]];	
	

	}

close(BLAST_SUMMARY);


#############################################

#extract information about the reciprocal blast
#target genes 

#############################################

my %reciprocal_blast_info; 

open(BLAST_INFO,$blast_target_info_file) or die ""; 

while(my $blast_info_input=<BLAST_INFO>){

	chomp($blast_info_input);

	my @blast_info_input=split("\t",$blast_info_input);

	my $prot_id=$blast_info_input[1];

	$reciprocal_blast_info{$prot_id}=[@blast_info_input[0,2..4]];

	}

close(BLAST_INFO);

#############################################

#Read metadata about swissprot protein

#############################################

#my $swissprot_info_file="swissprot_vertebrates_20190627_info.txt";

my %swissprot_info;


open(SWISSPROT_INFO,$swissprot_info_file) or die "cannot open swissprot info file $swissprot_info_file\n"; 

while(my $swissprot_info=<SWISSPROT_INFO>){

	next if $.==1;

	chomp($swissprot_info); 

	my ($entry,$entry_name,$status,$protein_names,$gene_names,$organism,$length)=split("\t",$swissprot_info); 

	#select the first gene name 
	
	my $gene_name;

	if($gene_names=~/(\S+)/){
		$gene_name=$1;
		$gene_name=~tr/[a-z]/[A-Z/;
		}
	else{
		$gene_name="uncharacterized";
		}

	

	$swissprot_info{$entry}=$gene_name

	}

close(SWISSPROT_INFO);

#############################################

#extract information from the swissprot blast

#############################################

open(SWISSPROT,$swissprot_blast) or die "\ncannot open the swissprot blast file $swissprot_blast\n$USAGE\n";

my %swissprot_blast_data;


while(my $swissprot_blast=<SWISSPROT>){

	chomp($swissprot_blast);

	my ($qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore)=split("\t",$swissprot_blast);

	$qseqid=~/([^.]+)/;
		
	my $gene=$1;

	#print "sseqid $sseqid\n";

	$sseqid=~/sp\|([^|]+)/;

	#print $1,"\t",$swissprot_info{$1},"\n";

	my $swissprot_id=$1;

	$swissprot_blast_data{$gene}=join(",",$swissprot_info{$swissprot_id},$sseqid,$pident,$length,$evalue,$bitscore);	
	
	
	}

close(SWISSPROT);


#############################################

#extract information from the kraken data

#############################################

my %best_overlap_transcript;
my %best_overlap_gene; 
my %best_overlap_perc;

my %target_genes;
my %target_transcripts;

open(KRAKEN_INTERSECT,$kraken_intersect_file) or die "cannot open kraken intersect file $kraken_intersect_file\n"; 

while(my $kraken_input=<KRAKEN_INTERSECT>){

	next if($.==1);

	chomp($kraken_input);

	my ($gene,$transcript,$scaffold,$cds_length,$cds_number,$n_overlap_genes,$overlap_id,$overlap_name,$overlap_transcripts,
	$transcripts,$overlap_transcript_cds_length,$overlap_transcript_cds_perc,$overlap_cds_count)=split("\t",$kraken_input);  	


	my @overlap_transcript_cds_perc=split(",",$overlap_transcript_cds_perc);

	my $max_overlap=0;
	my @genes;
	my @transcripts;

	if($n_overlap_genes>0){
		$max_overlap=max(@overlap_transcript_cds_perc);

		
		#add the genes and transcripts to a common list

		if(exists($target_genes{$gene})){
			$target_genes{$gene}.=",$overlap_id";
			$target_transcripts{$gene}.=",$transcripts";
			}
		else{
			$target_genes{$gene}=$overlap_id;
			$target_transcripts{$gene}=$transcripts;
			}		
		
			}

	#save the identity of the transcript with the best overlap, which gene
	#it overlaps with and what percentage the overlap is	
	
	if(exists($best_overlap_transcript{$gene})){
		if($best_overlap_perc{$gene}<$max_overlap){
			$best_overlap_perc{$gene}=$max_overlap;				
			$best_overlap_gene{$gene}="$overlap_id\t$overlap_name";
			$best_overlap_transcript{$gene}=$transcript;
			}
		
		}
	else{
		$best_overlap_transcript{$gene}=$transcript;
		$best_overlap_gene{$gene}="$overlap_id\t$overlap_name";
		$best_overlap_perc{$gene}=$max_overlap; 		
		}





	#print $input


	}

close(KRAKEN_INTERSECT); 



#############################################

#extract information from interpro

#############################################

open(INTERPRO,$interpro_summary) or die "\ncannot open the interpro summary file $interpro_summary\n$USAGE\n";

my %interpro_data;


while(my $interpro_input=<INTERPRO>){

	next if ($.==1);

	chomp($interpro_input);

	my ($gene,$annotation)=split("\t",$interpro_input);

	$interpro_data{$gene}=$annotation;

	}

close(INTERPRO);




#############################################

#Open gene info data and print the combined data

#############################################



#print "gene_model\tscaffold\tstart\tend\tstrand\tn_transcripts\tbest_transcript_hints_support_perc\ttranscript_hints_support\tcds_number\tcds_length\t";
my $new_column_heads="annotation\tannotation_source\tnew_gene_name\twarnings\t";
$new_column_heads.="kraken_gene_id\tkraken_gene_name\tkraken_best_overlap\tkraken_n_genes_overlap\tkraken_all_genes_overlap\t";
$new_column_heads.="blast_gene_name\tblast_gene_transcript\tblast_description\tblast_speceis\ttarget_start\ttarget_end\tperc_id\talign\tevalue\tbitscore\tn_genes_same_blast\treciprocity\ttarget_hit\ttarget_perc\ttarget_aln\t";
$new_column_heads.="swissprot_gene_name\tswissprot_id\tswissprot_perc\tswissprot_length\tswissprot_evalue\tswissport_score\tinterpro_domains\n";



#initialize a gene count variable
my $gene_count=1;


open(GENE_INFO,$gene_info_file) or die "cannot open $gene_info_file\n"; 

while(my $gene_info=<GENE_INFO>){

	chomp($gene_info); 

	if($.==1){
		print "$gene_info\t$new_column_heads";
		next;	
 		}

	my @gene_info=split("\t",$gene_info); 

	my $gene=$gene_info[0]; 

	print "$gene_info\t";



	#extract the reciprocal blast data for the gene model
	my @recip_blast_info;

	if(exists($blast_summary{$gene})){
		@recip_blast_info=@{$blast_summary{$gene}};
		}
	else{
		@recip_blast_info=("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA")
		}

	my ($blast_query_transcript,$blast_target_gene_name,$blast_target_transcript,$blast_target_start,$blast_target_end,
	   $perc_id,$aln_length,$evalue,$score,$n_hits_same_target,$reciprocity,$reciprocal_gene,$reciprocal_perc,$reciprocal_aln_length)=@recip_blast_info; 


	
	#extract the kraken data for the gene model 

	my ($kraken_gene_id,$kraken_gene_name)=split("\t",$best_overlap_gene{$gene});		


	my @genes;
	my $gene_number=0;
	my $unique_genes="NA"; 

	if(exists($target_genes{$gene})){
	
		@genes=split(",",$target_genes{$gene});
		$gene_number=scalar(uniq(@genes));
		$unique_genes=join(",",uniq(@genes));
		
		}
	
	#extract the swissprot data for the gene model

	my @swissprot_info;	

	if(exists($swissprot_blast_data{$gene})){
		@swissprot_info=split(",",$swissprot_blast_data{$gene});
		}
	else{
		@swissprot_info=("NA","NA","NA","NA","NA","NA")
		}


	my ($swissprot_gene,$swissprot_id,$swissprot_perc,$swissprot_length,$swissprot_evalue,$swissport_score)=@swissprot_info;

#	print "\n$swissprot_gene,$swissprot_perc,$swissprot_length,$swissprot_evalue,$swissport_score\n";



	#provide an annotation of the gene model and note the source
	#kraken > recip_blast > swissprot



	my ($gene_annotation,$annotation_source);

	
		if($kraken_gene_id ne "no_overlap" && $kraken_gene_id!~/^LOC/){
			$gene_annotation=$kraken_gene_name;
			$annotation_source="kraken";
			}
		else{
			if($blast_target_gene_name ne "NA" && $blast_target_gene_name!~/^LOC/){
				$gene_annotation=$blast_target_gene_name;
				$annotation_source="blast_reciprocal";
				}
			else{
				
				if($swissprot_gene ne "NA"){
					$gene_annotation=$swissprot_gene;
					$annotation_source="blast_swissprot";
					}
				else{
					$gene_annotation="unchar";
					$annotation_source="NA";
					}
				}
			}

	
	#create a new gene_id for the model

	my $gene_name="Phtr_g${gene_count}_${gene_annotation}";

	print "$gene_annotation\t$annotation_source\t$gene_name\t";



	#check for potential problems
	my $warnings;
	
	#if($model_type eq "aug"){
		if($gene_info[6]<50){$warnings.="low hint support ($gene_info[6]);"}
		#}	
	if($kraken_gene_name ne "NA" && $blast_target_gene_name){
		if($kraken_gene_name ne $blast_target_gene_name){
			$warnings.="different gene names kraken and blast;"
			}	
		}	
	if($annotation_source eq "blast" && $perc_id<50){$warnings.="low perc id to blast annotation ($perc_id);"}
	if($annotation_source eq "blast" && $reciprocity ne "RECIPROCAL"){$warnings.="not reciprocal blast hit;"}
	if($annotation_source eq "kraken" && $best_overlap_perc{$gene}<50){$warnings.="low overlap with kraken gene ($best_overlap_perc{$gene});"}
	if($annotation_source eq "kraken" && $gene_number>1){$warnings.="more than 1 kraken gene overlapping;"}
	
	if(!$warnings){
			$warnings="no warnings"
			}


	print "$warnings\t"; 
		


	print "$best_overlap_gene{$gene}\t$best_overlap_perc{$gene}\t$gene_number\t$unique_genes\t";
	print "$blast_target_gene_name\t$blast_target_transcript\t";
	
	if(exists($reciprocal_blast_info{$blast_target_transcript})){

		my ($species,$name,$symbol,$length) =@{$reciprocal_blast_info{$blast_target_transcript}};
	

		print "$name\t$species\t"	
		}
	else{
		print "NA\tNA\t"
		}

	print "$blast_target_start\t$blast_target_end\t$perc_id\t$aln_length\t$evalue\t$score\t$n_hits_same_target\t$reciprocity\t$reciprocal_gene\t$reciprocal_perc\t$reciprocal_aln_length\t";

	print "$swissprot_gene\t$swissprot_id\t$swissprot_perc\t$swissprot_length\t$swissprot_evalue\t$swissport_score\t";

	if(exists($interpro_data{$gene})){
		print $interpro_data{$gene};
		}
	else{
		print "NA"
		}

	print "\n";	

	$gene_count++;
	}


