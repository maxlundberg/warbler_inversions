#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

#Max Lundberg 2017

my ($intersection_file,$help); 


GetOptions ("intersection_file=s" =>\$intersection_file,
	    "help"	          =>\$help,	
	    );

my $USAGE = <<"USAGE";

Usage: summarize_gene_models_synteny_intersect.pl --intersection_file

--intersection_file		intersection between a gene model file (augustus gff)
				and syntenty-transferred ensembl gene models from kraken 
--help				display this message

USAGE

if($help){
	die $USAGE;
	}

if(!$intersection_file){die "\nNeed to specify an intersection file\n$USAGE"}

my (%cds_length, %cds_number, %overlap_length, %overlap_number, %overlap_start_exon, %overlap_end_exon); 

my %transcript_gene_id;
my %ensembl_gene_name;
my %kraken_genes;
my %exon_names;
my %target_gene_name;

my %scaffold;
 


open(INTERSECTION,$intersection_file) or die "Cannot open the intersection file\n"; 


while(my $input=<INTERSECTION>){

	chomp($input);

	my @input=split("\t",$input);

	my $scaffold=$input[0];
	my ($cds_start,$cds_end)=@input[3..4];
	my $length=$cds_end-$cds_start+1;	
	my $gene_info=$input[8];
	$gene_info=~/Parent=(.+)/;
	my $transcript_id=$1;

	my $overlap=$input[18];

	my $target_attributes=$input[17];


	$scaffold{$transcript_id}=$scaffold;


	#create a unique exon (cds) name from the positional information
	#summarize the number of cds exons and the total amount of CDS sequence for each
	#gene model. Since each transcript from the model gene could map to several
	#transcripts, we should be careful not to double-count 

	#Note that this should be restricted to each transcript	


	my $exon_name=join("-",$scaffold,$cds_start,$cds_end);

	if(!exists($exon_names{$transcript_id}{$exon_name})){
	
		if(exists($cds_length{$transcript_id})){
			$cds_length{$transcript_id}+=$length;	
			$cds_number{$transcript_id}++;
			}
		else{
			$cds_length{$transcript_id}=$length;
			$cds_number{$transcript_id}=1;
			}
	
		$exon_names{$transcript_id}{$exon_name}=1;		
		}


#	if(!exists(${$exon_name})){
#	
#		if(exists($cds_length{$transcript_id})){
#			$cds_length{$transcript_id}+=$length;	
#			$cds_number{$transcript_id}++;
#			}
#		else{
#			$cds_length{$transcript_id}=$length;
#			$cds_number{$transcript_id}=1;
#			}
#	
#		$exon_names{$exon_name}=1;		
#		}




	#if there is an overlap with the kraken genes
	
	if($overlap>0){	
	
	
		#extract target gene information from the kraken attribute field
		$target_attributes=~/gene_id "([^;]+)"\; transcript_id "([^;]+)"/;
		my $target_gene_id=$1;
		my $target_transcript_id=$2;
		$target_attributes=~/gene_name "([^;]+)"/;
		my $target_gene_name=$1;
		$target_gene_name=~tr/[a-z]/[A-Z]/;

		$target_gene_name{$target_gene_id}=$target_gene_name;
		
		###########################
		
		#count the number of times a specific kraken gene is overlapping with the gene model transcript
				
		if(exists($kraken_genes{$transcript_id}{$target_gene_id})){
			$kraken_genes{$transcript_id}{$target_gene_id}++;
			}
		else{
			$kraken_genes{$transcript_id}{$target_gene_id}=1;
			}	

		###########################	
			
		#also add statistics about the kraken transcript it is overlapping with 
		#cumulative length of overlap, number of exons and start-end exons
		
		#if there is any previously recorded overlap with the transcript (also implies that other data
		#about the transcript should be available 
		
		if(exists($overlap_length{$transcript_id}{$target_transcript_id})){
			
			$overlap_length{$transcript_id}{$target_transcript_id}+=$overlap;
			$overlap_number{$transcript_id}{$target_transcript_id}++;


			#if($overlap_start_exon{$transcript_id}{$target_transcript_id}>$cds_number{$transcript_id}){
			#	$overlap_start_exon{$transcript_id}{$target_transcript_id}=$cds_number{$transcript_i}
			#	}
			#if($overlap_end_exon{$transcript_id}{$target_transcript_id}<$cds_number{$transcript_id}){$overlap_end_exon{$transcript_id}{$target_transcript_id}=$cds_number{$transcript_id}}
			}
			
			
		else{
			$overlap_length{$transcript_id}{$target_transcript_id}=$overlap;
			$overlap_number{$transcript_id}{$target_transcript_id}=1; 
			$overlap_start_exon{$transcript_id}{$target_transcript_id}=$cds_number{$transcript_id};
			$overlap_end_exon{$transcript_id}{$target_transcript_id}=$cds_number{$transcript_id};		

			}
		
		} #end of kraken gene overlap


	}
	

close(INTERSECTION); 


	
###################################

#PROCESS THE OVERLAP OUTPUT

###################################


#print a table header

print "gene\ttranscript\tscaffold\tcds_length\tcds_number\tnumber_genes\tgene_ensembl_id\ttranscript_number\ttranscript_id\toverlap_length\tperc_overlap\tcds_overlap\n";


for my $transcript (keys %cds_length){

	#in case of augustus gene models
	$transcript=~/([^.]+)/; 
	
	#in case of evm-based models 
	#$transcript=~/(evm\.model\.[^.]+\.\d+)/;	
	
	my $gene_id=$1;

	my $scaffold=$scaffold{$transcript};

	print "$gene_id\t$transcript\t$scaffold\t$cds_length{$transcript}\t$cds_number{$transcript}\t";
	
	my $kraken_hits=scalar(keys %{$overlap_length{$transcript}});
	my $gene_hits=scalar(keys %{$kraken_genes{$transcript}});


	
	if($kraken_hits>0){

		my $kraken_genes;
		my $kraken_gene_names;

		my $kraken_transcripts;
		my $overlap_lengths;
		my $overlap_numbers;
		my $perc_overlaps;
		my $exon_bounds;

		my $gene_number=1;


		for my $gene_id (keys %{$kraken_genes{$transcript}}){

			#print "gene id is $gene_id\n";
		
			if($gene_number==1){
				$kraken_genes=$gene_id;
				$kraken_gene_names=$target_gene_name{$gene_id};
				}
			else{
				$kraken_genes.=",$gene_id";
				$kraken_gene_names.=",$target_gene_name{$gene_id}";
				}

			$gene_number++; 

			}




		my $transcript_number=1;


		for my $kraken_transcript (keys %{$overlap_length{$transcript}}){
			
			my $overlap_length=$overlap_length{$transcript}{$kraken_transcript}; 
			my $overlap_number=$overlap_number{$transcript}{$kraken_transcript};
			my $perc_overlap=sprintf("%.1f",$overlap_length/$cds_length{$transcript}*100);
			#my $start_exon=$overlap_start_exon{$transcript}{$kraken_transcript};
			#my $end_exon=$overlap_end_exon{$transcript}{$kraken_transcript};

			#my $kraken_genes; 
			
			#my $start_exon=1; my $end_exon=1;

			if($transcript_number==1){
				$kraken_transcripts=$kraken_transcript;
				$overlap_lengths=$overlap_length; 
				$overlap_numbers=$overlap_number;
				$perc_overlaps=$perc_overlap;
				#$exon_bounds="$start_exon-$end_exon";
				}	

			else{
				$kraken_transcripts.=",$kraken_transcript";
				$overlap_lengths.=",$overlap_length";
				$overlap_numbers.=",$overlap_number";
				$perc_overlaps.=",$perc_overlap";
				#$exon_bounds.=",$start_exon-$end_exon";
				#$end_exons.=",$end_exon";
				}


			#if($kraken_hits==1){print "$kraken_transcript\t$overlap_length\t$perc_overlap\t$overlap_number\t$start_exon-$end_exon"}
			#else{
			#print "$kraken_transcript:$overlap_length:$perc_overlap:$overlap_number:$start_exon-$end_exon;";
			#}
			$transcript_number++;
	
			}

			print "$gene_hits\t$kraken_genes\t$kraken_gene_names\t$kraken_hits\t$kraken_transcripts\t$overlap_lengths\t$perc_overlaps\t$overlap_numbers\n"
		}

	#if no gene overlap

	else{
		print "0\tno_overlap\tNA\tNA\tNA\tNA\tNA\tNA\n"
		}


	#if(exists($blastinfo{$transcript})){
	#	print $blastinfo{$transcript};
	#	}
	#else{
	#	print "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"		
		#}
		
	}
	

