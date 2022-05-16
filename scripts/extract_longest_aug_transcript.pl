#!/usr/bin/perl

use warnings;
use strict;

my %sequences;

my $seqname;

while(my $input=<STDIN>){

	chomp($input);

	if($input=~/^>(\S+)/){
		$seqname=$1;
		}

	else{

		#if sequence exists - concatenate the sequence that is found
		#on multiple rows 

		if(exists($sequences{$seqname})){
			$sequences{$seqname}=$sequences{$seqname}.=$input;
			}
		else{
			$sequences{$seqname}=$input;
			}

		}


	}
	

	
#check the sizes of transcripts

my %transcript_sizes;
my %longest_transcript;


for my $transcript (keys %sequences){

	$transcript=~/(g\d+)/;
	my $gene=$1;
	my $length=length($sequences{$transcript});
	
	#print "length is $length\n";
		
	if(exists($transcript_sizes{$gene})){
		if($transcript_sizes{$gene}<$length){
			$longest_transcript{$gene}=$transcript;
			$transcript_sizes{$gene}=$length;
			}
		}
	else{
		$transcript_sizes{$gene}=$length;
		$longest_transcript{$gene}=$transcript;
		}
	
	
	print $transcript,"\t";
	print $gene,"\t";
	print length $sequences{$transcript},"\n";
	}

for my $transcript (keys %sequences){

		$transcript=~/(g\d+)/;
		my $gene=$1;
		
		#print "longest transcript: $longest_transcript{$gene}\n";

		if($transcript eq $longest_transcript{$gene}){
			print ">$transcript\n$sequences{$transcript}\n";
			}


		}
