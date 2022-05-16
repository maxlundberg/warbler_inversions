#!/usr/bin/perl

use warnings;
use strict;

#This script will extract the gap information for contigs included in the supercontigs

#Usage: cat <file.gap> | get_gaps_in_supercontigs_from_bionano_scaffolding.pl supercontig.list


my %supercontig_info; 

open(SUPERCONTIGS,$ARGV[0]) or die "cannot open supercontig file\n";

while(my $supercontig_input=<SUPERCONTIGS>){

	chomp($supercontig_input); 

	my ($super_contig_id,$super_scaffold,$start,$end,$contig)=split("\t",$supercontig_input); 

	$supercontig_info{$contig}=$super_contig_id

	}

close(SUPERCONTIGS);


print "Supercontig\tNGSId1\tNGSId2\tSuperScaffoldId\tXmapGapLength\tAdjustedGapLength\tNGSLength1\tNGSLength2\n";

while(my $gap_input=<STDIN>){

	my ($NGSId1,$NGSId2,$SuperScaffoldId,$XmapGapLength,$AdjustedGapLength,$NGSLength1,$NGSLength2)=split("\t",$gap_input); 

	
	if(exists($supercontig_info{$NGSId1}) && exists($supercontig_info{$NGSId2})){

		if($supercontig_info{$NGSId1} eq $supercontig_info{$NGSId2}){


			if($NGSId1 ne $NGSId2){

				print $supercontig_info{$NGSId1},"\t",$gap_input

				}
			}
		}
	}


