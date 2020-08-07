#!/usr/bin/perl

use strict;
use warnings;

my $list = $ARGV[0]; #File with ID and filename
my $out = $ARGV[1];

my @inline;
my %file_list;
my @ordered_list;


#Read in all samples into an array
open (LIST, "$list") or die("ERROR: can not read file ($list): $!\n");
	while (<LIST>) {
		chomp;
		@inline = split(/\s+/);
		$file_list{$inline[0]}=$inline[1];
		push(@ordered_list,$inline[0]);
		}

my %sample;
my %counts;
my %lib;

my %oligo_id;
my $sample_ID;
my $key;

my $barcode;
my $bc_ct;
my $oligo;
my $library;

my $cur_file;

foreach $sample_ID (@ordered_list){
	$cur_file = $file_list{$sample_ID};
	print STDERR "Reading $sample_ID\n";

	open (COUNTS, "$cur_file") or die("ERROR: can not read file ($cur_file): $!\n");
	while (<COUNTS>) {
		chomp;
		@inline = split("\t");
		$barcode=$inline[0];
		$bc_ct=$inline[1];
		$oligo=$inline[2];
		$library=$inline[3]

		if($oligo !~ /\*/){
			die "Barcode & Sample combination seen twice\n" if(exists $counts{$barcode}{$sample_ID});

			$counts{$barcode}{$sample_ID}=$bc_ct;

			if(exists $oligo_id{$barcode}){
				die "Barcodes seen with different oligo IDs\n$barcode\n$oligo_id{$barcode}\n$cur_file\n$oligo\n" if($oligo_id{$barcode} ne $oligo);
			}
			else{
				$oligo_id{$barcode}=$oligo;
				$lib{$barcode}=$library;
			}
		}
	}
	close COUNTS;
}

print STDERR "Writing file\n";

open (OUT, ">$out") or die("ERROR: can not create $out: $!\n");
print OUT "Barcode\tOligo\tLibrary";
print OUT join ("\t",@ordered_list)."\n";

my $cur_bc;
my $cur_sample;

foreach $cur_bc (keys %counts)
	{
	print OUT "$cur_bc\t$oligo_id{$cur_bc}\t$lib{$cur_bc}";

	foreach $cur_sample (@ordered_list)
		{
		if(exists $counts{$cur_bc}{$cur_sample})
			{
			print OUT "\t$counts{$cur_bc}{$cur_sample}";
			}
		else
			{
			print OUT "\t0";
			}
		}
		print OUT "\n"

	}
close OUT
