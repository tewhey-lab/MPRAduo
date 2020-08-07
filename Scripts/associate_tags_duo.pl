#!/usr/bin/perl

##################
#
#052314 - Updated to account for multimapping oligos
#
####
#	FLAGS
##################

use strict;
use warnings;


my $tags = $ARGV[0]; #matched tag file
my $out = $ARGV[1];
my $read = $ARGV[2];

open (OUT, ">$out") or die("ERROR: can not create $out: $!\n");

my @inline;
my %tags; #[ct,orientation,location,tagseq(rc)]
my $revcomp;


open (TAGS, "$tags") or die("ERROR: can not read file ($tags): $!\n");
while (<TAGS>){

	chomp;
	@inline = split("\t");

  ${$tags{$inline[1]}}[0]++; #count
	${$tags{$inline[1]}}[2] = $inline[2]; #oligo
	${$tags{$inline[1]}}[4] = $inline[3]; #library assingment

}
close TAGS;

my $key;

foreach $key (keys %tags)
	{
	print OUT join("\t",$key,@{$tags{$key}})."\n";
	}
