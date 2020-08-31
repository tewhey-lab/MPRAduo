#!/usr/bin/perl

use strict;
use warnings;

my $fasta = $ARGV[0];
my $read = $ARGV[1];
my $out = $ARGV[2];
my $link_A_bc = $ARGV[3];
my $link_P_bc = $ARGV[4];
my $link_1_bc = $ARGV[5];
my $link_2_bc = $ARGV[6];
my $link_A_oligo = $ARGV[7];
my $link_P_oligo = $ARGV[8];
my $end_A_oligo = $ARGV[9];
my $end_P_oligo = $ARGV[10];
my $MIN_SEQ_SIZE = $ARGV[11];

open (FASTA, "$fasta") or die("ERROR: can not read file ($fasta): $!\n");
open (MATCH, ">$out".".match") or die("ERROR: can not create $out .matched: $!\n");
open (REJECT, ">$out".".reject") or die("ERROR: can not create $out .rejected: $!\n");

# my $link_A_bc = "TCTAGA";
# my $link_P_bc = "AGTCAG";
# my $link_1_bc = "CTAG";
# my $link_2_bc = "CTCA";
# my $link_A_oligo = "AGTG";
# my $link_P_oligo = "AGGG";
# my $end_A_oligo = "CGTC";
# my $end_P_oligo = "GCAA";
# my $MIN_SEQ_SIZE = 100;

my $barcode_seq;
my $oligo_seq;
my $barcode_start;
my $oligo_start;
my $oligo_end;
my $oligo_length;
my $umi_start;
my $id;
my $r1;
my $revcomp;
my $revcomp_oligo;
my $revcomp_barcode;
my $link_index;
my $library = "";
my $umi;

while (<FASTA>){
# Extract Sequence ID
  chomp;
  $id = $_;
  $id =~ s/^@/>/;
  $id =~ s/^>//;
  $id =~ s/\/1$//;

# Extract the sequence
  $r1 = <FASTA>;
  chomp $r1;
  if(length($r1) < $MIN_SEQ_SIZE | $id eq "+"){
    next;
  }

# If checking a Read 1 file take the reverse complement
  if($read == 1){
    $revcomp = reverse($r1);
    $revcomp =~ tr/ACGTNacgtn/TGCANtgcan/;
    $r1 = $revcomp;
  }
# Check for presence of linker P
  if(index(substr($r1, 18, 10), $link_P_bc) != -1){
    $link_index = index(substr($r1, 18, 10), $link_P_bc);
    $link_index += 18;
    $library = "P";
  }
# Check for presence of linker A
  elsif(index(substr($r1, 8, 10), $link_A_bc) != -1){
    $link_index = index(substr($r1, 18, 10), $link_A_bc);
    $link_index += 8;
    $library = "A";
  }
# Check for presence of linker P for Duo
  elsif(index(substr($r1, 40, 12), $link_P_bc) != -1){
    $link_index = index(substr($r1, 40, 12), $link_P_bc);
    $link_index += 40;
    $library = "P";
  }
# Check for presence of linker A for Duo
  elsif(index(substr($r1, 30, 12), $link_A_bc) != -1){
    $link_index = index(substr($r1, 30, 12), $link_A_bc);
    $link_index += 30;
    $library = "A";
  }
  elsif($library eq ""){
    print REJECT "$id\n";
    next;
  }

# Match libraries to barcodes
  if($library eq "P"){
    if($link_index - 20 < 0){
      $barcode_start = 0;
      $barcode_seq = substr($r1, $barcode_start, 20);
    }
    if($link_index - 20 >= 0){
      $barcode_start = $link_index - 20;
      $barcode_seq = substr($r1, $barcode_start, 20);
      $umi_start = index(substr($r1, 4, 10), $link_1_bc) + 4;
      $umi = substr($r1, 0, $umi_start);
    }
  }
  if($library eq "A"){
    if($link_index - 10 < 0){
      $barcode_start = 0;
      $barcode_seq = substr($r1, $barcode_start, 10);
    }
    if($link_index - 10 >= 0){
      $barcode_start = $link_index - 10;
      $barcode_seq = substr($r1, $barcode_start, 10);
      $umi_start = index(substr($r1, 4, 10), $link_2_bc) + 4;
      $umi = substr($r1, 0, $umi_start);
    }
  }

# Find start of oligo
  if($library eq "P"){
    $oligo_start = index(substr($r1, $link_index+30, 12), $link_P_oligo);
    if($oligo_start == -1){
      $oligo_start = 0;
    }
    $oligo_start += 34+$link_index;
  }
  if($library eq "A"){
    $oligo_start = index(substr($r1, $link_index+30, 19), $link_A_oligo);
    if($oligo_start == -1){
      $oligo_start = 0;
    }
    $oligo_start += 34+$link_index;
  }
# Find the end of the oligo
  if($library eq "P"){
    $oligo_end = index(substr($r1, -18, 12), $end_P_oligo);
    $oligo_end += -18;
  }
  if($library eq "A"){
    $oligo_end = index(substr($r1, -17, 12), $end_A_oligo);
    $oligo_end += -17;
  }

# Define the substring that is the oligo
  $revcomp_barcode = reverse($barcode_seq);
  $revcomp_barcode =~ tr/ACGTNacgtn/TGCANtgcan/;
  $barcode_seq = $revcomp_barcode;
  $oligo_length = length($r1) + $oligo_end - $oligo_start;
  $oligo_seq = substr($r1, $oligo_start, $oligo_length);
  $revcomp_oligo = reverse($oligo_seq);
  $revcomp_oligo =~ tr/ACGTNacgtn/TGCANtgcan/;
  $oligo_seq = $revcomp_oligo;
  $umi = $umi . ('B' x (8-length($umi)));
  print MATCH join("\t", $id, $barcode_seq, $oligo_seq, $oligo_length, $umi."\n");
}
