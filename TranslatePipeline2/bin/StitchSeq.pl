#!/bin/perl

# Concatenates distinct proteins, found in a same accession into a unique sequence
#
# Usage: perl StitchSeq.pl Infile Outfile
#
# WARNING: To be used only for genbank sequences already clustered into "gene families"
#          should not be used for sequences that WILL be clustered, where it is of interest to keep
#            distinct proteins as separate sequences.
## Nils Arrigo, UNil 2015
##

use File::Basename;

my $file1 = $ARGV[0];
my $bsn = basename($file1);
my $outfolder = $ARGV[1];

chomp($file1);
chomp($ID);

## import fasta file
open(FILE, "$file1");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

## load it into hash. Note that %fasta contains your sequences. Can be reused elsewhere.
my %fasta;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $acc = trim($tmp[0]);
 
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = trim($_); 

  if($seq){
    $fasta{"OK"}{$acc} = $seq;
    } else {
    $fasta{"NA"}{$acc} = 0;
    }
  }

my @content = keys(%fasta);
my $nacc = @content;
if($nacc == 0){
  die "Skipped file: $file1 is empty\n";
  }
if($content[0] =~ /^$/){
  die "Skipped file: $file1 is empty\n";
  }

## print sequences of interest into output
open(OUT, ">$outfolder/$bsn");
foreach $acc (keys(%{$fasta{"OK"}})){
  chomp($acc);
  my $seq = $fasta{"OK"}{$acc};
  print(OUT ">$acc\n$seq\n"); 
  }
close(OUT);

## Report lost sequences
open(OUT, ">$outfolder/$bsn.lost");
foreach $acc (keys(%{$fasta{"NA"}})){
  chomp($acc);
  print(OUT "$acc\tNA\n"); 
  }
close(OUT);


## remove whitespaces
sub trim($){
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
  }