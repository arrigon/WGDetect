#!/bin/perl

#####################################
#### Add a given prefix to all headers in a fasta file
####
####  perl PrefixFasta.pl fastafile prefix
#####################################

my $file1 = $ARGV[0];
my $prefix = $ARGV[1];

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
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  $fasta{$tmp[0]} = $seq; 
  }


## print sequences of interest into output
open(OUT, ">cnt\_$file1");
my $cnt = 0;
foreach $head (keys(%fasta)){
  chomp($head);
  my $seq = $fasta{$head};
  print(OUT ">$cnt\n$seq\n");
  $cnt++;
  }
close(OUT);

