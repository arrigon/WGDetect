#!/bin/perl

#####################################
#### Add a given prefix to all headers in a fasta file
####
####  perl PrefixFasta.pl fastafile prefix
#####################################

my $file1 = $ARGV[0];
my $outfile = $ARGV[1];

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
my $nr = 0;
my @speclist;
my $nbp;

foreach $input (@fields){
  chomp($input);
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r*|\n*//g;
  s/!/N/g;
  $seq = $_;
  
  $_ = $tmp[0];
  s/\r*|\n*//g;
  $acc = $_;
  
  $fasta{$acc} = $seq;
  push(@speclist, $acc);
  $nr++;
  }


## print sequences of interest into output
my $nseq = $#speclist + 1;
my $nbp = length($seq);

open(OUT, ">$outfile");
print(OUT "$nseq $nbp\n");
foreach $head (@speclist){
  chomp($head);
  my $seq = $fasta{$head};
  print(OUT "$head     $seq\n");
  }
close(OUT);

