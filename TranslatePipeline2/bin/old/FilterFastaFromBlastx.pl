#!/bin/perl

#####################################
#### Get specific sequences out of fasta file, using headers from headersfile
####
####  perl FilterFasta.pl fastafile headersfile keepdoublons
#####################################

my $file1 = $ARGV[0];
my $ID = $ARGV[1];
my $doublons = 0;

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
  @acc = split(/\s/, $tmp[0], 2);
  $acc = trim($acc[0]);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = trim($_);
  $fasta{$acc} = $seq; 
  }

## import targets
open(HD, "$ID");
my @ids;
my %idshs;
while (<HD>){
  chomp();
  my @tmp = split(/\s/, $_, 2);
  my $acc = $tmp[1];
  push (@ids, $acc);
  $idshs{$acc} = 1;
  }
close HD;

if($doublons == 0) {
  @ids = keys(%idshs);
  } else {
  @ids = @ids;
  }

## print sequences of interest into output
open(OUT, ">$file1.filter");
foreach $head (@ids){
  chomp($head);
  my $seq = $fasta{$head};
  print(OUT ">$head\n$seq\n");
  }
close(OUT);


## remove whitespaces
sub trim($){
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
  }