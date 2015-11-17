#!/bin/perl

# Cleans exonerate outputs (run with --ryo ">%ti [%tab - %tae] all%qi\n%tas\n>%ti [%tab - %tae] cds%qi\n%tcs\n")
# Adds all CDS found in a single fasta file, with distinct prefixes if from distinct proteins
#
# Usage: perl exoneratecontig.pl infile outfile type add
#
# - type = all or cds (to filter out either all / cds sequences)
# - add = 0 or 1 (add accessions to existing file / start from scratch)
#
###
use File::Basename;
use List::Util qw[min max];

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $type = $ARGV[2]; #either "all" or "cds"
my $add = $ARGV[3];
chomp($target);


chomp($file1);
chomp($ID);

## import Exonerate file
open(FILE, "$infile");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

## skip headers with parameters
shift @fields;

## load it into hash.
my %fasta; #hash with start / stop coordinates of CDS
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  my $acc = trim($tmp[0]);

  ## Get accessions info

#   print "testing $acc\n";

  if($acc =~ /(\w*)\s\[(\d*) - (\d*)\]\s(.*)\s$type/){
#     print "Import $acc\n";
    $target = $1;
    my $start = $2;
    my $stop = $3;
    my $range = abs($3 - $2);
    my $protein = $4;

    if($start > $stop){
      my @coo = ($stop, $start);
      my $tmp = $stop;
      $stop = $start;
      $start = $tmp;
      }
    ## Get sequence (clean it from end of file params)
    $_ = $tmp[1];
    s/\r|\n//g;
    s/\-\- completed exonerate analysis//g;
    my $seq = trim($_); 

    ## Store everything in CDS hash, with proteins being numbered (all cds from same protein share the same number)
    if($protein){
      $fasta{"$start"}{"$protein"}{"seq"} = $seq;
      $fasta{"$start"}{"$protein"}{"protein"} = $protein;
      $fasta{"$start"}{"$protein"}{"target"} = $target;
      $fasta{"$start"}{"$protein"}{"range"} = $range;
      $fasta{"$start"}{"$protein"}{"start"} = $start;
      $fasta{"$start"}{"$protein"}{"stop"} = $stop;
      }
    }
  }

#### Find start / stop of proteins
my %prts; #hash with start / stop coordinates of complete proteins (used for overlap test)
my $previous = "NULL";
my $prtnr = 0;
my @starts;
my @stops;

# print "Import CDS:\nProtnr\tStart\tStop\tProtein\n";
foreach $start (sort(keys(%fasta))){
  foreach $cds (keys(%{$fasta{$start}})){

    $start = $fasta{$start}{$cds}{"start"};    
    $stop = $fasta{$start}{$cds}{"stop"};

    if($cds ne $previous){ # we enter a new protein
      # open new protein
      $previous = $cds;
      $prtnr++;

      undef(@starts);
      undef(@stops);
      push @starts, $start;
      push @stops, $stop;

      } else { # we are within the same protein
      # evaluate coordinates
      push @starts, $start;
      push @stops, $stop;
      }

    # save ongoing protein in %prts
    $prts{$prtnr}{"protein"} = $cds;
    $prts{$prtnr}{"start"} = min(@starts);
    $prts{$prtnr}{"stop"} = max(@stops);
#     print "$prtnr\t$start\t$stop\t$cds\n";

    # add $prtnr in %fasta (in order to link it to %prts)
    $fasta{$start}{$cds}{"prtnr"} = $prtnr;        
    }
  }

#### FILTER PROTEINS
## Test of overlap among proteins
# print "\nOverlap test:\n";
my $ref = \%prts;
my %RES = %{ovlppw($ref)};

## Cleaning overlapping proteins (keeps longest ones)
foreach $start (sort(keys(%fasta))){
  foreach $cds (keys(%{$fasta{$start}})){
    foreach $test (keys(%RES)){
      $rmve = $RES{$test}{"rmve"};
      if($fasta{$start}{$cds}{"prtnr"} == $rmve){
	delete($fasta{$start}{$cds});
	}
      }
    }
  }

#### SPLIT CDS and save as Fasta
## Print fasta, suffix different cds of proteins with distinct codes
print "\nSave CDS:\nProtnr\tAnnotation\tSource\n";
if($add == 0){
  open(OUT, ">$outfile");
  } else {
  open(OUT, ">>$outfile");
  }
my $previous = 1e9;
my $cdsnr = 0;
foreach $cds (sort(keys(%fasta))){
  foreach $start (sort(keys(%{$fasta{$cds}}))){
    if($fasta{$cds}{$start}){
      my $protein = $fasta{$cds}{$start}{"protein"};
      if($protein){
	my $prtnr = $fasta{$cds}{$start}{"prtnr"};
	my $seq = $fasta{$cds}{$start}{"seq"};
	my $target = $fasta{$cds}{$start}{"target"};
	if($prtnr == $previous){
	  print OUT "$seq\n";
	  } else {
	  $previous = $prtnr;
	  print OUT ">$target\_PRT$prtnr\_$protein\n$seq\n";
	  }
	print "$prtnr\t$protein\t$target\n";
	}
      }
    }
  }
close(OUT);

########## SUBROUTINES
## remove whitespaces
sub trim($){
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
  }

## sub for testing CDS overlaps. 
# Input: hash containing start / end coordinates of all proteins found in focal sequence
#
# produces ref to hash
# where keys = pws comparison
# and contains keep / rmve status of each tested protein
#
###
sub ovlppw {
  my %coords = %{$_[0]};
  my @prots = keys(%coords);
  my %RES;
  my $pws = 0;
  foreach my $idx1 (0..$#prots){
    foreach my $idx2 (0..$#prots){
      if($idx2 > $idx1){
	my $prt1 = @prots[$idx1];
	my $prt2 = @prots[$idx2];

	my $s1 = $coords{$prt1}{"start"};
	my $e1 = $coords{$prt1}{"stop"};
	my $s2 = $coords{$prt2}{"start"};
	my $e2 = $coords{$prt2}{"stop"};
	
	my $test = ovlp($s1, $e1, $s2, $e2);
# 	print "test $prt1 vs $prt2: ovlp = $test\n";

	if($test == 1){
	  if(abs($e1 - $s1) >= abs($e2 - $s2)){
	    $RES{$pws}{"keep"} = $prt1;
	    $RES{$pws}{"rmve"} = $prt2;
	    } else {
	    $RES{$pws}{"keep"} = $prt2;
	    $RES{$pws}{"rmve"} = $prt1;
	    }
	  $pws++;
	  }
	}
      }
    }
  my $adrs = \%RES;
  return $adrs;
  }

# overlap test give start / end positions of two proteins
sub ovlp {
  $s1 = $_[0];
  $e1 = $_[1]; 
  $s2 = $_[2]; 
  $e2 = $_[3];
  if($s2 < $e1 && $s1 < $e2){
    return 1;
    } else {
    return 0;
    }
  }