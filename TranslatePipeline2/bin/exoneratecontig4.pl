#!/bin/perl

# Cleans exonerate outputs (run with --ryo ">%ti [%tab - %tae] all%qi\n%tas\n>%ti [%tab - %tae] cds%qi\n%tcs\n")
# Adds all CDS found in a single fasta file, with distinct prefixes if from distinct proteins
#
# Usage: perl exoneratecontig.pl infile outfile type add
#
# - type = all or cds (to filter out either all / cds sequences)
# - add = 0 or 1 (add accessions to existing file / start from scratch)
#
# Nils Arrigo, Unil 2015
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

### skip headers with parameters
shift @fields;

### load data into hash.
my %fasta; #hash with start / stop coordinates of CDS
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  my $acc = trim($tmp[0]);

  ## Get accessions info
  if($acc =~ /(\w*)\s\[(\d*) - (\d*)\]\s(.*)\s$type/){
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


### Organise CDS as non-overlapping stacks
my @allstarts = keys(%fasta);
@allstarts = sort{$a <=> $b} @allstarts;

# Second, visit hash of CDS, following that order 
# define non-overlapping stacks; and keep longest member of each stack.
my $stacknr = 0; 
my $rightlim = 0; #rightmost extent of ongoing stack
my %prtsII;
my $maxrge = 0;
my @keepprot; #array of protein numbers that will be kept after parsing

# Attribute CDS within stacks
my %fastaII;
my $cdsnr = 0;
# print "\nDefining Stacks:\nCDSnr\tSTACKnr\tStart\tStop\tRightLim\tProtein\n";
foreach $start (@allstarts){ #visit CDS hash following increasing order of start pos
  foreach $cds (keys(%{$fasta{$start}})){
    my %tmp = %{$fasta{$start}{$cds}};
    my $protein = $tmp{"protein"};
    my $start = $tmp{"start"};
    my $stop = $tmp{"stop"};
    my $range = abs($stop - $start);

    if($start <= $rightlim){ #startpos still < rightmost limit of stack, we are still in ongoing stack
      ## keep track of stop position and update rightlim if needed
      if($stop > $rightlim){ #start pos still in same stack, but stop pos implies to update rightmost limit
	$rightlim = $stop; #update rightmost limit of ongoing stack
	}

      ## store data
      $fastaII{$stacknr}{$cdsnr}{"protein"} = $fasta{$start}{$cds}{"protein"};
      $fastaII{$stacknr}{$cdsnr}{"start"} = $fasta{$start}{$cds}{"start"};
      $fastaII{$stacknr}{$cdsnr}{"stop"} = $fasta{$start}{$cds}{"stop"};
      $fastaII{$stacknr}{$cdsnr}{"target"} = $fasta{$start}{$cds}{"target"};
      $fastaII{$stacknr}{$cdsnr}{"cdsnr"} = $cdsnr;
      $fastaII{$stacknr}{$cdsnr}{"range"} = $fasta{$start}{$cds}{"range"};
      $fastaII{$stacknr}{$cdsnr}{"seq"} = $fasta{$start}{$cds}{"seq"};

      } else { #start pos > rightmost limit, we enter in a new stack
      ## open new stack
      $stacknr++; #update stack number
      $rightlim = $stop; #update rightmost limit of ongoing stack

      ## store data
      $fastaII{$stacknr}{$cdsnr}{"protein"} = $fasta{$start}{$cds}{"protein"};
      $fastaII{$stacknr}{$cdsnr}{"start"} = $fasta{$start}{$cds}{"start"};
      $fastaII{$stacknr}{$cdsnr}{"stop"} = $fasta{$start}{$cds}{"stop"};
      $fastaII{$stacknr}{$cdsnr}{"target"} = $fasta{$start}{$cds}{"target"};
      $fastaII{$stacknr}{$cdsnr}{"cdsnr"} = $cdsnr;
      $fastaII{$stacknr}{$cdsnr}{"range"} = $fasta{$start}{$cds}{"range"};
      $fastaII{$stacknr}{$cdsnr}{"seq"} = $fasta{$start}{$cds}{"seq"};
      }
    $cdsnr++;
#     print "$cdsnr\t$stacknr\t$start\t$stop\t$rightlim\t$protein\n";
    }
  }

### Organise CDS in stacks
my %prts;
foreach $stacknr (keys(%fastaII)){ #loop through all stacks
  $maxcount = 0; # limit to max 15 CDS per stack
  foreach $cdsnr (keys(%{$fastaII{$stacknr}})){ #loop through all cdsnr within protein name
    my $protein = $fastaII{$stacknr}{$cdsnr}{"protein"};
    my $start = $fastaII{$stacknr}{$cdsnr}{"start"};
    my $stop = $fastaII{$stacknr}{$cdsnr}{"stop"};
    $prts{$stacknr}{$cdsnr}{"protein"} = $cdsnr;
    $prts{$stacknr}{$cdsnr}{"start"} = $start;
    $prts{$stacknr}{$cdsnr}{"stop"} = $stop;
    $maxcount++;
    }
  }

### Find longest CDS of each stack
my %RES;
foreach $stacknr (keys(%prts)){ #loop over stacks
  my $maxrange = 0;
  my $best;
  foreach my $cdsnr (keys(%{$prts{$stacknr}})){ #check all CDS of stack
    my $s1 = $prts{$stacknr}{$cdsnr}{"start"};
    my $e1 = $prts{$stacknr}{$cdsnr}{"stop"};
    my $range = abs($e1 - $s1);
    if($range > $maxrange){
      $best = $cdsnr;
      $maxrange = $range;
      }
    }  #end of loop foreach cdsnr 
    $RES{$stacknr}{"keep"} = $best;
  } #end of loop foreach stacknr 


#### Keep only CDS that are OK
my %fastaIII;
foreach $stacknr (keys(%RES)){
  ## Check who was included in test
  $cdsnr = $RES{$stacknr}{"keep"};

  ## Store data
  $fastaIII{$stacknr}{$cdsnr}{"protein"} = $fastaII{$stacknr}{$cdsnr}{"protein"};
  $fastaIII{$stacknr}{$cdsnr}{"start"} = $fastaII{$stacknr}{$cdsnr}{"start"};
  $fastaIII{$stacknr}{$cdsnr}{"stop"} = $fastaII{$stacknr}{$cdsnr}{"stop"};
  $fastaIII{$stacknr}{$cdsnr}{"target"} = $fastaII{$stacknr}{$cdsnr}{"target"};
  $fastaIII{$stacknr}{$cdsnr}{"cdsnr"} = $cdsnr;
  $fastaIII{$stacknr}{$cdsnr}{"range"} = $fastaII{$stacknr}{$cdsnr}{"range"};
  $fastaIII{$stacknr}{$cdsnr}{"seq"} = $fastaII{$stacknr}{$cdsnr}{"seq"};  
  }

#### SPLIT CDS and save as Fasta
## Print fasta, suffix different cds of proteins with distinct codes
print "\nKeep best CDS of each stack:\nSTACKnr\tStart\tStop\tAnnotation\tSource\n";
if($add == 0){
  open(OUT, ">$outfile");
  open(ANN, ">$outfile.annot");
  } else {
  open(OUT, ">>$outfile");
  open(ANNOT, ">>$outfile.annot");
  }

$previous = 1e9;
foreach $stacknr (sort(keys(%fastaIII))){
  foreach $cdsnr (sort(keys(%{$fastaIII{$stacknr}}))){
    if($fastaIII{$stacknr}{$cdsnr}){
      my $protein = $fastaIII{$stacknr}{$cdsnr}{"protein"};
      if($protein){
	my $startpos = $fastaIII{$stacknr}{$cdsnr}{"start"};
	my $stoppos = $fastaIII{$stacknr}{$cdsnr}{"stop"};
	my $seq = $fastaIII{$stacknr}{$cdsnr}{"seq"};
	my $target = $fastaIII{$stacknr}{$cdsnr}{"target"};
	if($stacknr == $previous){
	  print OUT "$seq\n";
	  } else {
	  $previous = $stacknr;
	  print OUT ">$target\_PRT$stacknr\n$seq\n";
	  print ANNOT "$target\_PRT$stacknr\t$protein\n";
	  }
	print "$stacknr\t$startpos\t$stoppos\t$protein\t$target\n";
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

