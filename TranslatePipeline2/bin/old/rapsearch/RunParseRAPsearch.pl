#!/bin/perl

#####################################
#### Performs blastx between db (protein) and sequences (dna)
#### and extracts best protein hit for each dna accession
####
#### Produces ../tmp/besthits.tmp storing this information
####
#### WARNING: returns only best hit
#### Nils Arrigo, Unil 2015
#####################################
use File::Basename;

my $db = $ARGV[0];
my $target = $ARGV[1];
my $outfolder = $ARGV[2];
my $mode = $ARGV[3];
my $NCPU = $ARGV[4];
my $bsn = basename($target);
my $out = "$outfolder/$bsn.besthits.tmp";
$scriptname = "RunParseRAPsearch";

print "$NCPU";

chomp($mode);
chomp($ID);

## Run rapsearch
# $command = "blastall -p blastx -m 8 -d $db -i $target -e 1e-3 -P 90 -K 1 > $out";
# $command = "blastall -p blastx -m 8 -d $db -i $target -e 1e-3 -P 90 -K 1 | cut -f2,11 -s > $out";
# $command = "blastall -p blastx -m 8 -a $NCPU -d $db -i $target -e 1e-5 -P 50 | cut -f1,2,11 -s > $out.all"; #Version 2.2.26
$command = "./bin/rapsearch -q $target -d $db.rapdb -z $NCPU -e 1e-5 -u 1 | cut -f1,2,11 -s > $out.tmp";
print "### $scriptname : $command\n";
system("$command");

$command = "sed \'1d\' $out.tmp > $out.all";
print "### $scriptname : $command\n";
system("$command");


## Parse output and keep best candidate
if($mode eq "best"){
  print "### $scriptname : Parse rapsearch output, keep best hit of each accession\n";
  open(HD, "$out\.all") or die "No $out.all";
  my %best;
  while (<HD>){
    my @tmp = split(/\s/, $_, 3);
    my $query = $tmp[0];
    my $ref = $tmp[1];
    my $val = $tmp[2];
    if($best{$query}){
      $start = $best{$query}{"val"};
      if($val <= $start){
	$start = $val;
	$best{$query}{"ref"} = $ref;
	$best{$query}{"val"} = $val;
	}  
      } else {
      my $start = 10;
      if($val <= $start){
	$start = $val;
	$best{$query}{"ref"} = $ref;
	$best{$query}{"val"} = $val;
	}  
      }
    }
  close HD;
  ## Print results into final file
  open(OUT, ">$out");
  foreach $query (keys(%best)){
    my $ref = $best{$query}{"ref"};
    print OUT "$query\t$ref\n";
    }
  close(OUT);

} elsif($mode == "all") {
  print "### $scriptname : Parse rapsearch output, keep all hits\n";
  open(HD, "$out\.all") or die "No $out.all";
  my %best;
  while (<HD>){
    my @tmp = split(/\s/, $_, 3);
    my $query = $tmp[0];
    my $ref = $tmp[1];
    my $val = $tmp[2];
    $best{$query}{$ref} = $val;
    }
  close HD;

  ## Print results into final file
  open(OUT, ">$out");
  foreach $query (keys(%best)){
    foreach $ref (keys(%{$best{$query}})){
      print OUT "$query\t$ref\n";
      }
    }
  close(OUT);
  }

