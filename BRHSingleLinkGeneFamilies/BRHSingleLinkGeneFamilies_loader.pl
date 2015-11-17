#!/bin/perl

# Produces gene families from transcriptomes
# Usage:
# perl BRH_SingleLinkFamilies_build.pl minlen minsim NCPU MaxUse
#
# minlen = length over which sequences must match with each other
# minsim = minimum percentage of similarity between matching sequences
# NCPU = Number of cores to be used by each process
# MaxUse = Maximal load allowed on the computing server.
#
# expects that blastx (formatdb and blastall version) are already installed locally
#
# expects folders
#   data.in/ contains transcriptomes to cluster into gene families
#   config/ contains blastline.conf, leave multiple core parameter undefined, will be added by scripts
#   bin/ contains perl toolbox
#
# Produces folder
#   data.out/ contains gene families (fasta) as cDNA (codon framed) and peptides
#
# Nils Arrigo, Uni of Arizona 2012

# hard-coded params: will be passed to from command line

## Make sure we start fresh
$command = "rm -rf tmp/ data.out/* *.log";
print "### $scriptname : $command\n";
system("$command");

use Cwd;
use File::Basename;

my $rootdir = getcwd;

### Get script arguments
$minlen = $ARGV[0];
$minsim = $ARGV[1];
$NCPU = $ARGV[2];
$MaxUse = $ARGV[3];
$scriptname = "BRHSingleLinkGeneFamilies_build.pl";

# Start fresh
$command = "rm -rf data.out/families/* tmp/ *.log";
print "### $scriptname : $command\n";
system("$command");

# prepare projects folder
$command = "mkdir -p data.out/ data.out/families tmp/";
print "### $scriptname : $command\n";
system("$command");

#### Check contents of data.in and iterate over them
@list = `ls data.in/*.cds.fas`;
foreach $genome (@list){
  chomp($genome);
  $bsn = basename($genome);
  $bsn =~ s/\..*$//;

  ### Get CPU load
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + 1;

  print "### $scriptname : Current CPU load is $load\n\n";
  
  if($load < $MaxUse) { 
    # Label properly input fasta
    $command = "perl bin/PrefixFasta.pl $genome $bsn tmp/";
    print "### $scriptname : $command\n";
    system("$command");

    } else {
    print "### $scriptname : WAIT... CPU load is maximised\n\n";
    sleep(60);
    redo;
    }
  }

# Merge all these guys into a single fasta file
$command = "cat tmp/* > tmp/ALLDB.fas";
print "### $scriptname : $command\n";
system("$command"); 

# Perform blast search
# read blasline.conf
open BLA, "config/blastline.conf";
my @tmp = <BLA>;
my @config = ();
foreach my $val (@tmp){
  chomp($val);
  if($val !~ /^$/ & $val !~ /^\#/){
    push(@config, $val);
    }
  }
 
# makeblast db 
$command = "$config[0]";
print "### $scriptname : $command\n";
system("$command");

# run blast search
if($load < $MaxUse) { 
  $command = "$config[1] -num\_threads $NCPU";
  print "### $scriptname : $command\n";
  system("$command");
  } else {
  print "### $scriptname : WAIT... CPU load is maximised\n\n";
  sleep(60);
  redo;
  } 

# Parse blast output
$command = "perl bin/KeepBestBlastHits.pl tmp/AllvsAll.blast tmp/AllvsAll.best $minlen $minsim";
print "### $scriptname : $command\n";
system("$command");

# Produce gene families
$command = "perl bin/SLClusGeneFam.pl tmp/AllvsAll.best tmp/ALLDB.fas data.out";
print "### $scriptname : $command\n";
system("$command");

# Clean mess
# $command = "/bin/rm -rf tmp/";
# print "### $scriptname : $command\n";
# system("$command");

