#!/bin/perl

# Usage:
# perl Translate_single.pl NCPU
#
#
# To be used to launch one file at a time.
#   For parallel version, use instead
#
#   Usage:
#   perl Translate_loader.pl NCPU
#
# expects that blastx (formatdb and blastall version) are already installed locally
#
# expects folders
#   data.in/ contains gene families (fasta) to frame and translate
#   refs/ contains protein database to be used as reference
#   bin/ contains perl toolbox
#
# Produces folder
#   data.out/ contains gene families (fasta) as cDNA (codon framed) and peptides
#
# WARNING: current version keeps introns
#
# Nils Arrigo, Uni of Arizona 2012

use Cwd;
use File::Basename;

my $rootdir = getcwd;

#### Get script arguments
$CPU = $ARGV[0];
$scriptname = "Translate_build";

# start fresh
$command = "rm -rf tmp data.out";
print "### $scriptname : $command\n";
system("$command");

# prepare projects folder
$command = "mkdir -p data.out/ data.out/pep data.out/cdna tmp/ tmp/genewise tmp/fasta";
print "### $scriptname : $command\n";
system("$command");


# Check contents of refs
@list = `ls refs/`;
$db = @list[0];
chomp($db);

# Prepare blast db
$command = "formatdb -i refs/$db -p T";
print "### $scriptname : $command\n";
system("$command");

#### Check contents of data.in and iterate over them
@list = `ls data.in/`;
foreach $genome (@list){
  chomp($genome);
  $bsn = basename($genome);
  
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
  
  if($load < $CPU) {
    print "### $scriptname : RUNNING $genome\n";
    
    # run blast search, keep best hit of each queried accession
    $command = "perl bin/RunParseBlastx.pl refs/$db data.in/$genome tmp/genewise";
    print "### $scriptname : $command\n";
    system("$command");

    # run genewise and parse out fasta files
    $command = "perl bin/RunParseGenewise.pl refs/$db data.in/$bsn tmp/genewise/$genome.besthits.tmp tmp/fasta/$genome";
    print "### $scriptname : $command\n";
    system("$command");

    # further clean fasta files
    $command = "perl bin/CleanFastaGenewise_bestlength.pl tmp/fasta/$genome.pep.fas data.out/pep";
    print "### $scriptname : $command\n";
    system("$command");

    $command = "perl bin/CleanFastaGenewise_bestlength.pl tmp/fasta/$genome.cdna.fas data.out/cdna";
    print "### $scriptname : $command\n";
    system("$command");

    } else {
    print "### $scriptname : WAIT... CPU load is maximised\n\n";
    sleep(60);
    redo;
    }
  }
    
# clean mess
$command = "rm *.log";
print "### $scriptname : $command\n";
system("$command");

$command = "rm -rf tmp";
print "### $scriptname : $command\n";
system("$command");

