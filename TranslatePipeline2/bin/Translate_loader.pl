#!/bin/perl

# Usage:
# perl Translate_loader.pl NCPU MaxUse
#
# NCPU = Number of cores to be used by each process
# MaxUse = Maximal load allowed on the computing server.
#
# expects that blastx (formatdb and blastall version) are already installed locally
#
# expects folders
#   data.in/ contains gene families (fasta) to frame and translate
#   refs/ contains protein database to be used as reference
#   bin/ contains perl toolbox
#
# Results are in folder data.out/
#   *.cds.fas = framed DNA sequence (CDS only)
#   *.pep.fas = translated protein
#   *.all.fas = DNA sequence (CDS + introns), WARNING: not in frame over complete sequence.
#   *.all.fas.annot = correspondances between DNA sequence and accession from protein database
#
# NOTE: It is more productive to slice big input files into smaller chuncks 
#	and let the pipeline run in paralell on the wished number of cores
# You can slice your dataset using bin/SliceFasta.pl inputfile outfolder nslices
# Typically, use the following command:
# perl bin/SliceFasta.pl MyBigSingleInput.fasta data.in NCPU
#
# Nils Arrigo, Uni of Arizona 2012

use Cwd;
use File::Basename;

my $rootdir = getcwd;

#### Get script arguments
$CPU = $ARGV[0];
$MaxUse = $ARGV[1];
$scriptname = "Translate_loader";

# start fresh
$command = "rm -rf tmp data.out";
print "### $scriptname : $command\n";
system("$command");

$command = "rm ./data.in/*.idx";
print "### $scriptname : $command\n";
system("$command");

# prepare projects folder
# $command = "mkdir -p data.out/ data.out/pep data.out/cdna tmp/ tmp/genewise tmp/fasta";
$command = "mkdir -p data.out/ tmp/ tmp/genewise tmp/fasta";
print "### $scriptname : $command\n";
system("$command");

# Check contents of refs
@list = `ls refs/`;
$db = @list[0];
chomp($db);

# Prepare blast db for reference proteins
$command = "formatdb -i refs/$db -p T";
print "### $scriptname : $command\n";
system("$command");

# Prepare index of reference proteins
$command = "./bin/fastaindex refs/$db refs/$db.idx";
print "### $scriptname : $command\n";
system("$command");

#### Check contents of data.in and iterate over them
@list = `ls data.in/`;
foreach $genome (@list){
  chomp($genome);
  
  ### Get CPU load
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + $CPU;

  print "### $scriptname : Current CPU load is $load\n\n";
  
  if($load < $MaxUse) {
    #### Automated loader
    print "### $scriptname : RUNNING $genome\n";

    # Prepare index of DNA sequences
    $command = "./bin/fastaindex data.in/$genome tmp/fasta/$genome.idx";
    print "### $scriptname : $command\n";
    system("$command");

    # Start Translate pipepline
    $command = "perl Translate_pipeline_exonerate.pl $CPU $genome $db &";
    print "### $scriptname : $command\n";
    system("$command");
    sleep(10);

    } else {
    print "### $scriptname : WAIT... CPU load is maximised\n\n";
    sleep(10);
    redo;
    }
  }
    
# clean mess
# $command = "rm *.log";
# print "### $scriptname : $command\n";
# system("$command");

# $command = "rm -rf tmp";
# print "### $scriptname : $command\n";
# system("$command");

# $command = "rm ./refs/*.fa.*";
# print "### $scriptname : $command\n";
# system("$command");

# $command = "rm ./data.in/*.idx";
# print "### $scriptname : $command\n";
# system("$command");
