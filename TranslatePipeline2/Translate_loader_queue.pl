#!/bin/perl

# Usage:
# perl Translate_loader_queue.pl dbase NCPU MaxUse MinRam
#
# dbase = plant (use this spelling), pick protein reference database (downloaded from GenBank, May 1st 2015).
# NCPU = Number of cores to be used by each process
# MaxUse = Maximal load allowed on the computing server.
# MinRam = Minimal RAM to leave free (Gb).
#
# Depends on diamond for doing blastx searches. this soft is provided in the bin/ folder.
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
# Nils Arrigo, Uni of Arizona 2012, Uni of Lausanne 2015
# use Cwd;
use threads;
use threads::shared;
use POSIX;
use File::Basename;

my $rootdir = getcwd;

#### Get script arguments
my $dbase = $ARGV[0];
my $CPU = $ARGV[1];
my $MaxUse = $ARGV[2];
my $MinRam = $ARGV[3]; #Gb
chomp($MinRam);
my $scriptname = "Translate_loader";

# start fresh
$command = "/bin/rm -rf tmp data.out";
print "### $scriptname : $command\n";
system("$command");

$command = "/bin/rm ./data.in/*.idx";
print "### $scriptname : $command\n";
system("$command");

# prepare projects folder
# $command = "mkdir -p data.out/ data.out/pep data.out/cdna tmp/ tmp/genewise tmp/fasta";
$command = "mkdir -p data.out/ tmp/ tmp/genewise tmp/fasta";
print "### $scriptname : $command\n";
system("$command");

# Check contents of refs
# @list = `ls refs/`;
# $db = @list[0];
# chomp($db);
my $db = "Genbank\_$dbase.fas";
  
  
# Prepare diamond and idx db for reference proteins (only if needed) 
unless(-e "refs/$db.dmnd"){
  $command = "./bin/diamond makedb --in refs/$db -d refs/$db";
  print "### $scriptname : $command\n";
  system("$command");
  }
  
unless(-e "refs/$db.idx"){
  $command = "./bin/fastaindex refs/$db refs/$db.idx";
  print "### $scriptname : $command\n";
  system("$command");
  }

  
#### Check contents of data.in and iterate over them
@list = `ls data.in/`;

my @threads;
foreach $genome (@list){
  # get CPU load
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + $CPU;

  # get RAM load
  my $ram = `free -g | sed -n '3s/ \\+/\t/gp' | cut -f4`;
  my $dbload =-s "refs/$db.dmnd";
  $dbload = 4 * $dbload / 1e9;
  my $ramleft = floor($ram - $dbload);
  
  print "### $scriptname : Current CPU = $load and free RAM = $ramleft Gb\n\n";
  
  if($load < $MaxUse & $ramleft >= $MinRam) {
    chomp($genome);

    #### Automated loader
    print "### $scriptname : RUNNING $genome\n";

    my $t = threads->new(\&rungenome, $genome, $CPU, $MaxUse, $db, $scriptname);
    push(@threads,$t);
    sleep(30);

    } else {
    print "### $scriptname : $genome is waiting... CPU or RAM load is maximised\n\n";
    sleep(30);
    redo;
    }
  }

foreach (@threads) {
  my $genome = $_->join;
  print "done with $genome\n";
  }


sub rungenome {
  ### Get CPU load
  my ($genome, $CPU, $MaxUse, $db, $scriptname) = ($_[0], $_[1], $_[2], $_[3], $_[4]);

  # Prepare index of DNA sequences
  $command = "./bin/fastaindex data.in/$genome tmp/fasta/$genome.idx";
  print "### $scriptname : $command\n";
  system("$command");

  # Start Translate pipepline
  $command = "perl Translate_pipeline_exonerate2.pl $CPU $genome $db";
  print "### $scriptname : $command\n";
  system("$command");
  }

# # clean mess
# $command = "rm *.log";
# print "### $scriptname : $command\n";
# system("$command");
# 
# $command = "rm -rf tmp";
# print "### $scriptname : $command\n";
# system("$command");
# 
# $command = "rm ./refs/*.fa.*";
# print "### $scriptname : $command\n";
# system("$command");
# 
# $command = "rm ./data.in/*.idx";
# print "### $scriptname : $command\n";
# system("$command");
