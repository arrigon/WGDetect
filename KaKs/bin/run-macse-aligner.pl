#!/bin/perl

# Usage:
# perl run-macse-aligner.pl infolder outfolder MaxUse
#
# infolder = place where infiles are stored
# MaxUse = Maximal load allowed on the computing server.
#
# Nils Arrigo, Unil 2015
# use Cwd;

use threads;
use threads::shared;
use File::Basename;

my $rootdir = getcwd;

#### Get script arguments
my $infolder = $ARGV[0];
my $outfolder = $ARGV[1];
my $MaxUse = $ARGV[2];
my $scriptname = "macse_loader";


#### Check contents of infolder and iterate over them
@list = `ls $infolder/*.fas`;

my @threads;
foreach $aln (@list){
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
    chomp($aln);

    my $aln = basename($aln);
    
    #### Automated loader
    print "### $scriptname : RUNNING $aln\n";

    my $t = threads->new(\&runmacse, $aln, $infolder, $outfolder, $scriptname);
    push(@threads,$t);
    sleep(5);

    } else {
    print "### $scriptname : $aln is waiting... CPU load is maximised\n\n";
    sleep(5);
    redo;
    }
  }

foreach (@threads) {
  my $aln = $_->join;
  print "done with $aln\n";
  }


sub runmacse {
  ### params
  my ($aln, $infolder, $scriptname) = ($_[0], $_[1], $_[2], $_[3]);

  # Prepare index of DNA sequences    
  $command = "java -jar -Xmx500m bin/macse_v1.01.jar -prog alignSequences -seq $infolder/$aln -out_NT $outfolder/$aln.macse.aln";
  print "### $scriptname : $command\n";
  system("$command");
  }
