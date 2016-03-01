#!/bin/perl

##### Script to run mix models on several datasets; in parallel mode.
#
# Usage: perl RunMixModels.pl MaxCPU
# - MaxCPU = maximum CPU load allowed before launching the next job
# WARNING: Check out hard-coded params at lines 28-43. Adapt them as needed.
#
# Requirements:
# - perl module "File::Basename" must be installed locally
# - R package "mixtools" must be installed locally
# - inputs must be in data.in (*.scafSeq)
# - outputs will be generated in data.out
#
# the R script doing the job is bin/MixModels.r
#
# Outputs are in data.out/ 
# *.pdf = best fits
# *.mixmodels.txt = params of detected peaks 
#    (centers = mean, 
#    stdev = standard deviation, 
#    contribs = proportion of observations contributed by each peak, 
#    loglik = loglik of run)
# *.Rout = analysis logs, keep them to know what parameter values were used.
#
# Nils Arrigo, Unil 2015.
########################################################
use File::Basename;


#### Get script arguments
my $ksmin = $ARGV[0];
my $ksmax = $ARGV[1];
my $kmax = $ARGV[2];
my $boots = $ARGV[3];
my $epsilon = $ARGV[4];
my $CPU = $ARGV[5];


# start fresh
$command = "rm -rf data.out/*";
print "### $scriptname : $command\n";
system("$command");


## WARNING: Hard coded params
# I/O
my $infolder = "data.in"; #leave as is
my $outfolder = "data.out"; #leave as is

# # KS window
# my $ksmin = 1e-9; #min ks
# my $ksmax = 2; #max ks
# 
# # Mixture model params
# my $kmax = 5; #max number of peaks being expected, WARNING: analysis time increases with k
# my $boots = 5; #bootstrapping effort during search for optimal number of peaks, 
# 	        #WARNING: this is time consuming. Advised value is 1000
# my $epsilon = 1e-3; #convergence criterion; heuristics are stopped when loglik is improved by less than epsilon

# Graph params
my $breaks = 50; #number of breaks on histogram



##### Launch jobs on parallel cores
### List all files to analyse
my @list = `ls $infolder/*.txt`;
chomp(@list);

### loop over all datasets
foreach my $file (@list){
  chomp($file);
  my $bsn = basename($file);
  
  ### Get CPU load
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + 1;
  print "Current CPU load is $load \n";
  if($load < $CPU) {
    print "OK --- start $file\n";
    $command = "R CMD BATCH '--args infile=\"$file\" outfolder=\"$outfolder\" ksmin=$ksmin ksmax=$ksmax kmax=$kmax boots=$boots epsilon=$epsilon breaks=$breaks' bin/MixModelsKs.r $outfolder/$bsn.Rout";
    print "\n##########\n##########\nRUNNING:\n$command\n##########\n\n\n";
    system("$command");
    sleep(10);
    } else {
    print "WAIT... CPU load is maximised\n";
    sleep(10);
    redo;
    }
  }
