#!/bin/perl

#####################################
#### Performs Genewise analysis between db (protein) and sequences (dna)
#### and translated and framed sequences
####
#### perl prot dna output
#### 
#####################################
use File::Basename;

my $prot = $ARGV[0]; #give path to protein database (fasta)
my $dna = $ARGV[1]; #give path to dna fasta
my $hitpairs = $ARGV[2]; #give path to blastx pairs (best hits among dna <-> prot)
my $out = $ARGV[3]; #give path of final outputs
my $bsn = basename($dna);
my $tmp = "tmp/genewise";
$scriptname = "RunParseGenewise";

chomp($prot);
chomp($dna);

# open blastx pairs
open(PAIRS, $hitpairs);

# # prepare outfiles
# open(PROT, ">$out.cdna.fas");
# open(DNA, ">$out.pep.fas");

## Iterate over input sequences
while(<PAIRS>){
  chomp();
  my @tmp = split(/\s+/, $_, 2);
  $idxdna = $tmp[0];
  chomp($idxdna);
  $idxprot = $tmp[1];
  chomp($idxprot);

  open(TMP, ">$tmp/$bsn.dna.tmp");
  print TMP "$idxdna\n";
  close(TMP);

  open(TMP, ">$tmp/$bsn.prot.tmp");
  print TMP "$idxprot\n";
  close(TMP);

  # extract sequence from complete fasta, save it in tmp folder
  $command = "perl bin/FilterFasta.pl $dna $tmp/$bsn.dna.tmp tmp/genewise/$bsn.dnaquery";
  print "### $scriptname : $command\n";
  system("$command");

  # extract prot from complete ref, save it in tmp folder
  $command = "perl bin/FilterFasta.pl $prot $tmp/$bsn.prot.tmp tmp/genewise/$bsn.protref";
  print "### $scriptname : $command\n";
  system("$command");

  # Run Genewise, one sequence at a time
  $command = "genewisedb $tmp/$bsn.protref $tmp/$bsn.dnaquery -both -cdna -trans -pep -pretty -gener -alg 333 -silent > $tmp/$bsn.genewise";
  print "### $scriptname : $command\n";
  system("$command");

#   # print fasta headers
#   print(PROT "\n>$idxdna\n");
#   print(DNA "\n>$idxdna\n");
  
  # Do Mike's magic to parse Genewise outputs and produce pep and dna sequences; KEEP INTRONS and sequences POST STOP CODON
  $command = "perl bin/estcontig.pl < $tmp/$bsn.genewise > $tmp/$bsn.dirtycontig";
  print "### $scriptname : $command\n";
  system("$command");

  $command = "cat $tmp/$bsn.dirtycontig | perl bin/xout.pl dna $idxdna $out.cdna.fas";
  print "### $scriptname : $command\n";
  system("$command");

  $command = "cat $tmp/$bsn.dirtycontig | perl bin/xout.pl prot $idxdna $out.pep.fas";
  print "### $scriptname : $command\n";
  system("$command");
 
  # Do Mike's magic to parse Genewise outputs and produce pep and dna sequences REMOVE EVERYTHING AFTER FIRST STOP CODON (uncomment previous three commands if wish using this version)
#   $command = "perl bin/stopout.pl < $tmp/$bsn.dirtycontig > $tmp/$bsn.nostopcontig";
#   print "### $scriptname : $command\n";
#   system("$command");

#   $command = "cat $tmp/$bsn.nostopcontig | perl bin/xout.pl dna $idxdna $out.cdna.fas";
#   print "### $scriptname : $command\n";
#   system("$command");
# 
#   $command = "cat $tmp/$bsn.nostopcontig | perl bin/xout.pl prot $idxdna $out.pep.fas";
#   print "### $scriptname : $command\n";
#   system("$command");

  } #end of looping over accessions

#   close(PROT);
#   close(DNA);


## remove whitespaces
sub trim($){
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
  }