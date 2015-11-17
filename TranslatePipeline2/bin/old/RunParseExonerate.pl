#!/bin/perl

#####################################
#### Performs Genewise analysis between db (protein) and sequences (dna)
#### and translated and framed sequences
####
#### Usage: perl RunParseExonerate.pl prot dna hitpairs output
####
#### Assumes presence of tmp/genewise
#### 
#####################################
use File::Basename;

my $prot = $ARGV[0]; #give path to protein database (fasta)
my $dna = $ARGV[1]; #give path to dna fasta
my $hitpairs = $ARGV[2]; #give path to blastx pairs (best hits among dna <-> prot)
my $out = $ARGV[3]; #give path of final outputs
my $bsn = basename($dna);
my $tmp = "tmp/genewise";
$scriptname = "RunParseExonerate";

chomp($prot);
chomp($dna);

## Parse blastx pairs and load it into hash
open(PAIRS, $hitpairs);

my %blast;
while(<PAIRS>){
  chomp();
  my @tmp = split(/\s+/, $_, 2);
  $idxdna = $tmp[0];
  chomp($idxdna);
  $idxprot = $tmp[1];
  chomp($idxprot);
  $blast{$idxdna}{$idxprot} = 1;
  }

## Iterate over input sequences (loop through %blast)
foreach $idxdna (keys(%blast)){ #loop over DNA queries (one DNA vs multiple proteins, loop DNA)
  open(TMP, ">$tmp/$bsn.dna.tmp");
  print TMP "$idxdna\n";
  close(TMP);
  
  # save all corresponding accessions from protein DB
  open(PRT, ">$tmp/$bsn.prot.tmp");
  foreach $idxprot (keys(%{$blast{$idxdna}})){ #loop
    print PRT "$idxprot\n";
    }
    close(PRT);

  # extract sequence from complete fasta, save it in tmp folder
  $command = "perl bin/FilterFasta.pl $dna $tmp/$bsn.dna.tmp tmp/genewise/$bsn.dnaquery";
  print "### $scriptname : $command\n";
  system("$command");

  # extract prot from complete ref, save it in tmp folder
  $command = "perl bin/FilterFasta.pl $prot $tmp/$bsn.prot.tmp tmp/genewise/$bsn.protref";
  print "### $scriptname : $command\n";
  system("$command");

  ## Run Exonerate, one sequence at a time
  $command = "./bin/exonerate --showalignment no --showvulgar no --ryo \">\%ti [\%tab - \%tae] \%qi all\n\%tas\n>\%ti [\%tab - \%tae] \%qi cds\n\%tcs\n\" --model protein2genome $tmp/$bsn.protref $tmp/$bsn.dnaquery > $tmp/$bsn.exonerate";
  print "### $scriptname : $command\n";
  system("$command");

  ## CDS files
  # Parse CDS and create extra files if multiple proteins in target sequence
  $command = "perl bin/exoneratecontig4.pl $tmp/$bsn.exonerate $tmp/$bsn.dirtycontig cds 0";
  print "### $scriptname : $command\n";
  system("$command");
  
  # Clean and translate CDS (includes split different proteins among fasta files) and translate if asked to
  $command = "perl bin/CleanTranslateCodon.pl $tmp/$bsn.dirtycontig $out/$bsn";
  print "### $scriptname : $command\n";
  system("$command");

  ## *.ALL files
  # Parse complete sequences (CDS + introns) and create extra files if multiple proteins in target sequence
  $command = "perl bin/exoneratecontig4.pl $tmp/$bsn.exonerate $out/$bsn.all.fas all 1";
  print "### $scriptname : $command\n";
  system("$command");
  } #end of looping over DNAs


## remove whitespaces
sub trim($){
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
  }