#!/bin/perl

# Removes any sequences coming after stop codon and translates remainder (input = DNA)
# Saves pruned DNA and amino acid sequences to *.cds.fas and *.prot.fas
#
# WORKS ONLY WITH CDS. DO NOT USE IF SEQUENCE CONTAINS INTRONS
#
# Usage: perl CleanTranslateCodon.pl infile outfile
#
# Nils Arrigo, Unil 2015
###
use Bio::SeqIO;
use Bio::Seq; 

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

## import fasta file
open(FILE, "$infile");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

## clean from stop codons
open(DNA, ">>$outfile.cds.fas");
open(PRT, ">>$outfile.pep.fas");
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  my $acc = trim($tmp[0]);
  my $seq = $tmp[1];
  $seq =~ s/[\n*|\s*]//g;

  # translate
  my $seqobjnuc;
  $seqobjnuc = Bio::PrimarySeq->new(-seq=>$seq);
  my $seqobjprot = $seqobjnuc->translate;

  # check whether contains stop codon
  if ($seqobjprot->seq =~ /\D*/){ #yes: stop codon found in protein sequence
    my (@triplet) = ();

    #split sequence in codons (triplet), and visit it until finds stop codon
    my $tmpseq = $seqobjnuc->seq;
    while($tmpseq =~ /(\w{3})/g){
      my $codon = $1;
      next if $codon =~ /TAG/i || $codon =~ /TGA/i || $codon =~ /TAA/i; #end while loop if finds stop codon, /i stands for font-independant
      push @triplet, $codon;
      }

    #update Bio::Seq object
    $seqobjnuc->seq(join('', @triplet));
    $seqobjprot->seq($seqobjnuc->translate->seq);
    }

  # save cds
  my $seqout = $seqobjnuc->seq();
  print DNA ">$acc\n$seqout\n";

  # save nucleic acids
  $seqobjprot = $seqobjnuc->translate;
  my $seqout = $seqobjprot->seq;
  print PRT ">$acc\n$seqout\n";
  }
close(DNA);
close(PRT);

## remove whitespaces
sub trim($){
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
  }