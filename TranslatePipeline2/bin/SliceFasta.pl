#!/bin/perl

#####################################
#### Slice Fasta file into n subfiles
####
#### and produce an index of accessions in files.
####
#### perl SliceFasta.pl infile outfolder nchunks
#### Nils Arrigo, Unil 2015
#####################################
use File::Basename;
use POSIX;

my $infile = $ARGV[0];
my $outfolder = $ARGV[1];
my $splits = $ARGV[2];
my $indexfile = $ARGV[3];
chomp($indexfile);

my $bsn = basename($infile);

# Open fasta
open(FILE, "$infile");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

my %fasta;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  $fasta{$tmp[0]} = $seq; 
  }

# Save it into $splits chuncks
my @accs = keys(%fasta);
my $nseq = @accs;
my $nseqchunck = ceil($nseq / $splits);

print "Slicing $infile into $splits chunks of about $nseqchunck sequences\n";

my $cnt = 0;
my $seqcnt = 0;
open(LIST, ">$indexfile");
print(LIST "seqID\toriginalSeqHeader\tanalysisBatch\n");
while ( my @chunks = splice(@accs, 0, $nseqchunck) ) {
  open(OUT, ">$outfolder/$cnt.$bsn");
  foreach $acc (@chunks){
    $seq = $fasta{$acc};
    print(OUT ">$seqcnt\n$seq\n");
    print(LIST ">$seqcnt\t$acc\t$cnt\n");
    $seqcnt++;
    }  
  $cnt++;
  }

close(OUT);
close(LIST);