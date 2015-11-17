#!/bin/perl

# Usage:
# perl Translate_build.pl NCPU
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
$genome = $ARGV[1];
$db = $ARGV[2];
$bsn = basename($genome);

$scriptname = "Translate_build";

#### Pipeline
# run blast search, keep best hit of each queried accession
$command = "perl bin/RunParseBlastx.pl refs/$db data.in/$genome tmp/genewise best";
print "### $scriptname : $command\n";
system("$command");

# run genewise and parse out fasta files
$command = "perl bin/RunParseGenewise.pl refs/$db data.in/$bsn tmp/genewise/$genome.besthits.tmp tmp/fasta/$genome";
print "### $scriptname : $command\n";
system("$command");

# further clean fasta files
$command = "perl bin/CleanFasta.pl tmp/fasta/$genome.pep.fas data.out/pep";
print "### $scriptname : $command\n";
system("$command");

$command = "perl bin/CleanFasta.pl tmp/fasta/$genome.cdna.fas data.out/cdna";
print "### $scriptname : $command\n";
system("$command");

