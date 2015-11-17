#!/bin/perl

# Usage:
# perl Translate_build.pl NCPU
#
# expects that blastx (formatdb and blastall version) and exonerate are already installed locally
#
# expects folders
#   data.in/ contains gene families (fasta) to frame and translate
#   refs/ contains protein database to be used as reference
#   bin/ contains perl toolbox
#
# Produces folder
#   data.out/ with 
#	*.all.fas = CDS + INTRONS
#	*.cds.fas = CDS only
#	*.pep.fas = translated peptides
#
# Nils Arrigo, Uni of Arizona 2012

use Cwd;
use File::Basename;

my $rootdir = getcwd;

#### Get script arguments
$NCPU = $ARGV[0];
$genome = $ARGV[1];
$db = $ARGV[2];
$bsn = basename($genome);

$scriptname = "Translate_pipeline_exonerate";

#### Pipeline
# run blast search, keep best hit of each queried accession
$command = "perl bin/RunParseDiamond.pl refs/$db data.in/$genome tmp/genewise all $NCPU";
print "### $scriptname : $command\n";
system("$command");

# run exonerate and parse out fasta files
$command = "perl bin/RunParseExonerate.pl refs/$db data.in/$bsn tmp/genewise/$genome.besthits.tmp data.out/";
print "### $scriptname : $command\n";
system("$command");