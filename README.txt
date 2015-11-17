#####################################
###### WGD detection pipeline #######
#####################################
### usage: ./Run_Kaks_Pipeline.sh
###
### INPUT: One fasta file containing your DNA contigs (can be anything)
###        the only limitation might be with sequence names.
###        there might be some conflicts with codeml if your contig names are >15-20 characters
###
### Params : check out hard-coded params at lines >45 of this script.
### OUTPUT: all important outputs are stored in data.out/
###
### This pipeline runs the following analysis steps
### - Extracting exons from input sequences (provides annotations, check README in TranslatePipeline2/).
### - Clusters exons into gene families, using reciprocal best blast hits and single linkage clustering, check README in BRHSingleLinkGeneFamilies/
### - Aligns the obtained clusters (using macse) and compute DN/DS stats (using codeml). Check KaKs pipeline
### - run mixture models to detect peaks in Ks stats
###
### Dependencies:
## OS: linux (BioLinux would embark most of the needed softs).
## perl modules (these are standard modules, use CPAN to install them).
# Cwd 
# File::Basename
# threads
# threads::shared
# Bio::SeqIO
# Bio::Seq
# List::Util
# POSIX
# IO::File
# Getopt::Std
#
## R packages
# mixtools
# ape
#
## Stand-alone programs (most are already provided in bin/ folders of the sub-pipelines
# blast suite: makeblastdb and blastn 2.2.28+ (must be installed in your server)
#
#####################################
### Version 1.03
### Nils Arrigo (Unil, Uni. of Arizona), Stephan Zoller (ETHZ) and Celine Geiser (UniNE), 2013-2015
### some code sections from Erik R Hanschen, Katrina Dlugosch and Michael S. Barker, Uni. of Arizona
#####################################