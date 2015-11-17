#!/bin/bash

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


############# Params #################
### Job params
ncpu=4          # number of CPUs to use per job
minram=5	# RAM to leave free (Gb)
maxload=46      # maximum server load allowed before starting a new job
nslices=46      # slice input data into nslices subjobs (parallelize)

### Exon cleaning
dbase=vertebrate_other   # GenBank protein database to use. Check docs in TranslatePipeline2/
              # only plant is available for the moment. dbases are from ftp://ftp.ncbi.nlm.nih.gov/refseq/release/
              # consider using invertebrate, vertebrate_mammalian and vertebrate_other (use this exact spelling)
              # download scripts are available in TranslatePipeline2/bin, run them from TranslatePipeline2/refs
              # warning: the larger the DB, the longer the analysis.

### Gene families
minlen=450      # length over which sequences must match with each other
minsim=30       # minimum percentage of similarity between matching sequences

### Mixture model params
ksmin=1e-9 #min ks (discard lower values)
ksmax=2 #max ks (discard larger values)
kmax=5 #max number of peaks being expected in Ks distribution, WARNING: analysis time increases with k
boots=5 #bootstrapping effort during search for optimal number of peaks, 
	        #WARNING: this is time consuming. Advised value is 1000
epsilon=1e-3; #convergence criterion; heuristics are stopped when loglik is improved by less than epsilon



######### Working script, do not modify ########################
echo "#####################################"
echo "Preparing folders"  
### Preparing outfolders / start fresh
/bin/rm -rf data.out
mkdir -p data.out
echo "Preparing folders done"
echo "#####################################"
echo 


echo "#####################################"
echo "Extracting exons"  
echo 
cd TranslatePipeline2/
/bin/rm data.in/*
perl bin/SliceFasta.pl ../data.in/*.fas data.in $nslices ../data.out/indexAccessions.txt
perl Translate_loader_queue.pl $dbase $ncpu $maxload $minram
/bin/rm -rf tmp/
ln -s ../TranslatePipeline2/data.out ../data.out/Exons
cd ..
echo "Extracting exons done"
echo "#####################################"
echo 


echo "#####################################"
echo "Producing gene families"  
echo 
cd BRHSingleLinkGeneFamilies/
perl BRHSingleLinkGeneFamilies_loader.pl $minlen $minsim $ncpu $maxload
/bin/rm -rf tmp/ *.log
/bin/rm *.log
ln -s ../BRHSingleLinkGeneFamilies/data.out/families ../data.out/GeneFamilies
cd ..
echo "Producing gene families done"
echo "#####################################"
echo 


echo "#####################################"
echo "Computing DnDs"  
echo 
cd KaKs/
./pipeline-macse-to-codeml.sh $maxload
ln -s ../KaKs/tmp/macse_alignments ../data.out/Alignments
ln -s ../KaKs/tmp/codeml_runs ../data.out/CodemlRuns
ln -s ../KaKs/data.out ../data.out/KsValues
cd ..
echo "Computing DnDs done"
echo "#####################################"
echo 


echo "#####################################"
echo "Detecting WGDs"  
echo 
cd Ks_mixmodels
perl RunMixModels.pl $ksmin $ksmax $kmax $boots $epsilon $maxload
ln -s ../Ks_mixmodels/data.out ../data.out/MixModels
cd ..
echo "Detecting WGDs done"
echo "#####################################"
echo 


echo "The analysis is finished"
echo 
