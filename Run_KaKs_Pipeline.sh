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
### Job params and pipeline control
ncpu=4          # number of CPUs to use per job
minram=5	# RAM to leave free (Gb)
maxload=24      # maximum server load allowed before starting a new job
nslices=46      # slice input data into nslices subjobs (parallelize)

tasks=1234	# Flow control; run pipeline for a given set of tasks 
		# (to use in case of crash, or for running partial analyses).
		# 1 = Exons cleaning (input = fas -> output = *.cds.fas, *.pep.fas, *.all.fas, *.annot)
		# 2 = Clustering (input = *.cds.fas -> output = data.out/families/*.fas)
		# 3 = KaKs computing (input = *.fas -> output = aln, trees, phyml and Ks values)
		# 4 = Detection of WGDs in Ks distributions (input = all.dS.values.txt -> output = mixture models, pdf)

		# Examples: 
		# - tasks = 1234 runs the complete pipeline
		#
		# If running some subparts of the pipeline, 
		# make sure that relevant inputs (formatted as expected by pipeline) are available from data.in/
		# - tasks=234  resumes pipeline from step 2, takes *.cds.fas as input in data.in/
		# - tasks=34   resumes pipeline from step 3, takes gene families (*.fas) as input in data.in/
		# - tasks=4    resumes pipeline from step 4, takes ks values (all.dS.values.txt) as input in data.in/
		# - tasks=1	 runs only step 1, the same goes for running other unique steps

### Exon cleaning
dbase=vertebrate_other   # GenBank protein database to use. Check docs in TranslatePipeline2/
              # only plant is available for the moment. dbases are from ftp://ftp.ncbi.nlm.nih.gov/refseq/release/
              # consider using invertebrate, vertebrate_mammalian and vertebrate_other (use this exact spelling)
              # download scripts are available in TranslatePipeline2/bin, run them from TranslatePipeline2/refs
              # warning: the larger the DB, the longer the analysis.

              
### Gene families
minlen=300      # length over which sequences must match with each other
minsim=40       # minimum percentage of similarity between matching sequences

### Mixture model params
ksmin=1e-9 #min ks (discard lower values)
ksmax=2 #max ks (discard larger values)
kmax=4 #max number of peaks being expected in Ks distribution, WARNING: analysis time increases with k
boots=1 #bootstrapping effort during search for optimal number of peaks, 
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



################################################################
###### STEP 1 - Exons cleaning #################################
################################################################

if [[ $tasks =~ 1 ]]; then
  echo "#####################################"
  echo "Exon cleaning"  
  echo 
  cd TranslatePipeline2/
  /bin/rm data.in/*
  perl bin/SliceFasta.pl ../data.in/*.fas data.in $nslices ../data.out/indexAccessions.txt
  perl Translate_loader_queue.pl $dbase $ncpu $maxload $minram
  /bin/rm -rf tmp/
  ln -s ../TranslatePipeline2/data.out ../data.out/Exons
  cd ..
  echo "Exon cleaning done"
  echo "#####################################"
  echo 
else 
echo "#####################################"
  echo "Skipping exon cleaning"  
  echo 
  cd TranslatePipeline2/
  if [[ $tasks =~ 2 ]]; then
    /bin/rm data.out/*
    
    # check how many cds files we have in data.in
    nfls=$(find ../data.in/*.cds.fas -type f | wc -l)
    
    if [[ $nfls == 1 ]]; then #one file = typical input with CDS sequences from external source (one fasta containing all CDS)
      perl bin/SliceFasta.pl ../data.in/*.cds.fas data.out $nslices ../data.out/indexAccessions.txt
      else #several files = assuming we deal with an ongoing run of the WGDetect pipeline (and that cds have already been numbered)
      cp -l ../data.in/*.cds.fas data.out/
      fi
      
    fi
    
  cd ..
  echo "#####################################"
  echo 
fi



################################################################
###### STEP 2 - Gene families ##################################
################################################################

if [[ $tasks =~ 2 ]]; then
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
else 
  echo "#####################################"
  echo "Skipping gene families"  
  echo 
  cd BRHSingleLinkGeneFamilies/
  
  if [[ $tasks =~ 3 ]]; then
    /bin/rm -rf data.out/families
    mkdir data.out/ data.out/families
    cp -l ../data.in/*.aln data.out/families
    cp -l ../data.in/*.fas data.out/families
    cp -lr ../data.in/families data.out/
    fi
    
  cd ..
  echo "#####################################"
  echo 
fi



################################################################
###### STEP 3 - Ks values ######################################
################################################################

if [[ $tasks =~ 3 ]]; then
  echo "#####################################"
  echo "Computing DnDs"  
  echo 
  cd KaKs/
  
  # check if input files are normal or aligned fasta
  nfls=$(find data.in/*.aln -type f | wc -l)    
  if [[ $nfls == 0 ]]; then #no aln files, we are thus dealing with normal fasta files
    KAKSsteps=12
  else
    KAKSsteps=2
  fi
    
  ./pipeline-macse-to-codeml.sh $maxload $KAKSsteps
  ln -s ../KaKs/tmp/macse_alignments ../data.out/Alignments
  ln -s ../KaKs/tmp/codeml_runs ../data.out/CodemlRuns
  ln -s ../KaKs/data.out ../data.out/KsValues
  cd ..
  echo "Computing DnDs done"
  echo "#####################################"
  echo 
else
  echo "#####################################"
  echo "Skipping DnDs"  
  echo 
  cd KaKs/
  /bin/rm -rf data.out/
  mkdir data.out/

  if [[ $tasks =~ 4 ]]; then
    cp -l ../data.in/all.dS.values.txt data.out/all.dS.values.txt
    fi
    
  cd ..
  echo "#####################################"
  echo 
fi
  
  
  
################################################################
###### STEP 4 - Mixture models##################################
################################################################
  
if [[ $tasks =~ 4 ]]; then
  echo "#####################################"
  echo "Mixture models for detecting WGDs"  
  echo 
  cd Ks_mixmodels
  perl RunMixModels.pl $ksmin $ksmax $kmax $boots $epsilon $maxload
  ln -s ../Ks_mixmodels/data.out ../data.out/MixModels
  cd ..
  echo "Detecting WGDs done"
  echo "#####################################"
  echo 
else 
  echo "#####################################"
  echo "Skipping mixture models WGDs"  
  echo 
  echo "#####################################"
  echo 
fi



echo "The analysis is finished"
echo 
