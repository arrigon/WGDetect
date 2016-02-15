#!/bin/bash
#
# pipeline to 
# 1) run macse, 
# 2) transform fasta to paml
# 3) run fasttree 
# 4) prepare codeml input
# 5) run codeml

# dependencies
#
# programs and scripts:
#   run-macse-aligner.sh
#      macse_v0.9b1.jar 
#   fasta2paml-converter.sh
#      fasta2paml.SZ.pl
#   run-FastTree.sh
#      FastTree
#   codeml-ctl-generator.sh
#   run-codeml.sh
#     codeml
# files: 
#   codeml.ctl.body.template  

maxload=$1	
tasks=$2 	# Flow control; run pipeline for a given set of tasks 
		# (to use in case of crash, or for running partial analyses).
		#1 = start from alignments, 
		#2 = skip alignments and perform last steps (faster)
		
		# Examples: 
		# - tasks = 12 runs the complete pipeline
		# - tasks = 2 skips the alignment step


echo "#####################################"
echo "make sure we start fresh"  
echo 
/bin/rm -rf data.out tmp
echo "#####################################"
echo 


if [[ $tasks =~ 1 ]]; then
  echo "#####################################"
  echo "running macse alignments"  
  echo 
  mkdir tmp tmp/macse_alignments
  ln data.in/*.fas tmp/macse_alignments
  # ./bin/run-macse-aligner.sh
  perl ./bin/run-macse-aligner.pl data.in tmp/macse_alignments $maxload
  /bin/rm -f data.in/Fam.*_macse_AA.fasta
  echo "macse alignments done"
  echo "#####################################"
  echo 
else 
echo "#####################################"
  echo "Skipping alignments"  
  echo 
  mkdir tmp tmp/macse_alignments
  cp -l ../data.in/Fam.*aln tmp/macse_alignments
  echo "#####################################"
  echo 
fi

  
  

echo "#####################################"
echo "converting fasta to paml format"
echo 
mkdir tmp/paml_format
cd tmp/paml_format
ln -s ../macse_alignments/Fam.*aln .
# make paml file  
../../bin/fasta2paml-converter.sh
cd ../..
echo "converting done"
echo "#####################################"
echo 


echo "#####################################"
echo "running FastTree "
echo 
mkdir tmp/fasttree_runs
cd tmp/fasttree_runs
ln ../paml_format/Fam*paml .
../../bin/run-FastTree.sh
for file in `ls -1 Fam*paml` ; do unlink $file ; done
cd ../..
echo "FastTree runs done"
echo "#####################################"
echo 


echo "#####################################"
echo "preparing codeml input files"
echo 
mkdir tmp/codeml_runs
cd tmp/codeml_runs
ln -s ../paml_format/Fam*paml  .
ln -s ../fasttree_runs/Fam*tree  .
ln -s ../../bin/codeml.ctl.body.template .
../../bin/codeml-ctl-generator.sh
unlink codeml.ctl.body.template 
cd ../..
echo "done preparing"
echo "#####################################"
echo 


echo "#####################################"
echo "running codeml "
echo 
cd  tmp/codeml_runs
../../bin/run-codeml.sh
for file in `ls -1 Fam*paml ` ; do unlink $file ; done
for file in `ls -1 Fam*tree ` ; do unlink $file ; done
cd ../..
echo "done codeml"
echo "#####################################"
echo 



echo "#####################################"
echo "collecting Ks values"
echo 
mkdir data.out
R CMD BATCH '--args infolder="tmp/codeml_runs" outfile="data.out/all.dS.values.txt"' bin/ExtractKsValues.r
echo "done extracting Ks values"
echo "#####################################"
echo 
