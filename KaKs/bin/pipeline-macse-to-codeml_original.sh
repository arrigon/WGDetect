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


echo "#####################################"
echo "running macse alignments"  
echo 
mkdir macse_alignments
cd macse_alignments
ln -s ../../../BRHSingleLinkGeneFamilies/data.out/families/Fam*.fas .
run-macse-aligner.sh 
for file in `ls -1 Fam*.fas` ; do unlink $file ; done
cd ..
echo "macse alignments done"
echo "#####################################"
echo 

echo "#####################################"
echo "converting fasta to paml format"
echo 
mkdir paml_format
cd paml_format
ln -s ../macse_alignments/Fam.*_DNA.fasta .
# make paml file  
fasta2paml-converter.sh
for file in `ls -1 Fam.*_DNA.fasta` ; do unlink $file ; done
cd ..
echo "converting done"
echo "#####################################"
echo 

echo "#####################################"
echo "running FastTree "
echo 
mkdir fasttree_runs
cd fasttree_runs
ln -s ../paml_format/Fam*paml .
run-FastTree.sh
for file in `ls -1 Fam*paml` ; do unlink $file ; done
cd ..
echo "FastTree runs done"
echo "#####################################"
echo 


echo "#####################################"
echo "preparing codeml input files"
echo 
mkdir codeml_runs
cd codeml_runs
ln -s ../paml_format/Fam*paml  .
ln -s ../fasttree_runs/Fam*tree  .
ln -s ../codeml.ctl.body.template .
codeml-ctl-generator.sh
unlink codeml.ctl.body.template 
cd ..
echo "done preparing"
echo "#####################################"
echo 


echo "#####################################"
echo "running codeml "
echo 
cd  codeml_runs
run-codeml.sh
for file in `ls -1 Fam*paml ` ; do unlink $file ; done
for file in `ls -1 Fam*tree ` ; do unlink $file ; done
cd ..
echo "done codeml"
echo "#####################################"
echo 



