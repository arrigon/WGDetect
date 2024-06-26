#!/bin/bash

# run-macse-aligner.sh 
#
# expects: 1) a number of fasta files 
#          2) macse aligner in /usr/local/macse
#	      
#
# output:  aligned fasta sequences 

for aln in `ls data.in -1 *.fas` 
   do 
     base=`basename $aln .fas`
     echo "### running: $aln"
     #tre=`echo $base | perl -p -e 's/marker001\.Fam\.\d+/tre/'`
     #echo $tre
     #ls ${tre}
     #mlc=${base}.mlc
     
     java -jar -Xmx500m bin/macse_v1.01.jar -prog alignSequences -seq data.in/$aln -out_NT tmp/macse_alignments/$aln.macse.aln
   done
