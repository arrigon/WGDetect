#!/bin/bash

# creates codeml.ctl files for each pair of .paml and .tre files
#
# expects: 1) .paml file  
#          2) .tre file 
#          3) codeml.ctl.body.template file with all parameter settings 
#	      but without the first 3 lines (seqfile, treefile, outfile)
#
# output: .ctl file for each pair

for paml in `ls -1 *.paml` 
   do 
     base=`basename $paml .paml`
     echo $base
     #tree=`echo $base | perl -p -e 's/marker001\.Fam\.\d+/tre/'`
     tree="${base}.tree"
     #echo $tree
     #ls ${tree}
     mlc="${base}.mlc"

     echo "      seqfile = $paml"   >   ${base}.ctl
     echo "     treefile = $tree"  >>  ${base}.ctl
     echo "      outfile = $mlc"   >>  ${base}.ctl
     cat codeml.ctl.body.template  >>  ${base}.ctl 
   done
