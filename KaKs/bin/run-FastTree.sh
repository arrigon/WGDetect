#!/bin/bash

# runs FastTree on all Fam*paml files in the directory

for paml in `ls -1 Fam*paml` 
   do 
     base=`basename $paml .paml`
     echo $base
     ../../bin/FastTree  -gtr -gamma -nt $paml >  ${base}.tree 2> ${base}.log
   done

