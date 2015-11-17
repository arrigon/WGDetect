#!/bin/bash

# runs codeml on all Fam*ctl files in the directory
# and changes some file names

echo "${PWD##*/}"
     
for ctl in `ls -1 Fam*ctl` 
   do 
     base=`basename $ctl .ctl`
     echo $base
     mkdir ${base}.codmlrun
     cd ${base}.codmlrun
     ln -s ../$ctl .
     ln -s ../${base}.tree .
     ln -s ../${base}.paml .
     ./../../../bin/codeml $ctl
#      mv 2NG.dS ${base}.2NG.dS 
#      mv 2NG.dN ${base}.2NG.dN
     # rm -rf 4fold.nuc lnf rst 2NG.dS 2NG.dN 2NG.t rub rst1
#      unlink ${base}.tree 
#      unlink ${base}.paml
     cd ..
   done

