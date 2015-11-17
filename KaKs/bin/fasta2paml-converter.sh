#!/bin/bash
#
# convert all fasta files in the directory to paml format
# Fam.5611.macse.aln_DNA.fasta
# for fasta  in `ls -1 *.aln` ; do ../../bin/fasta2paml.SZ.pl -c 10 $fasta  > $fasta.paml ; done
for fasta  in `ls -1 *.aln` ; do perl ../../bin/Fas2Phy.pl $fasta $fasta.paml ; done
#