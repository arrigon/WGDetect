#!/bin/bash

## This script helps you to download all plant protein sequences from GenBank and build your reference database.
## To have it running properly:
## copy it into refs/
## ./Get_GenBankDB_vertebrate_other.sh

for i in {1..100} 
  do
  wget "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_other/vertebrate_other.$i.protein.faa.gz"
  done

mkdir tmp
mv *.faa.gz tmp
cd tmp
gunzip *.gz
cat *.faa > ../Genbank_vertebrate_other.fas
cd ..
rm -rf tmp