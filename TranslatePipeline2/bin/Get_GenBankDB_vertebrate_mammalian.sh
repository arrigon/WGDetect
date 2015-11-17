#!/bin/bash

## This script helps you to download all plant protein sequences from GenBank and build your reference database.
## To have it running properly:
## copy it into refs/
## ./Get_GenBankDB_vertebrate_mammalian.sh

for i in {1..433} 
  do
  wget "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian.$i.protein.faa.gz"
  done

mkdir tmp
mv *.faa.gz tmp
cd tmp
gunzip *.gz
cat *.faa > ../Genbank_vertebrate_mammalian.fas
cd ..
rm -rf tmp