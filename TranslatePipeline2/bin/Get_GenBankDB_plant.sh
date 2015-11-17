#!/bin/bash

## This script helps you to download all plant protein sequences from GenBank and build your reference database.
## To have it running properly:
## copy it into refs/
## ./GetPlantDB.sh

for i in {1..87} 
  do
  wget "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/plant.$i.protein.faa.gz"
  done

mkdir tmp
mv *.faa.gz tmp
cd tmp
gunzip *.gz
cat *.faa > ../Genbank_plant.fas
cd ..
rm -rf tmp