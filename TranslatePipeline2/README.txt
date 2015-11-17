### Purpose of the pipeline:
- Compare input DNA sequences to a reference protein database
- Keep only DNA sequences with hits in the protein DB
- Extracts, Frames and Translates the filtered DNA sequences 
  to produce outputs that are ready to be clustered into gene familes or aligned (if providing sequences from a single gene family).
- Sequences that span over several genes are sliced into single CDS (that hopefully span the complete protein)
  the fasta headers are updated accordingly.
- Produces three outputs:
  *.cds.fas = framed DNA sequence (CDS only)
  *.pep.fas = translated protein
  *.all.fas = DNA sequence (CDS + introns), WARNING: not in frame over complete sequence.


### Installation and usage
1. Install Exonerate, if the copy in ./bin is not working on your OS.


2. Install Blast
  The present pipeline uses the older version of blast (i.e. based on formatdb and not makeblastdb)
  It is super easy to run any newer version: 
  update the commands into ./Translate_build.pl ; at line 40 (use makeblastdb instead of formatdb)
  update the commands into bin/RunParseBlastx.pl; at line 24 (make sure to use blastx and get the tablular output)


3. Add the proper protein database in ./refs
  ./refs : Add a fasta file that includes all the candidate proteins (e.g. use the Phytozome annotation of your most closely related genome)
  This file MUST be named as follows: Genbank_dbase.fas, where "dbase" can be anything (and will be given as a parameter when launching the pipeline, see below)
  WARNING: the bigger the database, the better the results but the longer the search...
  N.B. you find a script to download all GenBank plant protins in bin/Get_GenBankDB_*.sh, run it from refs/ to build your database if needed.
  
  Another point: these reference sequences contain the "J" character. This causes exonerate to crash (and to loose those cds that hit against such a reference); 
  since J is not an amino acid. J actually stands for "I or L". 
  
  One work around is to replace J with I or L, using a sed command, e.g. 
  
  sed -i -- 's/J/I/g' Genbank_vertebrate.fas
  
  We do not need exact matches with exonerate, but rather look for appropriate codon phasing. Hence, replacing J with I or L is probably fine.
  
  
4. Add your input sequences in ./data.in 
  Either add one global fasta file or several of them, the script will iterate over everything and produce separate outputs if given several files

  NOTE: It is more productive to slice big input files into smaller chuncks and let the pipeline run in parallel on the wished number of cores
  You can slice your dataset using bin/SliceFasta.pl inputfile outfolder nslices

  Typically, use the following command (assuming that the inputfile is in the root of the pipeline folder):
  perl bin/SliceFasta.pl MyBigSingleInput.fasta data.in 2*NbCPUs


5. Run the pipeline (run everything from the root folder)
  perl Translate_build.pl dbase NbCPUS MaxUse
  - dbase = database to use (either plant, invertebrate, vertebrate_mammalian, vertebrate_other or any name of your own: use this exact spelling)
  - NbCPUS = Number of CPUs to allocate per blast search
  - MaxUse = Maximum allowed load on the server 


6. Get the results in ./data.out (added on the fly by the pipeline)
  *.cds.fas = framed DNA sequence (CDS only)
  *.pep.fas = translated protein
  *.all.fas = DNA sequence (CDS + introns), WARNING: not in frame over complete sequence.


### Misc
(Old version, not in use anymore but maybe usefull for GeneWise users)
1. Get Genewise at work; 
  Wise2, can be obtained via synaptic manager, you have to set up environment variables though.
  Add this line to your ~/.bashrc : (change /home/arrigon/soft/ for you actual install path)
  export WISECONFIGDIR=/home/arrigon/soft/wise2.2.0/wisecfg/
