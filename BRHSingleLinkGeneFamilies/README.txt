### Purpose of the pipeline:
- Define sequence clusters (a.k.a. gene families if provided with CDS) using Single Linkage Clustering
- Based on all-vs-all blast searches
- The granularity of clusters can be adjusted in blast search parameters

### Installation and usage
1. Get blast at work

2. Add all input sequences in ./data.in
  Add either a sinlge (e.g. all genbank accessions) or several (e.g. if analysing a set of transcriptomes) FASTA files (DNA by default)
  IF WORKING with proteins: you need to adjust the blast search parameters (see below)

3. Check the blast config file in ./config/blastline.conf : 
  Specify the makeblastdb and blastall lines according to your OS / analysis (DNA or protein)
  Keep input/output and number of cores parameters as undefined, will be taken care of the by pipeline

  NOTE: you can adapt the blast line according to the desired granularity, the defaults are e-value = 1e-3 and %similarity = 30%

4. Run the pipeline
perl BRH_SingleLinkFamilies_build.pl minlen minsim NCPU MaxUse

minlen = length over which sequences must match with each other
minsim = minimum percentage of similarity between matching sequences
NCPU = Number of cores to be used by each process
MaxUse = Maximal load allowed on the computing server.


### Details
The pipeline will:
1. Prefix the accession names, according to fasta name (if using one fasta per transcriptome) or one unique fasta including all sequences to be clustered (e.g. if clustering all data from genbank). The prefix is added to ensure we get unique sequence headers (e.g. when analysing transcriptomes where sequences are just numbered).
2. Merge all fasta into a single file
3. Run blast as follows:
  3.a	blast all-vs-all (blast line can be configured using config/blastline.conf; useful for adapting among blast flavors)
  3.b	parse output to keep best hits (according to minsim and minlen)
4. Perform single linkage clustering using pairwise best hits
5. Produce fasta files of each gene family

IMPORTANT:
Sequence headers are modified with a prefix; might need to be cleaned in further steps.

NEXT STEP: OrthoMCL (to refine this clustering to keep Orthologs)
Note that gene family analyses (e.g. detection of ancioent genome duplications, etc) can be initiated using the current clusters.