
# Selection scan and Genome-wide association analysis for nitrogen use efficiency related trait in maize.

This repository contains all the files and codes necessary for selection scan (Fst, nucleotide diversity, Tajima's D) and GWAS analysis.

There are three main directories:

**1. cache**/
Contains the plots and csv files obtained from GWAS analysis.

**2. profilling**/
Contains an RMarkdown (.Rmd) file with both Bash and R code used to:

  Create VCF files from FASTQ files

  Perform selection scan and GWAS analysis
  
  A PDF version of the code is also included in this directory.

**3. slurm-script**/
Contains the SLURM Bash scripts used for SNP calling from FASTQ files.