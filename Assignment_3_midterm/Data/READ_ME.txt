# README
## GENERAL INFORMATION 

DATA TITLE: Phenotypic data of nitrogen use efficiency related traits and genotypic data of BGEM lines

PROJECT TITLE: Phenotypic and genome-wide association analyses for nitrogen use efficiency related traits in maize (Zea mays L.) exotic introgression lines.

DATA ABSTRACT: These data set (phenotypic and genotypic data) can be used in the genetics and plant breeding area. Currently, many projects have been using analyzes that associate the phenotype of interest with the genotype. Thus, it allows identifying possible regions of the genome that have the greatest effect on the trait. Our data comprise 181 doubled haploid (DH) lines derived from crosses between landraces from the Germplasm Enhancement of Maize (BGEM lines) project and two inbreds, PHB47 and PHZ51. These DH lines were genotyped using 62,077 molecular markers (single nucleotide polymorphisms - SNPs) using genotyping-by-sequencing (GBS). The same lines from the per se trials were used as parental lines for the testcross field trials. Plant height (PHT, cm), days to female flowering (SLK, days), days to male flowering (ANT, days), anthesis to silking interval (ASI), and grain yield (YLD, tons per hectare) were collected from high (HN) and low N (LN) conditions in two environments for testcrosses (Ames 2015B, Nashua 2015) and three environments for per se trials (Ames 2014, Ames 2015 and Nashua 2015). Our data might be used for Genome Wide Association Studies (GWAS) and Selection (GWS). Furthermore, molecular characterization and genetic diversity might be performed. 

AUTHORS:

	Author: Darlene L. Sanchez
	Institution: Texas A&M University
	Email: darlene.sanchez@ag.tamu.edu

    Author: Alice Silva Santana
	Institution: Iowa State University
	Email: a.santana@iastate.edu

    Author: Thomas Lubberstedt
	ORCID: https://orcid.org/0000-0002-0526-0798
	Institution: Iowa State University
	Email: thomasl@iastate.edu

	Corresponding author: Thomas Lubberstedt

ASSOCIATED PUBLICATIONS:
        - Publication coming soon

COLLECTION INFORMATION:

	Time period(s):  During the summer of 2013 and 2014
	Location(s): Iowa State University


### FILE DIRECTORY 
#### FILE LIST

1. Data pheno and Geno
	a.  BGEM_MMC_imputed_180entries_62077markers_testcross.hmp, is the file of genotypic data from the tescrosses.
	b.  BGEM_MMC_imputed_181entries_62077markers_perse.hmp, is the file of genotypic data derived from the BGEM lines.
        c.  Phenotypic data, is the phenotypic data of perse and tescross genotypes.


### CODEBOOK 

a. BGEM_MMC_imputed_180entries_62077markers_testcross.hmp
Number Of Variables/Columns: 187
Number Of Cases/Rows: 62077

b.  BGEM_MMC_imputed_181entries_62077markers_perse.hmp
Number Of Variables/Columns: 227
Number Of Cases/Rows: 62077

c. Phenotypic data
Number Of Variables/Columns: 13
Number Of Cases/Rows: 905


#### VARIABLES

a. BGEM_MMC_imputed_180entries_62077markers_testcross.hmp
| All columns | is the DH lines | 
| All Rows | is the SNPs |  

b.  BGEM_MMC_imputed_181entries_62077markers_perse.hmp
| All columns | is the DH lines | 
| All Rows | is the SNPs |

c.Phenotypic data.csv
Environment - column with the genetic background (per se or testcross) plus environment where the genotypes where evaluated, for example: "Testcross_Ames2015A" means the results from the testcross trial evaluated in Ames in the year of 2015 | YLD_HN - Yield trait under HN| PHT_HN - Plant height under HN | ASI_HN - Anthesis-silking interval under HN | ANT_HN - Anthesis-silking interval under HN | SLK_HN - days to female flowering under HN | YLD_LN - Yield under LN | PHT_LN - Plant Height under LN | ASI_LN - Antesis-silking interval | ANT_LN - days to male flowering under LN | SLK_LN - days to female flowering under LN | 
| All Rows | is the testcross hybrids or per se lines


## METHODS AND MATERIALS
### DATA COLLECTION METHODS 

a. BGEM_MMC_imputed_180entries_62077markers_testcross.hmp
The DNA extraction was done using the standard International Maize and Wheat Improvement Center (CIMMYT) laboratory protocol CIMMYT, 2005). Genotyping was carried out using the genotyping-by-sequencing (GBS) markers. GBS data were generated at the Cornell Institute for Genomic Diversity (IGD) laboratory. After filtering out markers with more than 25% missing data, below 2.5% minor allele frequency, and monomorphic markers, 247,775 markers were left for further analyses. For markers at the same genetic position (0 cM distance), only one marker was randomly selected. The final number of markers used for further analyses was 62,077 markers distributed across all 10 chromosomes.

CIMMYT. (2005). Laboratory Protocols: CIMMYT Applied Molecular Genetics Laboratory (Third Edit). Mexico, D.F: CIMMYT. Retrieved from papers2://publication/uuid/B00DBE67-68B3-4845-BA5F-F687141F2DF7

b.  BGEM_MMC_imputed_181entries_62077markers_perse.hmp
The DNA extraction was done using the standard International Maize and Wheat Improvement Center (CIMMYT) laboratory protocol CIMMYT, 2005). Genotyping was carried out using the genotyping-by-sequencing (GBS) markers. GBS data were generated at the Cornell Institute for Genomic Diversity (IGD) laboratory. After filtering out markers with more than 25% missing data, below 2.5% minor allele frequency, and monomorphic markers, 247,775 markers were left for further analyses. For markers at the same genetic position (0 cM distance), only one marker was randomly selected. The final number of markers used for further analyses was 62,077 markers distributed across all 10 chromosomes.

CIMMYT. (2005). Laboratory Protocols: CIMMYT Applied Molecular Genetics Laboratory (Third Edit). Mexico, D.F: CIMMYT. Retrieved from papers2://publication/uuid/B00DBE67-68B3-4845-BA5F-F687141F2DF7



c. Phenotypic data.csv
The phenotypic data were collected in the field by measuring the BGEM lines and their testcrosses.
Male flowering (days) - days until male flowering;
Female flowering (days) - days until female flowering;
Anthesis-silking interval (days) - flowering interval between male and female;
Plant height (cm) - Plant height from ground to tassel;
Yield (tons per hectare) - computed after moisture content was adjusted to 15.50%.


### DATA PROCESSING METHODS 

a. BGEM_MMC_imputed_180entries_62077markers_testcross.hmp
Across the samples assessed a total of 955,690 SNP markers were generated, of these 62,077 SNP markers were successfully aligned to the B73 RefGen_v4 and passed throught the quality control (filtering out markers with more than 25% missing data, below 2.5% minor allele frequency, and monomorphic markers).

b.  BGEM_MMC_imputed_181entries_62077markers_perse.hmp
Across the samples assessed a total of 955,690 SNP markers were generated, of these 62,077 SNP markers were successfully aligned to the B73 RefGen_v4 and passed throught the quality control (filtering out markers with more than 25% missing data, below 2.5% minor allele frequency, and monomorphic markers).

c.  Phenotypic data
Only the most discrepant data were removed (outliers).


### SOFTWARE
 
Name: TASSEL software
Version: v.5.2.64
URL: https://www.maizegenetics.net/tassel
Developer: Cornell University

Name: R software
Version: 4.3.1
URL: https://www.r-project.org/
Developer: Bell Laboratories (formerly AT&T, now Lucent Technologies) by John Chambers and colleagues



### LICENSING 

This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).



