

install.packages("rMVP")
library(rMVP)

MVP.Data(fileHMP="Data/BGEM_MMC_imputed_181entries_62077markers_perse.hmp.txt",
         filePhe="Data/Phenotypic_data.csv",
         sep.phe=";",
         fileKin=FALSE,
         filePC=FALSE,
         #maxLine=10000,
         out="Niranjan.mvp.hmp"
)

pheno<- read.csv("Assignment_3_midterm/Data/Phenotypic_data.csv")


genotype <- attach.big.matrix("/home/agro932/niranjan27/Genetics-Assignment/Assignment_3_midterm/Niranjan.mvp.hmp.geno.desc")
View(genotype)
phenotype <- read.table("/home/agro932/niranjan27/Genetics-Assignment/Assignment_3_midterm/Niranjan.mvp.hmp.phe", head=TRUE)
View(phenotype)
map <- read.table("/home/agro932/niranjan27/Genetics-Assignment/Assignment_3_midterm/Niranjan.mvp.hmp.geno.map" , head = TRUE)


head(phenotype)
# Set output directory
setwd("/home/agro932/niranjan27/Genetics-Assignment/Final_term/Plots")

trait <- phenotype[c("Taxa", "SIL_LN")] #can change the trait here!!!!!
imMVP <- MVP(
  phe=trait,          #NA is acceptable in phenotype
  geno=genotype,
  map=map,             #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  #CV.GLM=Covariates,     #if you have environmental covariates, please keep all 'CV.*' open
  #CV.MLM=Covariates,
  #CV.FarmCPU=Covariates,
  nPC.GLM=5,              #if you have added PCs into covariates, please keep there closed
  nPC.MLM=3,              #if you don't want to add PCs as covariates, please comment out the parameter instead of setting it to 0.
  nPC.FarmCPU=3,
  maxLine=10000,          #smaller value would reduce the memory cost
  #ncpus=10,
  vc.method="BRENT",      #only works for MLM
  method.bin="static",    # "FaST-LMM", "static" (#only works for FarmCPU)
  threshold=0.05,
  method=c("FarmCPU"),   #can adjust GLM and MLM here
  file.output=c("pmap", "pmap.signal", "plot", "log")
)


