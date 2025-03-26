# === Load required package ===
library(data.table)

# === Read arguments from SLURM (chunk input/output filenames) ===
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_csv <- args[2]
output_plot <- args[3]

# === Read SNP chunk ===
geno <- fread(input_file, header = FALSE)

# === Assign column names ===
names(geno) <- c("chr", "pos", "ref", "alt", "quality", "depth", paste0("l", 1:20))

# === Extract alleles for each individual ===
geno <- as.data.frame(geno)
for(i in 7:26){
  allele1 <- gsub("/.*", "", geno[, i])
  allele2 <- gsub(".*/", "", geno[, i])
  nm <- names(geno)[i]
  geno[[paste0(nm, "_a1")]] <- allele1
  geno[[paste0(nm, "_a2")]] <- allele2
}

# === Replace missing values ===
geno[geno == "."] <- NA

# === Calculate allele frequencies ===
geno$p <- apply(geno[, 27:66], 1, function(x) sum(x == 0, na.rm=TRUE)) / 40
geno$p1 <- apply(geno[, 27:46], 1, function(x) sum(x == 0, na.rm=TRUE)) / 20
geno$p2 <- apply(geno[, 47:66], 1, function(x) sum(x == 0, na.rm=TRUE)) / 20

# === Compute Fst ===
geno$fst <- with(geno, ((p1 - p)^2 + (p2 - p)^2) / (2 * p * (1 - p)))
geno$fst[is.nan(geno$fst)] <- NA  # clean up invalid values

# === Save Fst values to CSV ===
write.csv(geno[, c("chr", "pos", "fst")], output_csv, row.names = FALSE)

# === Plot Fst ===
png(output_plot, width=1200, height=600)
plot(geno$pos, geno$fst, pch=20, col="blue", xlab="Position", ylab="Fst", main=paste("Fst -", basename(input_file)))
dev.off()
