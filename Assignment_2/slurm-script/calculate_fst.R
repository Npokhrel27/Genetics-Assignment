# Load required package
library(data.table)

# Read SNP file (adjust path if needed)
geno <- fread("/work/agro932/niranjan27/snp_calls_whole_genome/snp_calls.txt", header = FALSE)

# Assign column names
names(geno) <- c("chr", "pos", "ref", "alt", "quality", "depth", paste0("l",1:20))

# Split genotype fields into alleles
geno <- as.data.frame(geno)
for(i in 7:26){
  allele1 <- gsub("/.*", "", geno[, i])
  allele2 <- gsub(".*/", "", geno[, i])
  nm <- names(geno)[i]
  geno[[paste0(nm, "_a1")]] <- allele1
  geno[[paste0(nm, "_a2")]] <- allele2
}

# Replace "." with NA
geno[geno == "."] <- NA

# Calculate allele frequencies
geno$p <- apply(geno[, 27:66], 1, function(x) sum(x == 0, na.rm=TRUE)) / 40
geno$p1 <- apply(geno[, 27:46], 1, function(x) sum(x == 0, na.rm=TRUE)) / 20
geno$p2 <- apply(geno[, 47:66], 1, function(x) sum(x == 0, na.rm=TRUE)) / 20

# Compute Fst
geno$fst <- with(geno, ((p1 - p)^2 + (p2 - p)^2) / (2 * p * (1 - p)))

# Save to CSV
write.csv(geno, "/work/agro932/niranjan27/fst_output/fst_values.csv", row.names = FALSE)

# Plot
png("/work/agro932/niranjan27/fst_output/fst_plot.png", width=1200, height=600)
plot(geno$pos, geno$fst, pch=20, col="blue", xlab="Position", ylab="Fst", main="Per-site Fst")
dev.off()
