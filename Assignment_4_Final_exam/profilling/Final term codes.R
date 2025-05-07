
# Bash codes
module load bcftools vcftools
REF=/work/agro932/niranjan27/Reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa

mkdir -p /work/agro932/niranjan27/Final_exam/merged_fst
cd       /work/agro932/niranjan27/Final_exam/merged_fst


# merge temperate bcf into one vcf
bcftools merge /work/agro932/niranjan27/Final_exam/bcf_files/*.bcf \
-Oz -o temperate.raw.vcf.gz
bcftools index temperate.raw.vcf.gz

# Normalize temperate data
bcftools norm -f "$REF" -m -both \
temperate.raw.vcf.gz -Oz -o temperate.norm.vcf.gz
bcftools index -f temperate.norm.vcf.gz

# Normalize tropical data

# bgzip + index the raw file once
bcftools view -Oz -o BGEM_trop.vcf.gz \
/work/agro932/niranjan27/Final_exam/Tropical_vcf/BGEM_15founders.recode.vcf
bcftools index -f BGEM_trop.vcf.gz

# normalise
bcftools norm -f "$REF" -m -both \
BGEM_trop.vcf.gz -Oz -o tropical.norm.vcf.gz
bcftools index -f tropical.norm.vcf.gz


# Keep only 10 vcf files of tropical one
# list the 15 IDs
bcftools query -l /work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical.norm.vcf.gz \
> tropical_all.txt

# edit this file down to 10 lines you want → tropical10.txt
nano tropical_all.txt            # delete any 5 IDs, save as tropical10.txt

# slice the VCF down to those 10 founders
bcftools view -S tropical10.txt -Oz \
-o /work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical10.norm.vcf.gz \
/work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical.norm.vcf.gz
bcftools index -f /work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical10.norm.vcf.gz


# Merging the tropical and temperate individuals in same vcf file

mkdir -p /work/agro932/niranjan27/Final_exam/merged_fst

bcftools merge -Oz \
-o /work/agro932/niranjan27/Final_exam/merged_fst/maize_merged.vcf.gz \
/work/agro932/niranjan27/Final_exam/temperate_vcf/temperate.norm.vcf.gz \
/work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical10.norm.vcf.gz

bcftools index /work/agro932/niranjan27/merged_fst/maize_merged.vcf.gz


# Calculate per SNP Fst 

vcftools --gzvcf maize_merged.vcf.gz \
--weir-fst-pop temperate.txt \
--weir-fst-pop tropical.txt \
--out temp_vs_trop_10vs10

# After filtering, kept 20 out of 20 Individuals
# Outputting Weir and Cockerham Fst estimates.
# Weir and Cockerham mean Fst estimate: 0.39787
# Weir and Cockerham weighted Fst estimate: 0.52008
# After filtering, kept 50278096 out of a possible 50278096 Sites



# Plotting Fst

library(data.table)
library(ggplot2)


# Load data
fst_data <- fread("/work/agro932/niranjan27/Final_exam/merged_fst/temp_vs_trop_win50k.windowed.weir.fst")

# Plot histogram
hist(fst_data$MEAN_FST,
     breaks = 50,
     col = "steelblue",
     main = "Histogram of FST values",
     xlab = "Mean FST",
     ylab = "Number of windows")

###### Plotting fst

# Clean: replace NA and negative FST values
fst_data <- fst_data %>%
  filter(!is.na(MEAN_FST)) %>%
  mutate(MEAN_FST = ifelse(MEAN_FST < 0, 0, MEAN_FST))

# Ensure CHROM is treated as an ordered factor
fst_data$CHROM <- factor(fst_data$CHROM, levels = mixedsort(unique(as.character(fst_data$CHROM))))

# Calculate chromosome lengths and offsets
chr_lengths <- fst_data %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(BIN_START)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

# Merge back to data and compute cumulative position
fst_data <- fst_data %>%
  left_join(chr_lengths, by = "CHROM") %>%
  mutate(cum_pos = BIN_START + chr_start)

# Midpoint for chromosome labels
axis_df <- chr_lengths %>%
  mutate(mid = chr_start + chr_len / 2)


ggplot(fst_data, aes(x = cum_pos, y = MEAN_FST)) +
  geom_point(size = 0.3, alpha = 0.6, color = "blue") +
  scale_x_continuous(
    label = axis_df$CHROM,
    breaks = axis_df$mid,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Sliding Window FST (Tropical vs Temperate Maize)",
    x = "Chromosome",
    y = "Mean FST"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, size = 8),
    plot.title = element_text(hjust = 0.5)
  )









#### CALCULATING NUCLEOTIDE DIVERSSITY

## Run VCF tool for nucleotide diversity

# Tropical
vcftools --gzvcf /work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical10.norm.vcf.gz \
--keep tropical.txt \
--window-pi 50000 \
--window-pi-step 25000 \
--out tropical_pi

# Temperate
vcftools --gzvcf /work/agro932/niranjan27/Final_exam/temperate_vcf/temperate.norm.vcf.gz \
--keep temperate.txt \
--window-pi 50000 \
--window-pi-step 25000 \
--out temperate_pi


#Calculate nucleotide diversity. 

trop <- read.table("/work/agro932/niranjan27/Final_exam/merged_fst/tropical_pi.windowed.pi", header=TRUE)
temp <- read.table("/work/agro932/niranjan27/Final_exam/merged_fst/temperate_pi.windowed.pi", header=TRUE)

mean(trop$PI, na.rm=TRUE)  # Tropical π
mean(temp$PI, na.rm=TRUE)  # Temperate π






# Plotting nucleotide diversity in R

library(data.table)
library(ggplot2)
library(tidyverse)

# Load data
trop <- fread("/work/agro932/niranjan27/Final_exam/merged_fst/tropical_pi.windowed.pi")
temp <- fread("/work/agro932/niranjan27/Final_exam/merged_fst/temperate_pi.windowed.pi")


# Add population info
trop$pop <- "Tropical"
temp$pop <- "Temperate"

# Combine
pi_all <- rbind(trop, temp)

# Convert CHROM to numeric
pi_all$CHROM <- as.numeric(as.character(pi_all$CHROM))

# Sort by chromosome and position
pi_all <- pi_all %>% arrange(CHROM, BIN_START)

# Calculate cumulative position for plotting
chr_lengths <- pi_all %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(BIN_END)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

pi_all <- left_join(pi_all, chr_lengths, by = "CHROM")
pi_all$position <- pi_all$BIN_START + pi_all$chr_start

# Plot as a continuous line graph
ggplot(pi_all, aes(x = position, y = PI, color = pop)) +
  geom_line(alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("Tropical" = "red", "Temperate" = "blue")) +
  labs(x = "Genomic Position", y = "Nucleotide Diversity (π)",
       title = "Nucleotide Diversity Across Genome",
       color = "Population") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())










# Calculating Tajima's D value


# Temperate
vcftools --gzvcf /work/agro932/niranjan27/Final_exam/temperate_vcf/temperate.norm.vcf.gz \
--TajimaD 10000 \
--out temperate_tajima

# Tropical
vcftools --gzvcf /work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical10.norm.vcf.gz \
--TajimaD 10000 \
--out tropical_tajima


# Plotting tajima's D values 

library(data.table)
library(ggplot2)
library(dplyr)
library(gtools)

# TROPICAL LINES 
# Load Tajima's D data
tajima <- fread("/work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical_tajima.Tajima.D")
colnames(tajima) <- c("CHROM", "BIN_START", "N_SNPS", "TajimaD")

# Remove extreme values for better visualization (optional)
tajima <- tajima %>% filter(TajimaD > -5, TajimaD < 5)

# Ensure CHROM is treated as a factor and sorted naturally
tajima$CHROM <- factor(tajima$CHROM, levels = mixedsort(unique(tajima$CHROM)))

# Calculate chromosome lengths and offsets
chr_lengths <- tajima %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(BIN_START)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

# Merge back to get cumulative positions
tajima <- tajima %>%
  left_join(chr_lengths, by = "CHROM") %>%
  mutate(cum_pos = BIN_START + chr_start)

# Midpoint of each chromosome for axis labeling
axis_df <- chr_lengths %>%
  mutate(mid = chr_start + chr_len / 2)

# Plot
ggplot(tajima, aes(x = cum_pos, y = TajimaD)) +
  geom_point(size = 0.4, color = "darkred", alpha = 0.6) +
  scale_x_continuous(label = axis_df$CHROM, breaks = axis_df$mid) +
  labs(
    x = "Chromosome",
    y = "Tajima's D",
    title = "Tajima's D Across Genome (Tropical Maize)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, size = 8),
    plot.title = element_text(hjust = 0.5)
  )



# TEMPERATE LINES

# Load Tajima's D data for temperate maize
tajima <- fread("/work/agro932/niranjan27/Final_exam/temperate_vcf/temperate_tajima.Tajima.D")
colnames(tajima) <- c("CHROM", "BIN_START", "N_SNPS", "TajimaD")

# Remove extreme values (optional)
tajima <- tajima %>% filter(TajimaD > -5, TajimaD < 5,N_SNPS >= 10)

# Ensure chromosomes are treated as ordered factors
tajima$CHROM <- factor(tajima$CHROM, levels = mixedsort(unique(tajima$CHROM)))

# Calculate cumulative genome positions
chr_lengths <- tajima %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(BIN_START)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

# Merge and compute cumulative positions
tajima <- tajima %>%
  left_join(chr_lengths, by = "CHROM") %>%
  mutate(cum_pos = BIN_START + chr_start)

# Axis labels
axis_df <- chr_lengths %>%
  mutate(mid = chr_start + chr_len / 2)

# Plot
ggplot(tajima, aes(x = cum_pos, y = TajimaD)) +
  geom_point(size = 0.4, color = "blue", alpha = 0.6) +
  scale_x_continuous(label = axis_df$CHROM, breaks = axis_df$mid) +
  labs(
    x = "Chromosome",
    y = "Tajima's D",
    title = "Tajima's D Across Genome (Temperate Maize)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, size = 8),
    plot.title = element_text(hjust = 0.5)
  )




######

# Read the VCFtools output
tajima <- fread("/work/agro932/niranjan27/Final_exam/temperate_vcf/temperate_tajima.Tajima.D")

# Rename columns for clarity
colnames(tajima) <- c("CHROM", "BIN_START", "N_SNPS", "TajimaD")

# Filter extreme D values and low-SNP windows (optional but recommended)
tajima <- tajima %>%
  filter(TajimaD > -5, TajimaD < 5, N_SNPS >= 10)

# Keep only chromosomes 1–10
tajima <- tajima %>%
  filter(CHROM %in% as.character(1:10))

# Sort chromosomes naturally (1,2,...10 instead of 1,10,2,...)
tajima$CHROM <- factor(tajima$CHROM, levels = mixedsort(unique(tajima$CHROM)))

# Get length of each chromosome based on max BIN_START
chr_lengths <- tajima %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(BIN_START)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

# Merge to get cumulative position for each window
tajima <- tajima %>%
  left_join(chr_lengths, by = "CHROM") %>%
  mutate(cum_pos = BIN_START + chr_start)

# For chromosome labels on the plot
axis_df <- chr_lengths %>%
  mutate(mid = chr_start + chr_len / 2)


ggplot(tajima, aes(x = cum_pos, y = TajimaD)) +
  geom_point(size = 0.4, color = "blue", alpha = 0.6) +
  scale_x_continuous(
    label = axis_df$CHROM,
    breaks = axis_df$mid,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Chromosome",
    y = "Tajima's D",
    title = "Tajima's D Across Genome (Temperate Maize)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, size = 8),
    plot.title = element_text(hjust = 0.5)
  )





