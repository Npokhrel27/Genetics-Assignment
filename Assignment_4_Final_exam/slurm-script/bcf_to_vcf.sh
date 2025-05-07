# Aligning the FASTQ files with reference genome and indexing the resulting BAM files

#!/usr/bin/env bash
#SBATCH -D /work/agro932/niranjan27/Final_exam          # working dir
#SBATCH -o /work/agro932/niranjan27/logs/align-%A_%a.out
#SBATCH -e /work/agro932/niranjan27/logs/align-%A_%a.err
#SBATCH -J align_wg
#SBATCH -p schnablelab
#SBATCH --time=96:00:00
#SBATCH --array=0-9                      # 10 samples indexes 0-9
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --mail-user=npokhrel3@huskers.unl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module load bwa samtools

# paths
REF=/work/agro932/niranjan27/Reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa
FASTQ_DIR=/work/agro932/niranjan27/Final_exam/fastq
OUT_DIR=/work/agro932/niranjan27/Final_exam/aligned
mkdir -p "$OUT_DIR"

# sample list (ten SRR IDs)
SAMPLES=(
  SRR5725631 SRR5725632 SRR5725633 SRR5725634 SRR5725635
  SRR5725636 SRR5725637 SRR5725641 SRR5725642 SRR5725643
)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# input FASTQs
R1=${FASTQ_DIR}/${SAMPLE}.sralite.1_1.fastq
R2=${FASTQ_DIR}/${SAMPLE}.sralite.1_2.fastq
[[ -f $R1 && -f $R2 ]] || { echo "missing FASTQ for $SAMPLE"; exit 1; }

# reference index (only task 0 checks/builds)
if [[ $SLURM_ARRAY_TASK_ID -eq 0 && ! -f ${REF}.bwt ]]; then
  echo "$(date)  Indexing reference"
  bwa index "$REF"
fi
wait   # guarantee index exists before others start

# align, sort, index BAM
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"
SAMTOOLS_THREADS=$(( SLURM_CPUS_PER_TASK - 1 ))

bwa mem -t "$SLURM_CPUS_PER_TASK" -R "$RG" "$REF" "$R1" "$R2" \
  | samtools sort -@ "$SAMTOOLS_THREADS" -o "${OUT_DIR}/${SAMPLE}.sorted.bam" -

samtools index "${OUT_DIR}/${SAMPLE}.sorted.bam"

echo "$(date) finished $SAMPLE"
```

## SNP calling from the bam file created before and creating the bcf files.

```{bash, eval=FALSE}


#!/usr/bin/env bash
#SBATCH -D /work/agro932/niranjan27/Final_exam
#SBATCH -o /work/agro932/niranjan27/Final_exam/logs/snp-%A_%a.out
#SBATCH -e /work/agro932/niranjan27/Final_exam/logs/snp-%A_%a.err
#SBATCH -J snp_call_array                 #job name changed
#SBATCH -p schnablelab
#SBATCH --time=72:00:00
#SBATCH --array=0-9                       # 10 BAMs indexes 0-9
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --mail-user=npokhrel3@huskers.unl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module load samtools bcftools

# Directories
REF=/work/agro932/niranjan27/Reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa
BAMS_DIR=/work/agro932/niranjan27/Final_exam/aligned
OUTDIR=/work/agro932/niranjan27/snp_calls_whole_genome
LOGDIR=/work/agro932/niranjan27/Final_exam/logs

mkdir -p "$OUTDIR" "$LOGDIR"

# BAM list & sample name
BAMS=(${BAMS_DIR}/*.sorted.bam)
BAM=${BAMS[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename "$BAM" .sorted.bam)

# index BAM if missing
[[ -f ${BAM}.bai ]] || samtools index "$BAM"

# SNP calling
bcftools mpileup -Ou -f "$REF" -a AD,DP "$BAM" \
| bcftools call -mv -Ob -o "$OUTDIR/${SAMPLE}.bcf"

bcftools index "$OUTDIR/${SAMPLE}.bcf"

echo "$(date) Finished SNP calling for ${SAMPLE}"

```


## Converting the bcf file of temperate maize lines to vcf files

```{bash, eval=FALSE}
module load bcftools vcftools
REF=/work/agro932/niranjan27/Reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa

mkdir -p /work/agro932/niranjan27/Final_exam/merged_fst
cd       /work/agro932/niranjan27/Final_exam/merged_fst

bcftools merge /work/agro932/niranjan27/Final_exam/bcf_files/*.bcf \
-Oz -o temperate.raw.vcf.gz
bcftools index temperate.raw.vcf.gz

# Normalize temperate to same referene genome as tropical
bcftools norm -f "$REF" -m -both \
temperate.raw.vcf.gz -Oz -o temperate.norm.vcf.gz
bcftools index -f temperate.norm.vcf.gz

# Normalizing the tropical data with same reference genome to remmove any uncertainty

# bgzip + index the raw file once
bcftools view -Oz -o BGEM_trop.vcf.gz \
/work/agro932/niranjan27/Final_exam/Tropical_vcf/BGEM_15founders.recode.vcf
bcftools index -f BGEM_trop.vcf.gz

# normalise
bcftools norm -f "$REF" -m -both \
BGEM_trop.vcf.gz -Oz -o tropical.norm.vcf.gz
bcftools index -f tropical.norm.vcf.gz


# Since we had 15 tropical files, keeping only 10 vcf files of tropical one
# list the 15 IDs
bcftools query -l /work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical.norm.vcf.gz \
> tropical_all.txt

nano tropical_all.txt            # delete any 5 IDs, save as tropical10.txt

bcftools view -S tropical10.txt -Oz \
-o /work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical10.norm.vcf.gz \
/work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical.norm.vcf.gz
bcftools index -f /work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical10.norm.vcf.gz


# Merging the tropical and temperate individuals to create single vcf file

mkdir -p /work/agro932/niranjan27/Final_exam/merged_fst

bcftools merge -Oz \
-o /work/agro932/niranjan27/Final_exam/merged_fst/maize_merged.vcf.gz \
/work/agro932/niranjan27/Final_exam/temperate_vcf/temperate.norm.vcf.gz \
/work/agro932/niranjan27/Final_exam/Tropical_vcf/tropical10.norm.vcf.gz

bcftools index /work/agro932/niranjan27/merged_fst/maize_merged.vcf.gz
