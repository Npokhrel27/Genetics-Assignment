#!/bin/bash
#SBATCH -D /work/agro932/niranjan27/
#SBATCH -o /work/agro932/niranjan27/logs/align-wg-stdout-%A_%a.txt
#SBATCH -e /work/agro932/niranjan27/logs/align-wg-stderr-%A_%a.txt
#SBATCH -J align_wholegenome
#SBATCH -t 96:00:00
#SBATCH -p schnablelab
#SBATCH --array=1-20
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --mail-user=npokhrel3@huskers.unl.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# === Load required tools ===
module load bwa samtools

# === Paths ===
REF_WG="/work/agro932/niranjan27/Whole_reference_genome/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa"
MAIZE_DIR="/work/agro932/niranjan27/maize"
MEXICANA_DIR="/work/agro932/niranjan27/mexicana"
OUTPUT_DIR="/work/agro932/niranjan27/aligned_wholegenome"
SAMPLE_LIST="/work/agro932/niranjan27/sample_list.txt"

# === Index the reference genome if not already done ===
if [[ ! -f "${REF_WG}.bwt" ]]; then
    echo "🧬 Indexing reference genome..."
    bwa index "$REF_WG"
fi

# === Create output directory if it doesn't exist ===
mkdir -p "$OUTPUT_DIR"

# === Get sample info for this array task ===
SAMPLE_INFO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
PREFIX=$(echo "$SAMPLE_INFO" | awk '{print $1}')
SAMPLE=$(echo "$SAMPLE_INFO" | awk '{print $2}')

# === Set FASTQ read paths ===
if [[ "$PREFIX" == "maize" ]]; then
    READ_DIR="$MAIZE_DIR"
else
    READ_DIR="$MEXICANA_DIR"
fi

R1="${READ_DIR}/${SAMPLE}.lite.1_1.fastq"
R2="${READ_DIR}/${SAMPLE}.lite.1_2.fastq"

# === Confirm input files exist ===
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "❌ Missing FASTQ files for $SAMPLE"
    exit 1
fi

# === Define output BAM filenames ===
OUT_BAM="${OUTPUT_DIR}/${PREFIX}_${SAMPLE}_wg.bam"
SORTED_BAM="${OUTPUT_DIR}/${PREFIX}_${SAMPLE}_wg_sorted.bam"

# === Alignment process ===
echo "🔄 Aligning $SAMPLE..."
bwa mem -t 16 "$REF_WG" "$R1" "$R2" | samtools view -bSh - > "$OUT_BAM"
samtools sort -@ 8 "$OUT_BAM" -o "$SORTED_BAM"
samtools index "$SORTED_BAM"
rm "$OUT_BAM"

echo "✅ Done aligning $SAMPLE"
