#!/usr/bin/env bash
#SBATCH -D /work/agro932/niranjan27/Final_exam          # working dir
#SBATCH -o /work/agro932/niranjan27/logs/align-%A_%a.out
#SBATCH -e /work/agro932/niranjan27/logs/align-%A_%a.err
#SBATCH -J align_wg
#SBATCH -p schnablelab
#SBATCH --time=96:00:00
#SBATCH --array=0-9                      # 10 samples  → indexes 0-9
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --mail-user=npokhrel3@huskers.unl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module load bwa samtools

# ──────── paths ────────────────────────────────────────────────────
REF=/work/agro932/niranjan27/Reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa
FASTQ_DIR=/work/agro932/niranjan27/Final_exam/fastq
OUT_DIR=/work/agro932/niranjan27/Final_exam/aligned
mkdir -p "$OUT_DIR"

# ──────── sample list (exactly your ten SRR IDs) ───────────────────
SAMPLES=(
  SRR5725631 SRR5725632 SRR5725633 SRR5725634 SRR5725635
  SRR5725636 SRR5725637 SRR5725641 SRR5725642 SRR5725643
)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# ──────── input FASTQs ─────────────────────────────────────────────
R1=${FASTQ_DIR}/${SAMPLE}.sralite.1_1.fastq
R2=${FASTQ_DIR}/${SAMPLE}.sralite.1_2.fastq
[[ -f $R1 && -f $R2 ]] || { echo "❌ missing FASTQ for $SAMPLE"; exit 1; }

# ──────── reference index (only task 0 checks/builds) ──────────────
if [[ $SLURM_ARRAY_TASK_ID -eq 0 && ! -f ${REF}.bwt ]]; then
  echo "$(date)  Indexing reference…"
  bwa index "$REF"
fi
wait   # guarantee index exists before others start

# ──────── align, sort, index BAM ──────────────────────────────────
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"
SAMTOOLS_THREADS=$(( SLURM_CPUS_PER_TASK - 1 ))

bwa mem -t "$SLURM_CPUS_PER_TASK" -R "$RG" "$REF" "$R1" "$R2" \
  | samtools sort -@ "$SAMTOOLS_THREADS" -o "${OUT_DIR}/${SAMPLE}.sorted.bam" -

samtools index "${OUT_DIR}/${SAMPLE}.sorted.bam"

echo "$(date)  ✅  finished $SAMPLE"
