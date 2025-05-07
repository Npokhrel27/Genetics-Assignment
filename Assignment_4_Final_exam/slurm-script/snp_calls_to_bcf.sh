#!/usr/bin/env bash
#SBATCH -D /work/agro932/niranjan27/Final_exam
#SBATCH -o /work/agro932/niranjan27/Final_exam/logs/snp-%A_%a.out
#SBATCH -e /work/agro932/niranjan27/Final_exam/logs/snp-%A_%a.err
#SBATCH -J snp_call_array                 # ← job name changed
#SBATCH -p schnablelab
#SBATCH --time=72:00:00
#SBATCH --array=0-9                       # 10 BAMs → indexes 0-9
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --mail-user=npokhrel3@huskers.unl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module load samtools bcftools

# ───── directories ────────────────────────────────────────────────
REF=/work/agro932/niranjan27/Reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa
BAMS_DIR=/work/agro932/niranjan27/Final_exam/aligned
OUTDIR=/work/agro932/niranjan27/snp_calls_whole_genome
LOGDIR=/work/agro932/niranjan27/Final_exam/logs

mkdir -p "$OUTDIR" "$LOGDIR"

# ───── BAM list & sample name ─────────────────────────────────────
BAMS=(${BAMS_DIR}/*.sorted.bam)
BAM=${BAMS[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename "$BAM" .sorted.bam)

# ───── index BAM if missing ───────────────────────────────────────
[[ -f ${BAM}.bai ]] || samtools index "$BAM"

# ───── SNP calling ────────────────────────────────────────────────
bcftools mpileup -Ou -f "$REF" -a AD,DP "$BAM" \
| bcftools call -mv -Ob -o "$OUTDIR/${SAMPLE}.bcf"

bcftools index "$OUTDIR/${SAMPLE}.bcf"

echo "$(date)  ✅  Finished SNP calling for ${SAMPLE}"
