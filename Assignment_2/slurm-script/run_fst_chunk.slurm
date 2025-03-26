#!/bin/bash
#SBATCH -D /work/agro932/niranjan27/
#SBATCH -J fst_chunk
#SBATCH -o /work/agro932/niranjan27/fst_output/fst-stdout-%A_%a.txt
#SBATCH -e /work/agro932/niranjan27/fst_output/fst-stderr-%A_%a.txt
#SBATCH -t 24:00:00
#SBATCH -p schnablelab
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --array=1-90
#SBATCH --mail-type=FAIL

# === Load R ===
module load R

# === Define chunk filenames ===
CHUNK_DIR="/work/agro932/niranjan27/snp_calls_whole_genome/chunks"
CHUNK_LIST=($(ls $CHUNK_DIR/snp_chunk_*.txt))
INPUT_FILE=${CHUNK_LIST[$SLURM_ARRAY_TASK_ID-1]}

# === Output filenames ===
OUT_DIR="/work/agro932/niranjan27/fst_output"
BASENAME=$(basename "$INPUT_FILE" .txt)
CSV_OUT="${OUT_DIR}/${BASENAME}_fst.csv"
PLOT_OUT="${OUT_DIR}/${BASENAME}_fst.png"

# === Run R script ===
Rscript /work/agro932/niranjan27/scripts/calculate_fst_chunk.R "$INPUT_FILE" "$CSV_OUT" "$PLOT_OUT"
