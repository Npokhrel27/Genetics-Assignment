#!/usr/bin/env bash
# tidy_final_exam.sh
# Put FASTQ in fastq/, SRAlite in sralite/, everything else in logs/

set -euo pipefail
shopt -s nullglob          # “mv” silently skips if nothing matches

BASE="/work/agro932/niranjan27/Final_exam"
cd "$BASE"

# ── create target dirs only if needed ────────────────────────────────
mkdir -p fastq sralite logs

# ── move FASTQ (plain or gzipped) ────────────────────────────────────
mv *.fastq *.fastq.gz *.fq *.fq.gz fastq/ 2>/dev/null || true

# ── move SRAlite archives ────────────────────────────────────────────
mv *.sralite* sralite/ 2>/dev/null || true

# ── move logs, SLURM outputs, scripts, lists, etc. ───────────────────
mv *.err *.out *.log *.txt *.slurm logs/ 2>/dev/null || true

echo "✅  FASTQ ⇒ fastq/  |  SRAlite ⇒ sralite/  |  logs & misc ⇒ logs/"
