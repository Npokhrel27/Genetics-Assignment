#!/bin/bash
#SBATCH -D /work/agro932/niranjan27/
#SBATCH -o /work/agro932/niranjan27/fst_output/fst-stdout-%j.txt
#SBATCH -e /work/agro932/niranjan27/fst_output/fst-stderr-%j.txt
#SBATCH -J fst_calc
#SBATCH -t 48:00:00
#SBATCH -p schnablelab
#SBATCH -c 16
#SBATCH --mem=130G
#SBATCH --mail-user=npokhrel3@huskers.unl.edu
#SBATCH --mail-type=END,FAIL

module load R

Rscript /work/agro932/niranjan27/scripts/calculate_fst.R
