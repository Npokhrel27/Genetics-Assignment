#!/bin/bash
#SBATCH -D /work/agro932/niranjan27/fst_output/
#SBATCH -o /work/agro932/niranjan27/fst_output/fst_combine_stdout.txt
#SBATCH -e /work/agro932/niranjan27/fst_output/fst_combine_stderr.txt
#SBATCH -J combine_fst_plot
#SBATCH -t 01:00:00
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --mail-user=npokhrel3@huskers.unl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load R

Rscript combine_plot_fst.R
