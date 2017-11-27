#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/maize_BSFG/analysis/rep-12-09-2017
#SBATCH -o /home/caryn89/Projects/maize_BSFG/logs/R_analysis_%j.out
#SBATCH -e /home/caryn89/Projects/maize_BSFG/logs/R_analysis_%j.out
#SBATCH -J analysis
#SBATCH -t 10:00:00
#SBATCH --mem=40000

set -e
set -o

Rscript analysis.R
echo $?


