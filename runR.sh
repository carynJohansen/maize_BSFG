#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/maize_BSFG
#SBATCH -o /home/caryn89/Projects/maize_BSFG/logs/R_out_%j.out
#SBATCH -e /home/caryn89/Projects/maize_BSFG/logs/R_out_%j.out
#SBATCH -J R-BSFG
#SBATCH -t 35:00:00
#SBATCH --mem=60000

set -e
set -o

Rscript scripts/farm_BSFG.R
echo $?


