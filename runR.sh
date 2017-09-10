#!/bin/bash -l

#SBATCH -D /home/caryn89/maize_BSFG
#SBATCH -o /home/caryn89/maize_BSFG/logs/R_out_%j.out
#SBATCH -e /home/caryn89/maize_BSFG/logs/R_out_%j.out
#SBATCH -J R-BSFG
#SBATCH -t 24:00:00

set -e
set -o

Rscript scripts/farm_BSFG.R
echo $?


