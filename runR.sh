#!/bin/bash -l

#SBATCH -D /home/caryn89/Projects/maize_BSFG
#SBATCH -o /home/caryn89/Projects/maize_BSFG/logs/bsfg_out_%j.out
#SBATCH -e /home/caryn89/Projects/maize_BSFG/logs/bsfg_out_%j.out
#SBATCH -J R-BSFG_4
#SBATCH -t 35:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=40000

set -e
set -o

Rscript scripts/farm_BSFG.R
echo $?


