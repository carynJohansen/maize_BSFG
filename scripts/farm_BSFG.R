#Purpose: run BSFG for the C. Hirsch maize data with the goal to generate gene modules and compare the modules to WGCNA results

# working directory

setwd("/home/caryn89/maize_BSFG")

# Libraries

.libPaths(c("~/R/x86_64-pc-linux-gnu-library/3.3", .libPaths()))
library(BSFG)
library(MCMC)

# --------------
# Data 

K <- read.table("/home/caryn89/rotations/maize_genomics/data/BSFG_input/ZeaGBSv27_sF.sXX.txt")
Y <- read.table("/home/caryn89/rotations/maize_genomics/data/BSFG_input/maize_414_expressedFPKM.bb", sep=",")
k_lines <- ("/home/caryn89/rotations/maize_genomics/data/BSFG_input/ZeaGBSv27_414Exp_lines.txt")

# organize data
colnames(k_lines) <- "Line"
rownames(K) <- k_lines[,1]
rownames(Y) <- k_lines[,1]



