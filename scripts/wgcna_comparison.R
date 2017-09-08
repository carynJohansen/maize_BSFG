
# ------------------ #
# Libraries
# ------------------ #

library(WGCNA)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(flashClust)

# ------------------ 
# Data

exp <- read.table("~/Box Sync/UCD/RRI/maize_genomics_project/data/eQTL/maize_414_expressedFPKM_names.txt", sep=",", header=T, row.names = 1)

lsdim(exp)
exp[1:5,1:5]

# explore the distributions of expression values for each sample
exp2 <- melt(as.matrix(exp))
head(exp2)
colnames(exp2) <- c("Line", "Gene", "FPKM")
str(exp2)

ggplot(exp2, aes(x=FPKM, color=Line)) + geom_density() + theme(legend.position="none")

#yeah, not normal

# ------------------ 
# Normalize


# convert to log scores

expL <- log10(exp)

#just to look
expLm <- melt(as.matrix(expL))
head(expLm)
colnames(expLm) <- c("Line", "Gene", "log")
ggplot(expLm, aes(x=log, color=Line)) + geom_density() + theme(legend.position="none")
# much more normal!

# normalize data using z-score

zscore <- function(row) {
	z_row <- (row - mean(row))/sd(row)
}

exp_z <- apply(expL, 1, zscore)
dim(exp_z)
exp_z[1:5, 1:5]
range(exp_z[,1])
exp_z_m <- melt(as.matrix(t(exp_z)))
colnames(exp_z_m) <- c("Line", "Gene", "z_score")

ggplot(exp_z_m, aes(x=z_score, color=Line)) + geom_density() + theme(legend.position="none") 

# apply transformed the direction of the expression data frame, so switch back.
expZ <- t(exp_z)
dim(expZ)
expZ[1:5,1:5]

# ------------------ 
# WGCNA


#start with subset of data, just to get this going. Come back and run all the data

d <- expZ

#begin WGCNA analysis 
#set powers to determine scale free network
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(d, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.95,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


softPower <- 4 #this is based on those plots.
adjacency_matrix <- adjacency(d, power = softPower)
system.time(TOM <- TOMsimilarity(adjacency_matrix))
dissTOM <- 1 - TOM

# Gene tree based on the dissimilarity matrix from TOM
geneTree <- flashClust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04)
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = 30)
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module Colors")

#or maybe...
#gene_network <- WGCNA::pickSoftThreshold(data=d)

#Using the one-step blockwise approach

system.time( net <- blockwiseModules(d, maxBlockSize = 1000, power = 4, 
                          TOMType = "signed", networkType = "signed",
                          minModuleSize = 30, reassignThreshold = 0, 
                          mergeCutHeight = 0.25, numericLabels = TRUE,
                          saveTOMs = TRUE, saveTOMFileBase ="cann_network_TOM",
                          nThreads = 3, verbose = 3))

str(net)
table(net$colors)
#mergedColors <- bwnet$colors
moduleColors <- labels2colors(net$colors)
#dynamicColors ?
table(moduleColors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels <- net$colors
molduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree = net$dendrograms[[1]];


# Merging the modules with similar expression profiles
MEs <- moduleEigengenes(d, moduleColors)$eigengenes
MEsO <- orderMEs(MEs)
dim(MEsO)
head(MEsO)
dim(MEs)
head(MEs)

#MEs <- net$MEs
#dim(MEs)
#head(MEs)

colorOrder <- c("grey", standardColors(50))
colorOrder
MEDiss = 1 - cor(MEs)
METree <- flashClust(as.dist(MEDiss), method = "average")
sizeGrWindow(7,6)
plot(METree, main = "Clustering of module Eigengenes")
#add threshold of 0.85 cutoff for merging modules
abline(h=0.15, col = "red")


merge <- mergeCloseModules(d, moduleColors, cutHeight = 0.15, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

sizeGrWindow(9, 5)
plotDendroAndColors(METree, cbind(moduleColors, mergedColors),
	c("Dynamic Tree Cut", "Merged dynamics"),
	dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

#merge some modules
MEList <- moduleEigengenes(d, colors = dynamicColors)
MEList
MEs <- MEList$eigengenes
MEDiss = 1 - cor(MEs)                    