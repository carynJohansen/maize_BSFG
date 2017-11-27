# Analysis of BSFG run that began on September 12, 2017

# Pathway for this analysis on farm, where is ran, is: /home/caryn89/Projects/maize_BSFG/analysis/rep-12-09-2017

setwd("/home/caryn89/Projects/maize_BSFG/analysis/rep-12-09-2017")

# Libraries
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/3.3", .libPaths()))
library(BSFG)
library(ggplot2)
library(reshape2)

# load data

load("current_loop.RData")

# reload the whole database of posterior samples

#This is too big for farm with a membory limit. Going to have to find a way around this.
#BSFG_state$Posterior <- reload_Posterior(BSFG_state)

#load Posterior lambda array
#Lambda <- BSFG_state$Posterior$Lambda
#dim(Lambda)
# examine the factor loadings across all the Posterior

#g1 <- Lambda[,1,] #the factor loadings for one gene across all Posteriors
#g1_m <- melt(g1)

#summary(g1_m)
#colnames(g1_m) <- c("posterior", "factor", "loading")

#pdf(file="gene1.pdf")
#ggplot(g1_m, aes(x=factor, y=loading, color=posterior)) + geom_point()
#dev.off()

# Visualize first Posterior

#f <- Lambda[1,,]
#f <- melt(f)
#colnames(f) <- c("gene", "factor", "loading")
#f$factor <- as.factor(f$factor)

#pdf("posterior1_lambda.pdf")
#ggplot(f, aes(x=factor, y=loading)) + geom_boxplot()
#dev.off()

# Visualize last Posterior

n <- dim(Lambda)[1] #the number of Posterior data available
p <- Lambda[n,,] #this will be the last posterior

p <- melt(p)
colnames(p) <- c("gene", "factor", "loading")
p$factor <- as.factor(p$factor)

pdf(sprintf("posterior%s_lambda.pdf", n))
ggplot(p, aes(x=factor, y=loading)) + geom_boxplot()
dev.off()

# modules

modules <- p[p$loading >= 0.2,]

table(modules$factor)

saveRDS(modules, file="modules.rds")
