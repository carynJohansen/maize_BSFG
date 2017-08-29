
#-------------------------
# Installing BSFG

#library(git2r)
#library(devtools)
#devtools::install_github('deruncie/SparseFactorMixedModel',ref='develop',subdir='BSFG')
library(BSFG)
library(MCMCpack)

#----------------
# set working directory
setwd("~/Box Sync/Projects/BSFG/maize/")

#----------------
# Load Data

K <- read.table("~/Box Sync/UCD/RRI/maize_genomics_project/data/Exp_RelMat/ZeaGBSv27_sF.sXX.txt")
Y <- read.table("~/Box Sync/UCD/RRI/maize_genomics_project/data/eQTL/maize_414_expressedFPKM.bb", sep=",")
k_lines <- read.table("~/Box Sync/UCD/RRI/maize_genomics_project/data/Exp_RelMat/ZeaGBSv27_503Exp_lines.txt")
colnames(k_lines) <- "Line"
rownames(K) <- k_lines[,1]
rownames(Y) <- k_lines[,1]


# ----------------------
# Make K into a Matrix using the Matrix package


library(Matrix)
K <- Matrix(as.matrix(K), sparse = TRUE)
K_mod = K
K_mod[K_mod < 1e-10] = 0 #modify the K matrix to remove extremely small values, and any negative values. make them zero
K_mod = forceSymmetric(drop0(K_mod,tol = 1e-10)) #force square matrix to being a symmetric Matrix, while keeping the Matrix sparse (no explicit zeros)
rownames(K_mod) = rownames(K)

#------------------------
# set daily working directory

rep = Sys.Date()
folder = sprintf('Rep_%s',rep)
try(dir.create(folder))
setwd(folder)

#-------------------
# Initialize the priors

run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  #sampler = 'general_BSFG',
  simulation   = FALSE,
  scale_Y = TRUE, # internally rescales all the Y to have a mean value, z-score
  h2_divisions = 2,
  h2_step_size = NULL,
  burn = 100
)

#--------------------
#Set the prior hyperparameters of the BSFG model

priors = BSFG_priors(
  fixed_var = list(V = 1,     nu = 3),
  # tot_Y_var = list(V = 0.5,   nu = 3),
  tot_Y_var = list(V = 0.5,   nu = 5),
  tot_F_var = list(V = 18/20, nu = 20),
  delta_1   = list(shape = 2.1,  rate = 1/20),
  delta_2   = list(shape = 3, rate = 1),
  Lambda_df = 3,
  B_df      = 3,
  B_F_df    = 3,
  h2_priors_resids_fun = function(h2s,n) 1,#pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1#ifelse(h2s == 0,n,n/(n-1))
)

#----------------------
# Construct model

#?BSFG_init
BSFG_state = BSFG_init(Y,
                       model=~1 + (1|Line), # This model has one fixed and one random term.
                       factor_model_fixed = NULL,
                       priors=priors,
                       run_parameters=run_parameters, 
                       #factor_model_fixed = ~0, # we could specify a different fixed effect model for the factors
                       data = k_lines, # the data.frame with information for constructing the model matrices
                       K_mats = list(Line = K_mod) # covariance matrices for the random effects. If not provided, assume uncorrelated
)

#------------------
# MCMC


n_samples = 1;  # how many samples to collect at once?
for(i  in 1:70) {
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1)  # run MCMC chain n_samples iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
  
  # set of commands to run during burn-in period to help chain converge
  #if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn - 100) {
  #  BSFG_state = reorder_factors(BSFG_state) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
  #  # BSFG_state$current_state = update_k(BSFG_state) # use to drop insignificant factors
  #  BSFG_state$run_parameters$burn = max(BSFG_state$run_parameters$burn,BSFG_state$current_state$nrun+100) # if you made changes, set a new burn-in period
  #  print(BSFG_state$run_parameters$burn)
  #}
  BSFG_state = save_posterior_chunk(BSFG_state)  # save any accumulated posterior samples in the database to release memory
  print(BSFG_state) # print status of current chain
  plot(BSFG_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
}

#------------------------
# Working with the Posterior

# reload the whole database of posterior samples
BSFG_state$Posterior = reload_Posterior(BSFG_state)

# What is the difference between this BSFG_state object, and if I were to load(current_state.RData)?
# How do I pick up a project from where I left off? It seems like I can't just load the RData file and work with the posterior.

# all parameter names in Posterior
BSFG_state$Posterior$posteriorSample_params
BSFG_state$Posterior$posteriorMean_params  # these ones only have the posterior mean saved, not individual posterior samples

# instead, load only a specific parameter
Lambda = load_posterior_param(BSFG_state,'Lambda')
image(as.matrix(Lambda))

library(ggplot2)
library(reshape2)
g1 <- Lambda[,1,] # the factor loadings for the samples in the posterior
g1.m <- melt(g1)
summary(g1.m)
colnames(g1.m) <- c("sample", "factor", "loading")

ggplot(g1.m, aes(x=factor, y=loading, color=sample)) + geom_point()
ggplot(g1.m, aes(x=sample, y=loading, color=as.factor(factor))) + geom_line()

# boxplots are good ways to visualize Posterior distributions on sets of related parameters
boxplot(BSFG_state$Posterior$F_h2[,1,])

# get posterior distribution on a function of parameters
G_samples = get_posterior_FUN(BSFG_state,Lambda %*% diag(F_h2[1,]) %*% t(Lambda) + resid_h2[1,]/tot_Eta_prec[1,])

# get posterior mean of a parameter
G = get_posterior_mean(G_samples)

# get Highest Posterior Density intervals for paramters
F_h2_HPD = get_posterior_HPDinterval(BSFG_state,F_h2)




