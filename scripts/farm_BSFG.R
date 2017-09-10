#Purpose: run BSFG for the C. Hirsch maize data with the goal to generate gene modules and compare the modules to WGCNA results

# working directory

setwd("/home/caryn89/maize_BSFG")

# set daily working directory

rep = Sys.Date()
folder = sprintf('analysis/Rep_%s',rep)
try(dir.create(folder))
setwd(folder)

# Libraries

.libPaths(c("~/R/x86_64-pc-linux-gnu-library/3.3", .libPaths()))
library(BSFG)
#install.packages("MCMC")
#library(MCMC)

# --------------
# Data 

K <- read.table("/home/caryn89/Projects/maize_BSFG/data/processed/ZeaGBSv27_sF.sXX.txt")
Y <- read.table("/home/caryn89/Projects/maize_BSFG/data/processed/maize_414_expressedFPKM.bb", sep=",")
k_lines <- read.table("/home/caryn89/Projects/maize_BSFG/data/processed/ZeaGBSv27_503Exp_lines.txt")

# organize data
colnames(k_lines) <- "Line"
rownames(K) <- k_lines[,1]
rownames(Y) <- k_lines[,1]

# Make K into a Matrix using the Matrix package

library(Matrix)
K <- Matrix(as.matrix(K), sparse = TRUE)
K_mod = K
K_mod[K_mod < 1e-10] = 0 #modify the K matrix to remove extremely small values, and any negative values. make them zero
K_mod = forceSymmetric(drop0(K_mod,tol = 1e-10)) #force square matrix to being a symmetric Matrix, while keeping the Matrix sparse (no explicit zeros)
rownames(K_mod) = rownames(K)

#-------------------
# Initialize the priors

run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  #sampler = 'general_BSFG',
  simulation   = FALSE,
  scale_Y = TRUE, # internally rescales all the Y to have a mean value, z-score
  h2_divisions = 50,
  h2_step_size = NULL,
  burn = 1000
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


n_samples = 100;  # how many samples to collect at once?
for(i  in 1:70) {
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1)  # run MCMC chain n_samples iterations. grainSize is a parameter for parallelization (smaller = more parallelization)
  
  # set of commands to run during burn-in period to help chain converge
  if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn - 100) {
    BSFG_state = reorder_factors(BSFG_state) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    # BSFG_state$current_state = update_k(BSFG_state) # use to drop insignificant factors
    BSFG_state$run_parameters$burn = max(BSFG_state$run_parameters$burn,BSFG_state$current_state$nrun+100) # if you made changes, set a new burn-in period
    print(BSFG_state$run_parameters$burn)
  }
  BSFG_state = save_posterior_chunk(BSFG_state)  # save any accumulated posterior samples in the database to release memory
  print(BSFG_state) # print status of current chain
  plot(BSFG_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
}

save(BSFG_state, "postrun.RData")

