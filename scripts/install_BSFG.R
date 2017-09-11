.libPaths(c("~/R/x86_64-pc-linux-gnu-library/3.3", .libPaths()))

library(git2r)
library(devtools)
#devtools::install_github('deruncie/SparseFactorMixedModel',ref='develop',subdir='BSFG')
BSFG_path <- 'https://github.com/deruncie/SparseFactorMixedModel'
my_local_library <- '~/R/x86_64-pc-linux-gnu-library/3.3'
withr::with_libpaths(my_local_library,install_git(BSFG_path,branch = 'develop',subdir = 'BSFG'))


