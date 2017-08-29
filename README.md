# Running BSFG with maize RNA-Seq expression data

## Installing BSFG locally

BSFG: https://github.com/deruncie/SparseFactorMixedModel

```
library(git2r)
library(devtools)
devtools::install_github('deruncie/SparseFactorMixedModel',ref='develop',subdir='BSFG')
```

Also requires `MCMC` package

## Running BSFG on Farm

You must modify your library paths to run BSFG from a local repository:


Option 1: add the local path to .libPaths()
```
.libPaths(c(.libPaths(),"~/R/x86_64-pc-linux-gnu-library/3.3"))
library(BSFG)
```

Option 2: (have not tried)
```
library("BSFG", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
```

The packages `git2r` and `devtools` appear to be installed.

This could perhaps be done on a more permanent basis by altering an .Rprofile file, but I have not yet looked into it.

