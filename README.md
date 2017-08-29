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
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/3.3", .libPaths()))
library(BSFG)
```

The string `"~/R/x86_64-pc-linux-gnu-library/3.3"` is my path from my home directory to my locally installed R packages. It goes before `.libPaths()` so R looks there first.
The benefit to this is that it's version controlled - at least, because I am lazy and won't update the local until something breaks. (Note to self: be more diligent about versions of software used in R code. Explore the session.info())

Option 2: (have not tried)
```
library("BSFG", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
```

The packages `git2r` and `devtools` appear to be installed.

This could perhaps be done on a more permanent basis by altering an .Rprofile file, but I have not yet looked into it.

