# Running BSFG with maize RNA-Seq expression data

## Installing BSFG locally

BSFG: https://github.com/deruncie/SparseFactorMixedModel

```
library(git2r)
library(devtools)
devtools::install_github('deruncie/SparseFactorMixedModel',ref='develop',subdir='BSFG')
```

Also requires `MCMC` package.

## Installing BSFG on Farm

On Farm (UC Davis hpc), I (obviously) do not have permission to install R packages willy-nilly, and you probably don't either.
But, you do get to install whatever R packages you want into your local directories and then use those.
My home directory has an `R/` directory, where I have packages I wanted to install. 
The full path to those packages is: `/home/caryn89/R/x86_64-pc-linux-gnu-library/3.3`
As of this moment in time I don't remember how I set it up, but I will one day add that information.

So to install BSFG on Farm is similar to my local setup, but I have to change the library path to preferentially save to my personal `R/` directory:

```
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/3.3", .libPaths()))

library(git2r)
library(devtools)
devtools::install_github('deruncie/SparseFactorMixedModel',ref='develop',subdir='BSFG')
```

This requires testing still.

## Loading BSFG on Farm

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

## Updating BSFG

Dan pushes changes quite a bit, so it's a good idea to check with him to see what's he had changed before updating.

