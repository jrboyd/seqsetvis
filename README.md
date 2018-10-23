[![Coverage Status](https://img.shields.io/codecov/c/github/jrboyd/seqsetvis/master.svg)](https://codecov.io/github/jrboyd/seqsetvis?branch=master)
[![Travis-CI Build Status](https://travis-ci.org/jrboyd/seqsetvis.svg?branch=master)](https://travis-ci.org/jrboyd/seqsetvis)

# Installation

## From Bioconductor

```{r bioc install, eval=FALSE}
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("seqsetvis")
```

## From github
One dependency from bioconductor isn't getting installed automatically
```{r bioC missed dependency, eval=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomeInfoDbData")
```

Install devtools if needed
```{r devtools check, eval=FALSE}
if(!require(devtools)){
    install.packages("devtools")    
}
```

Install *seqsetvis* from github
```{r install seqsetvis from github, eval=FALSE}
devtools::install_github("jrboyd/seqsetvis")
#if you're still on R3.4
#devtools::install_github("jrboyd/seqsetvis", ref = "R3.4")
```
