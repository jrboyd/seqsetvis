[![Coverage Status](https://img.shields.io/codecov/c/github/jrboyd/seqsetvis/master.svg)](https://codecov.io/github/jrboyd/seqsetvis?branch=master)
[![Travis-CI Build Status](https://travis-ci.org/jrboyd/seqsetvis.svg?branch=master)](https://travis-ci.org/jrboyd/seqsetvis)

# Installation

## From Bioconductor

```{r bioc install, eval=FALSE}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("seqsetvis")
```

## From github
One dependency from bioconductor isn't getting installed automatically
```{r bioC missed dependency, eval=FALSE}
source("https://bioconductor.org/biocLite.R")
biocLite("GenomeInfoDbData")
```

Install devtools if needed
```{r devtools check, eval=FALSE}
if(!require(devtools)){
    install.packages("devtools")    
}
```

Install *seqsetvis* from github
```{r install seqsetvis from github, eval=FALSE}
devtools::install_github("jrboyd/seqsetvis", build_vignettes = TRUE)
```
