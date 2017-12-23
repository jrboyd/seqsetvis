pkgname <- "seqsetvis"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('seqsetvis')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("overlapIntervalSets")
### * overlapIntervalSets

flush(stderr()); flush(stdout())

### Name: overlapIntervalSets
### Title: Intersect a list of GRanges to create a single GRanges object of
###   merged ranges including metadata describing overlaps per input
###   GRanges
### Aliases: overlapIntervalSets

### ** Examples

library(GenomicRanges)
a = GRanges("chr1", IRanges(1:7*10, 1:7*10))
b = GRanges("chr1", IRanges(5:10*10, 5:10*10))
overlapIntervalSets(list(a, b))



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
