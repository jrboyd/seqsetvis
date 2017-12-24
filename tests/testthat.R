
## load dependencies
library(testthat)
library(seqsetvis)
library(data.table)

## test package
test_check(package = "seqsetvis", path = "seqsetvis.Rcheck/00_pkg_src/seqsetvis/tests/testthat/", filter = "test_")
