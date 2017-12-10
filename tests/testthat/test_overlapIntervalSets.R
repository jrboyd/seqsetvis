library(GenomicRanges)
library(testthat)
library(assertthat)
a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1))
b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1))
c = GRanges("chr1", IRanges(8:10*10+5, 8:10*10+6))

test_that("overlapIntervalSets grs input are valid", {
  expect_warning(overlapIntervalSets(list(a, b))) #unnamed list throws warning
  expect_error(overlapIntervalSets(a))
  expect_error(overlapIntervalSets(1))
  expect_s4_class(overlapIntervalSets(GRangesList("a" = a, "b" = b)), class = "GRanges")
  expect_s4_class(overlapIntervalSets(list("a" = a, "b" = b)), class = "GRanges")
})

test_that("gg_venn grs input are valid", {
  olap = overlapIntervalSets(list("a" = a, "b" = b, "c" = c))
  p = ggVenn(olap, circle.col = c("red", "blue", "green"))
  p
  expect_s3_class(p, class = "ggplot")


  olap = overlapIntervalSets(list("a" = a, "b" = b))
  p = ggVenn(elementMetadata(olap))
  p
  expect_s3_class(p, class = "ggplot")
  ggVenn(olap, circle.col = c("red", "blue", "green"), labels_size = 50, counts_size = 10, show_outside_count = F) +
    coord_fixed()
})
