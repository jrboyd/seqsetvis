# library(peakvisr)
library(testthat)
a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1))
b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1))
c = GRanges("chr1", IRanges(8:10*10+5, 8:10*10+6))

test_that("overlapIntervalSets grs input are valid", {
  # expect_warning(overlapIntervalSets(list(a, b))) #unnamed list throws warning
  expect_error(overlapIntervalSets(a))
  expect_error(overlapIntervalSets(1))
  expect_s4_class(overlapIntervalSets(GRangesList("a" = a, "b" = b)), class = "GRanges")
  expect_s4_class(overlapIntervalSets(list("a" = a, "b" = b)), class = "GRanges")
})

test_that("setPlotVenn various other paramters don't throw error", {
  olap = overlapIntervalSets(list("a" = a, "b" = b, "c" = c))
  p = setPlotVenn(olap, circle_color = c("red", "blue", "green"), fill_alpha = .1,
             counts_txt_size = 10, show_outside_count = T,
             counts_as_labels = T)
  expect_s3_class(p, class = "ggplot")
  p = setPlotVenn(olap, circle_color = c("red", "blue", "green"), fill_circles = F,
             counts_txt_size = 10, show_outside_count = T,
             counts_as_labels = T)
  expect_s3_class(p, class = "ggplot")
})

test_that("col2hex", {
  string_colors = c("red", "blue", "green")
  rgb_colors = col2rgb(string_colors)
  hex_colors = c("#FF0000", "#0000FF", "#00FF00")
  expect_equal(col2hex(string_colors), hex_colors)
  expect_error(col2hex(rgb_colors))
  expect_error(col2hex("asdf"))
})
