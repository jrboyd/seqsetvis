library(seqsetvis)
library(GenomicRanges)
library(testthat)
gr_a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1))
gr_b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1))
gr_c = GRanges("chr1", IRanges(8:10*10+5, 8:10*10+6))

test_that("overlapIntervalSets grs input are valid. single GR fails.", {
  expect_error(overlapIntervalSets(gr_a))
})

test_that("overlapIntervalSets grs input are valid. single numeric fails.", {
    expect_error(overlapIntervalSets(1))
})

test_that("overlapIntervalSets grs input are valid. GRangesList OK.", {
    expect_s4_class(overlapIntervalSets(GRangesList("gr_a" = gr_a, "gr_b" = gr_b)), class = "GRanges")
})

test_that("overlapIntervalSets grs input are valid. List of Granges OK.", {
    expect_s4_class(overlapIntervalSets(list("gr_a" = gr_a, "gr_b" = gr_b)), class = "GRanges")
})

test_that("setPlotVenn various other paramters don't throw error", {
  olap = overlapIntervalSets(list("gr_a" = gr_a, "gr_b" = gr_b, "gr_c" = gr_c))
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
