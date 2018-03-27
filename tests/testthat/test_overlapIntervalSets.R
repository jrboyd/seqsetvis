library(seqsetvis)
library(GenomicRanges)
library(testthat)
gr_a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1))
gr_b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1))
gr_c = GRanges("chr1", IRanges(8:10*10+5, 8:10*10+6))

test_that("ssvOverlapIntervalSets grs input are valid. single GR fails.", {
  expect_error(ssvOverlapIntervalSets(gr_a))
})

test_that("ssvOverlapIntervalSets grs input are valid. single numeric fails.", {
    expect_error(ssvOverlapIntervalSets(1))
})

test_that("ssvOverlapIntervalSets grs input are valid. GRangesList OK.", {
    expect_s4_class(ssvOverlapIntervalSets(GRangesList("gr_a" = gr_a, "gr_b" = gr_b)), class = "GRanges")
})

test_that("ssvOverlapIntervalSets grs input are valid. List of GRanges OK.", {
    expect_s4_class(ssvOverlapIntervalSets(list("gr_a" = gr_a, "gr_b" = gr_b)), class = "GRanges")
})

test_that("ssvOverlapIntervalSets grs input are valid. GRangesList OK.", {
    expect_s4_class(ssvOverlapIntervalSets(GRangesList(list("gr_a" = gr_a, "gr_b" = gr_b))), class = "GRanges")
})

test_that("ssvOverlapIntervalSets grs input are valid. List of GRanges with non-matching mcols OK.", {
    gr_amod = gr_a
    gr_amod$score = "scored"
    expect_s4_class(ssvOverlapIntervalSets(list("gr_a" = gr_amod, "gr_b" = gr_b)), class = "GRanges")
})

test_that("ssvOverlapIntervalSets useFirst.", {
    ol = ssvOverlapIntervalSets(list("gr_a" = gr_a, "gr_b" = gr_b), use_first = T)
    expect_s4_class(ol, class = "GRanges")
    expect_equal(length(ol), length(gr_a))
})

test_that("ssvOverlapIntervalSets ext reduces number of results", {
    ol_ext1 = ssvOverlapIntervalSets(list("gr_a" = gr_a, "gr_b" = gr_b), ext = 1)
    ol_ext10 = ssvOverlapIntervalSets(list("gr_a" = gr_a, "gr_b" = gr_b), ext = 10)
    expect_equal(length(ol_ext1), 10)
    expect_equal(length(ol_ext10), 1)
})
