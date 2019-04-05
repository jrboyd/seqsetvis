testthat::context("PrepareFetch")
# flipping viewGranges
library(seqsetvis)
library(GenomicRanges)
library(testthat)
library(data.table)

qgr = CTCF_in_10a_overlaps_gr[1:5]
start(qgr) = seq_len(5)
end(qgr) = 3*seq_len(5)
qgr = GenomicRanges::shift(qgr, -2)
strand(qgr) = c("+", "-", "-", "+", "-")

#sampling
test_that("prepare_fetch_GRanges start never less than 1", {
    expect_warning({
        tgr = prepare_fetch_GRanges(qgr, win_size = 5, target_size = 10)
    }, regexp = "Some out of bounds")
    expect_true(all(width(tgr) == 10))
    expect_true(all(start(tgr) >= 1))
})

