testthat::context("FetchGranges")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
library(data.table)

qgr = CTCF_in_10a_overlaps_gr[1:10]
qgr = centerFixedSizeGRanges(qgr, 500)
#bed used to intersect bam
# rtracklayer::export.bed(qgr, con = "ctcf_5.bed")

# ssvFetchBam(bam_file, qgr)
res50 = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 50)
expect_equal(300, length(res50))

test_that("ssvFetchGRanges win_size", {
    res10 = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 10)
    expect_equal(300*5, length(res10))
})


test_that("ssvFetchGRanges attribs", {
    res10 = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 10, file_attribs = data.frame("a" = 1:3, "bc" = 4:6))
    expect_true(all(res10$a %in% 1:3))
    expect_true(all(res10$bc %in% 4:6))
})

test_that("ssvFetchGRanges win_method", {
    res5sum = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 5, win_method = "summary")
    expect_equal(length(qgr)*length(CTCF_in_10a_narrowPeak_grs)*5, length(res5sum))
    expect_gte(min(res5sum$x), -.5)
    expect_lte(max(res5sum$x), .5)
})


test_that("ssvFetchGRanges anchor", {
    res5left = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 5,
                              win_method = "summary", anchor = "left")
    expect_equal(length(qgr)*length(CTCF_in_10a_narrowPeak_grs)*5, length(res5left))
    expect_gte(min(res5left$x), 0)
    expect_lte(max(res5left$x), 1)
})

test_that("ssvFetchGRanges target_strand", {
    res50star = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr,
                              target_strand = "*", return_data.table = TRUE)
    expect_equal(max(res50star$y), 1)
    res50minus = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr,
                               target_strand = "-", return_data.table = TRUE)
    expect_equal(max(res50minus$y), 0)
    res50plus = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr,
                               target_strand = "+", return_data.table = TRUE)
    expect_equal(max(res50plus$y), 0)
    res50both = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr,
                               target_strand = "both", return_data.table = TRUE)
    expect_equal(max(res50both$y), 0)

    np_gr = CTCF_in_10a_narrowPeak_grs
    strand(np_gr$MCF10A_CTCF[1:3]) = "+"
    strand(np_gr$MCF10AT1_CTCF[4:6]) = "-"

    res50minus = ssvFetchGRanges(np_gr, qgr,
                                 target_strand = "-", return_data.table = TRUE)
    expect_equal(max(res50minus$y), 1)
    res50plus = ssvFetchGRanges(np_gr, qgr,
                                target_strand = "+", return_data.table = TRUE)
    expect_equal(max(res50plus$y), 1)
    res50both = ssvFetchGRanges(np_gr, qgr,
                                target_strand = "both", return_data.table = TRUE)
    expect_equal(max(res50both$y), 1)

    expect_equal(nrow(res50both), nrow(res50minus)*2)
})
