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
expect_equal(300, nrow(res50))

test_that("ssvFetchGRanges win_size", {
    res10 = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 10)
    expect_equal(300*5, nrow(res10))
})

test_that("ssvFetchGRanges win_method", {
    res5sum = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 5, win_method = "summary")
    expect_equal(length(qgr)*length(CTCF_in_10a_narrowPeak_grs)*5, nrow(res5sum))
    expect_gte(min(res5sum$y), -.5)
    expect_lte(max(res5sum$y), .5)
})


test_that("ssvFetchGRanges anchor", {
    res5sum = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 5,
                              win_method = "summary", anchor = "left")
    expect_equal(length(qgr)*length(CTCF_in_10a_narrowPeak_grs)*5, nrow(res5sum))
    expect_gte(min(res5sum$y), 0)
    expect_lte(max(res5sum$y), 1)
})
