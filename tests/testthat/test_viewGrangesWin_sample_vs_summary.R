# flipping viewGranges
library(seqsetvis)
library(GenomicRanges)
library(testthat)
library(data.table)

qgr = CTCF_in_10a_overlaps_gr[1:5]
strand(qgr) = c("+", "-", "-", "+", "-")
# qgr = centerFixedSizeGRanges(qgr, 500)
#bed used to intersect bam
# rtracklayer::export.bed(qgr, con = "ctcf_5.bed")
bam_file = system.file("extdata/test.bam", package = "seqsetvis")

bam_gr = fetchBam(bam_file, qgr)

#test varying size regions
vgr = qgr
end(vgr) = end(vgr) + 0:4*100

test_that("viewGrangeWinSample_dt ids match input", {
    names(qgr) = paste0("peak_", seq_along(qgr))
    sample_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 100)
    expect_equal(sort(unique(names(qgr))), sort(unique(sample_dt$id)))
})

test_that("viewGrangeWinSummary_dt ids match input", {
    names(qgr) = paste0("peak_", seq_along(qgr))
    summary_dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 100)
    expect_equal(sort(unique(names(qgr))), sort(unique(summary_dt$id)))
})

test_that("viewGrangeWinSample_dt unnamed qgr still creates id", {
    names(qgr) = NULL
    sample_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 100)
    expect_true(!is.null(sample_dt$id))
})

test_that("viewGrangeWinSummary_dt unnamed qgr still creates id", {
    names(qgr) = NULL
    summary_dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 100)
    expect_true(!is.null(sample_dt$id))
})

test_that("viewGrangeWinSample_dt ids match input", {
    names(qgr) = paste0("peak_", seq_along(qgr))
    summary_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 100)
    expect_equal(sort(unique(names(qgr))), sort(unique(summary_dt$id)))

    names(qgr) = seq_along(qgr)
    summary_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 100)
    expect_equal(sort(unique(names(qgr))), sort(unique(summary_dt$id)))
})

test_that("viewGrangeWinSummary_dt ids match input", {
    names(qgr) = paste0("peak_", seq_along(qgr))
    summary_dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 100)
    expect_equal(sort(unique(names(qgr))), sort(unique(summary_dt$id)))

    names(qgr) = seq_along(qgr)
    summary_dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 100)
    expect_equal(sort(unique(names(qgr))), sort(unique(summary_dt$id)))
})


test_that("viewGRangesWinSample_dt sizes vary, viewGRangesWinSummary_dt don't", {
    sample_dt = viewGRangesWinSample_dt(bam_gr, vgr, window_size = 100, anchor = "left")
    expect_gt(length(unique(sample_dt[, .N, by = id]$N)), 1)
    summary_dt = viewGRangesWinSummary_dt(bam_gr, vgr, n_tiles = 10, anchor = "left")
    summary_dt[, .N, by = id]
    expect_equal(length(unique(summary_dt[, .N, by = id]$N)), 1)
})


test_that("can fetch bam summary", {
    #sample is bp scale
    sample_dt = ssvFetchBam(bam_file, vgr, return_data.table = TRUE)
    expect_gte(max(abs(range(sample_dt$x))), 100)
    #summaries have range of x < .5 for center anchor
    summary_dt = ssvFetchBam(bam_file, vgr, win_method = "summary", return_data.table = TRUE)
    expect_lte(max(abs(range(summary_dt$x))), .5)
    summary_dt.med = ssvFetchBam(bam_file, vgr[5], win_method = "summary",
                                 summary_FUN = function(x, w){min(w)}, return_data.table = TRUE)
    expect_lte(max(abs(range(summary_dt.med$x))), .5)
    #summary fun did something
    expect_false(all(summary_dt$y == summary_dt.med$y))
    #win size did something
    summary_dt10 = ssvFetchBam(bam_file, vgr, win_method = "summary",
                               win_size = 10,
                               return_data.table = TRUE)
    expect_true(all(summary_dt10[, .N, by = id]$N == 10))
    expect_true(all(summary_dt[, .N, by = id]$N == 50))
})

test_that("can fetch bam summary", {
    bams = c("A" = bam_file, "B" = bam_file)
    bams_sample_dt = ssvFetchBam(bams, vgr)
    expect_gte(max(abs(range(bams_sample_dt$x))), 100)
    bams_summary_dt = ssvFetchBam(bams, vgr, win_method = "summary")
    expect_lte(max(abs(range(bams_summary_dt$x))), .5)
    bams_summary_dt10 = ssvFetchBam(bams, vgr, win_method = "summary",
                                    win_size = 10, return_data.table = TRUE)
    bams_summary_grFUN = ssvFetchBam(bams, vgr, win_method = "summary",
                                     summary_FUN = function(x, w){min(w)})
    expect_false(all(bams_summary_dt$y == bams_summary_grFUN$y))
    expect_true(all(bams_summary_dt10[, .N, by = id]$N == 10*length(bams)))
    expect_true(all(as.data.table(bams_summary_dt)[, .N, by = id]$N == 50*length(bams)))
})

bigwig_file = system.file("extdata/MCF10A_CTCF_FE_random100.bw", package = "seqsetvis", mustWork = TRUE)

test_that("can fetch bigwig summary", {
    skip_on_os("windows")
    #sample is bp scale
    sample_dt = ssvFetchBigwig(bigwig_file, vgr, return_data.table = TRUE)
    expect_gte(max(abs(range(sample_dt$x))), 100)
    #summaries have range of x < .5 for center anchor
    summary_dt = ssvFetchBigwig(bigwig_file, vgr, win_method = "summary", return_data.table = TRUE)
    expect_lte(max(abs(range(summary_dt$x))), .5)
    summary_dt.med = ssvFetchBigwig(bigwig_file, vgr[5], win_method = "summary",
                                    summary_FUN = function(x, w){min(w)}, return_data.table = TRUE)
    expect_lte(max(abs(range(summary_dt.med$x))), .5)
    #summary fun did something
    expect_false(all(summary_dt$y == summary_dt.med$y))
    #win size did something
    summary_dt10 = ssvFetchBigwig(bigwig_file, vgr, win_method = "summary",
                                  win_size = 10,
                                  return_data.table = TRUE)
    expect_true(all(summary_dt10[, .N, by = id]$N == 10))
    expect_true(all(summary_dt[, .N, by = id]$N == 50))
})

test_that("can fetch bigwig summary", {
    skip_on_os("windows")
    bigwigs = c("A" = bigwig_file, "B" = bigwig_file)
    bigwigs_sample_dt = ssvFetchBigwig(bigwigs, vgr)
    expect_gte(max(abs(range(bigwigs_sample_dt$x))), 100)
    bigwigs_summary_dt = ssvFetchBigwig(bigwigs, vgr, win_method = "summary")
    expect_lte(max(abs(range(bigwigs_summary_dt$x))), .5)
    bigwigs_summary_dt10 = ssvFetchBigwig(bigwigs, vgr, win_method = "summary",
                                          win_size = 10, return_data.table = TRUE)
    bigwigs_summary_dtFUN = ssvFetchBigwig(bigwigs, vgr, win_method = "summary",
                                           summary_FUN = function(x, w){min(w)})
    expect_false(all(bigwigs_summary_dt$y == bigwigs_summary_dtFUN$y))
    expect_true(all(bigwigs_summary_dt10[, .N, by = id]$N == 10*length(bigwigs)))
    expect_true(all(as.data.table(bigwigs_summary_dt)[, .N, by = id]$N == 50*length(bigwigs)))
})

# fetchBam and fetchBw params
# win_method
# summary_FUN
