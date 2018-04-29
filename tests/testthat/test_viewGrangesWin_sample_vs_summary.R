# flipping viewGranges
library(seqsetvis)
library(GenomicRanges)
library(testthat)

qgr = CTCF_in_10a_overlaps_gr[1:5]
strand(qgr) = c("+", "-", "-", "+", "-")
# qgr = centerFixedSizeGRanges(qgr, 500)
#bed used to intersect bam
# rtracklayer::export.bed(qgr, con = "ctcf_5.bed")
bam_file = system.file("extdata/test.bam", package = "seqsetvis")

bam_gr = fetchBam(bam_file, qgr)

#sampling
test_that("viewGRangesWinSample_dt center", {
    dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, x0 = "center")
    # dt$group = "stranded"
    dt_us = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, x0 = "center_unstranded")
    # dt_us$group = "unstranded"
    # ssvSignalLineplot(rbind(dt, dt_us), sample_ = "id", color_ = "strand", group_ = "group")
    expect_true(all(dt[strand == "+"][order(x)]$y == dt_us[strand == "+"][order(x)]$y))
    expect_false(all(dt[strand == "-"][order(x)]$y == dt_us[strand == "-"][order(x)]$y))
    expect_false(all(dt$x > 0))
    expect_false(all(dt$x < 0))
    expect_false(all(dt_us$x > 0))
    expect_false(all(dt_us$x < 0))
})

test_that("viewGRangesWinSample_dt left", {
    dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, x0 = "left")
    # dt$group = "stranded"
    dt_us = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, x0 = "left_unstranded")
    # dt_us$group = "unstranded"
    # ssvSignalLineplot(rbind(dt, dt_us), sample_ = "id", color_ = "strand", group_ = "group")
    expect_equal(dt[strand == "+"][order(x)]$y, dt_us[strand == "+"][order(x)]$y)
    expect_failure(expect_equal(dt[strand == "-"][order(x)]$y, dt_us[strand == "-"][order(x)]$y))
    expect_true(all(dt$x > 0))
    expect_true(all(dt_us$x > 0))
})

#sumamry
test_that("viewGRangesWinSummary_dt center", {
    dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, x0 = "center")
    # dt$group = "stranded"
    dt_us = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, x0 = "center_unstranded")
    # dt_us$group = "unstranded"
    # ssvSignalLineplot(rbind(dt, dt_us), sample_ = "id", color_ = "strand", group_ = "group")
    expect_equal(dt[strand == "+"][order(x)]$y, dt_us[strand == "+"][order(x)]$y)
    expect_failure(expect_equal(dt[strand == "-"][order(x)]$y, dt_us[strand == "-"][order(x)]$y))
    expect_false(all(dt$x > 0))
    expect_false(all(dt$x < 0))
    expect_false(all(dt_us$x > 0))
    expect_false(all(dt_us$x < 0))
})

test_that("viewGRangesWinSummary_dt left", {
    dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, x0 = "left")
    # dt$group = "stranded"
    dt_us = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, x0 = "left_unstranded")
    # dt_us$group = "unstranded"
    # ssvSignalLineplot(rbind(dt, dt_us), sample_ = "id", color_ = "strand", group_ = "group")
    expect_equal(dt[strand == "+"][order(x)]$y, dt_us[strand == "+"][order(x)]$y)
    expect_failure(expect_equal(dt[strand == "-"][order(x)]$y, dt_us[strand == "-"][order(x)]$y))
    expect_true(all(dt$x > 0))
    expect_true(all(dt_us$x > 0))
})

#test varying size regions
vgr = qgr
end(vgr) = end(vgr) + 0:4*100

test_that("viewGRangesWinSample_dt sizes vary, viewGRangesWinSummary_dt don't", {
    sample_dt = viewGRangesWinSample_dt(bam_gr, vgr, window_size = 100, x0 = "left")
    expect_gt(length(unique(sample_dt[, .N, by = id]$N)), 1)
    summary_dt = viewGRangesWinSummary_dt(bam_gr, vgr, n_tiles = 10, x0 = "left")
    summary_dt[, .N, by = id]
    expect_equal(length(unique(summary_dt[, .N, by = id]$N)), 1)
})


test_that("can fetch bam summary", {
    #sample is bp scale
    sample_dt = fetchWindowedBam(bam_file, vgr, return_data.table = TRUE)
    expect_gte(max(abs(range(sample_dt$x))), 100)
    #summaries have range of x < .5 for center x0
    summary_dt = fetchWindowedBam(bam_file, vgr, win_method = "summary", return_data.table = TRUE)
    expect_lte(max(abs(range(summary_dt$x))), .5)
    summary_dt.med = fetchWindowedBam(bam_file, vgr[5], win_method = "summary",
                                      summary_FUN = function(x, w){min(w)}, return_data.table = TRUE)
    expect_lte(max(abs(range(summary_dt.med$x))), .5)
    #summary fun did something
    expect_false(all(summary_dt$y == summary_dt.med$y))
    #win size did something
    summary_dt10 = fetchWindowedBam(bam_file, vgr, win_method = "summary",
                                  win_size = 10,
                                  return_data.table = TRUE)
    expect_true(all(summary_dt10[, .N, by = id]$N == 10))
    expect_true(all(summary_dt[, .N, by = id]$N == 50))
})

test_that("can fetch bamList summary", {
    bams = c("A" = bam_file, "B" = bam_file)
    bams_sample_dt = fetchWindowedBamList(bams, vgr)
    expect_gte(max(abs(range(bams_sample_dt$x))), 100)
    bams_summary_dt = fetchWindowedBamList(bams, vgr, win_method = "summary")
    expect_lte(max(abs(range(bams_summary_dt$x))), .5)
    bams_summary_dt10 = fetchWindowedBamList(bams, vgr, win_method = "summary",
                                             win_size = 10)
    bams_summary_dtFUN = fetchWindowedBamList(bams, vgr, win_method = "summary",
                                              summary_FUN = function(x, w){min(w)})
    expect_false(all(bams_summary_dt$y == bams_summary_dtFUN$y))
    expect_true(all(bams_summary_dt10[, .N, by = id]$N == 10*length(bams)))
    expect_true(all(bams_summary_dt[, .N, by = id]$N == 50*length(bams)))
})

bigwig_file = system.file("extdata/MCF10A_CTCF_FE_random100.bw", package = "seqsetvis", mustWork = TRUE)

test_that("can fetch bigwig summary", {
    #sample is bp scale
    sample_dt = fetchWindowedBigwig(bigwig_file, vgr, return_data.table = TRUE)
    expect_gte(max(abs(range(sample_dt$x))), 100)
    #summaries have range of x < .5 for center x0
    summary_dt = fetchWindowedBigwig(bigwig_file, vgr, win_method = "summary", return_data.table = TRUE)
    expect_lte(max(abs(range(summary_dt$x))), .5)
    summary_dt.med = fetchWindowedBigwig(bigwig_file, vgr[5], win_method = "summary",
                                      summary_FUN = function(x, w){min(w)}, return_data.table = TRUE)
    expect_lte(max(abs(range(summary_dt.med$x))), .5)
    #summary fun did something
    expect_false(all(summary_dt$y == summary_dt.med$y))
    #win size did something
    summary_dt10 = fetchWindowedBigwig(bigwig_file, vgr, win_method = "summary",
                                    win_size = 10,
                                    return_data.table = TRUE)
    expect_true(all(summary_dt10[, .N, by = id]$N == 10))
    expect_true(all(summary_dt[, .N, by = id]$N == 50))
})

test_that("can fetch bigwigList summary", {
    bigwigs = c("A" = bigwig_file, "B" = bigwig_file)
    bigwigs_sample_dt = fetchWindowedBigwigList(bigwigs, vgr)
    expect_gte(max(abs(range(bigwigs_sample_dt$x))), 100)
    bigwigs_summary_dt = fetchWindowedBigwigList(bigwigs, vgr, win_method = "summary")
    expect_lte(max(abs(range(bigwigs_summary_dt$x))), .5)
    bigwigs_summary_dt10 = fetchWindowedBigwigList(bigwigs, vgr, win_method = "summary",
                                             win_size = 10)
    bigwigs_summary_dtFUN = fetchWindowedBigwigList(bigwigs, vgr, win_method = "summary",
                                              summary_FUN = function(x, w){min(w)})
    expect_false(all(bigwigs_summary_dt$y == bigwigs_summary_dtFUN$y))
    expect_true(all(bigwigs_summary_dt10[, .N, by = id]$N == 10*length(bigwigs)))
    expect_true(all(bigwigs_summary_dt[, .N, by = id]$N == 50*length(bigwigs)))
})

# fetchBam and fetchBw params
# win_method
# summary_FUN
