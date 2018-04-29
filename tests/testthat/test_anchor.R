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

#sampling
test_that("viewGRangesWinSample_dt center", {
    dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, anchor = "center")
    # dt$group = "stranded"
    dt_us = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, anchor = "center_unstranded")
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
    dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, anchor = "left")
    # dt$group = "stranded"
    dt_us = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, anchor = "left_unstranded")
    # dt_us$group = "unstranded"
    # ssvSignalLineplot(rbind(dt, dt_us), sample_ = "id", color_ = "strand", group_ = "group")
    expect_equal(dt[strand == "+"][order(x)]$y, dt_us[strand == "+"][order(x)]$y)
    expect_failure(expect_equal(dt[strand == "-"][order(x)]$y, dt_us[strand == "-"][order(x)]$y))
    expect_true(all(dt$x > 0))
    expect_true(all(dt_us$x > 0))
})

#summary
test_that("viewGRangesWinSummary_dt center", {
    dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, anchor = "center")
    # dt$group = "stranded"
    dt_us = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, anchor = "center_unstranded")
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
    dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, anchor = "left")
    # dt$group = "stranded"
    dt_us = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, anchor = "left_unstranded")
    # dt_us$group = "unstranded"
    # ssvSignalLineplot(rbind(dt, dt_us), sample_ = "id", color_ = "strand", group_ = "group")
    expect_equal(dt[strand == "+"][order(x)]$y, dt_us[strand == "+"][order(x)]$y)
    expect_failure(expect_equal(dt[strand == "-"][order(x)]$y, dt_us[strand == "-"][order(x)]$y))
    expect_true(all(dt$x > 0))
    expect_true(all(dt_us$x > 0))
})

# bams
bams = c("A" = bam_file, "B" = bam_file)

test_that("ssvFetchBam anchor", {
    gr = ssvFetchBam(bams, qgr, win_size = 50, anchor = "center")
    gr_uns = ssvFetchBam(bams, qgr, win_size = 50, anchor = "center_unstranded")
    is_pos = as.character(strand(gr)) == "+"
    # plot(gr_uns$x[is_pos], gr$x[is_pos])
    # plot(gr_uns$x[!is_pos], gr$x[!is_pos])
    expect_equal(gr_uns$x[is_pos], gr$x[is_pos])
    expect_failure(expect_equal(gr_uns$x[!is_pos], gr$x[!is_pos]))
    expect_lt(min(gr$x), 0)
    expect_lt(min(gr_uns$x), 0)

    dt = ssvFetchBam(bams, qgr, win_size = 50, anchor = "left", return_data.table = TRUE)
    dt_uns = ssvFetchBam(bams, qgr, win_size = 50, anchor = "left_unstranded", return_data.table = TRUE)

    is_pos = dt$strand == "+"
    # plot(dt_uns$x[is_pos], dt$x[is_pos])
    # plot(dt_uns$x[!is_pos], dt$x[!is_pos])
    expect_equal(dt_uns$x[is_pos], dt$x[is_pos])
    expect_failure(expect_equal(dt_uns$x[!is_pos], dt$x[!is_pos]))
    expect_gte(min(dt$x), 0)
    expect_gte(min(dt_uns$x), 0)

})

#bigwigs
bigwigs = dir(system.file("extdata", package = "seqsetvis"), pattern = "random100.bw$", full.names = TRUE)
test_that("ssvFetchBigwig anchor", {
    gr = ssvFetchBigwig(bigwigs, qgr, win_size = 50, anchor = "center")
    gr_uns = ssvFetchBigwig(bigwigs, qgr, win_size = 50, anchor = "center_unstranded")
    is_pos = as.character(strand(gr)) == "+"
    # plot(gr_uns$x[is_pos], gr$x[is_pos])
    # plot(gr_uns$x[!is_pos], gr$x[!is_pos])
    expect_equal(gr_uns$x[is_pos], gr$x[is_pos])
    expect_failure(expect_equal(gr_uns$x[!is_pos], gr$x[!is_pos]))
    expect_lt(min(gr$x), 0)
    expect_lt(min(gr_uns$x), 0)

    dt = ssvFetchBigwig(bigwigs, qgr, win_size = 50, anchor = "left", return_data.table = TRUE)
    dt_uns = ssvFetchBigwig(bigwigs, qgr, win_size = 50, anchor = "left_unstranded", return_data.table = TRUE)

    is_pos = dt$strand == "+"
    # plot(dt_uns$x[is_pos], dt$x[is_pos])
    # plot(dt_uns$x[!is_pos], dt$x[!is_pos])
    expect_equal(dt_uns$x[is_pos], dt$x[is_pos])
    expect_failure(expect_equal(dt_uns$x[!is_pos], dt$x[!is_pos]))
    expect_gte(min(dt$x), 0)
    expect_gte(min(dt_uns$x), 0)
})


# bams
bams = c("A" = bam_file, "B" = bam_file)

test_that("ssvFetchBam anchor - summary", {
    gr = ssvFetchBam(bams, qgr, win_size = 50, anchor = "center", win_method = "summary")
    gr_uns = ssvFetchBam(bams, qgr, win_size = 50, anchor = "center_unstranded", win_method = "summary")
    is_pos = as.character(strand(gr)) == "+"
    # plot(gr_uns$x[is_pos], gr$x[is_pos])
    # plot(gr_uns$x[!is_pos], gr$x[!is_pos])
    expect_equal(gr_uns$x[is_pos], gr$x[is_pos])
    expect_failure(expect_equal(gr_uns$x[!is_pos], gr$x[!is_pos]))
    expect_lt(min(gr$x), 0)
    expect_lt(min(gr_uns$x), 0)

    dt = ssvFetchBam(bams, qgr, win_size = 50, anchor = "left", return_data.table = TRUE, win_method = "summary")
    dt_uns = ssvFetchBam(bams, qgr, win_size = 50, anchor = "left_unstranded", return_data.table = TRUE, win_method = "summary")

    is_pos = dt$strand == "+"
    # plot(dt_uns$x[is_pos], dt$x[is_pos])
    # plot(dt_uns$x[!is_pos], dt$x[!is_pos])
    expect_equal(dt_uns$x[is_pos], dt$x[is_pos])
    expect_failure(expect_equal(dt_uns$x[!is_pos], dt$x[!is_pos]))
    expect_gte(min(dt$x), 0)
    expect_gte(min(dt_uns$x), 0)

})

#bigwigs
bigwigs = dir(system.file("extdata", package = "seqsetvis"), pattern = "random100.bw$", full.names = TRUE)
test_that("ssvFetchBigwig anchor - summary", {
    gr = ssvFetchBigwig(bigwigs, qgr, win_size = 50, anchor = "center", win_method = "summary")
    gr_uns = ssvFetchBigwig(bigwigs, qgr, win_size = 50, anchor = "center_unstranded", win_method = "summary")
    is_pos = as.character(strand(gr)) == "+"
    # plot(gr_uns$x[is_pos], gr$x[is_pos])
    # plot(gr_uns$x[!is_pos], gr$x[!is_pos])
    expect_equal(gr_uns$x[is_pos], gr$x[is_pos])
    expect_failure(expect_equal(gr_uns$x[!is_pos], gr$x[!is_pos]))
    expect_lt(min(gr$x), 0)
    expect_lt(min(gr_uns$x), 0)

    dt = ssvFetchBigwig(bigwigs, qgr, win_size = 50, anchor = "left", return_data.table = TRUE, win_method = "summary")
    dt_uns = ssvFetchBigwig(bigwigs, qgr, win_size = 50, anchor = "left_unstranded", return_data.table = TRUE, win_method = "summary")

    is_pos = dt$strand == "+"
    # plot(dt_uns$x[is_pos], dt$x[is_pos])
    # plot(dt_uns$x[!is_pos], dt$x[!is_pos])
    expect_equal(dt_uns$x[is_pos], dt$x[is_pos])
    expect_failure(expect_equal(dt_uns$x[!is_pos], dt$x[!is_pos]))
    expect_gte(min(dt$x), 0)
    expect_gte(min(dt_uns$x), 0)
})

