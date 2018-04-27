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
    expect_equal(dt[strand == "+"], dt_us[strand == "+"])
    expect_failure(expect_equal(dt[strand == "-"], dt_us[strand == "-"]))
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
    expect_equal(dt[strand == "+"], dt_us[strand == "+"])
    expect_failure(expect_equal(dt[strand == "-"], dt_us[strand == "-"]))
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
    expect_equal(dt[strand == "+"], dt_us[strand == "+"])
    expect_failure(expect_equal(dt[strand == "-"], dt_us[strand == "-"]))
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
    expect_equal(dt[strand == "+"], dt_us[strand == "+"])
    expect_failure(expect_equal(dt[strand == "-"], dt_us[strand == "-"]))
    expect_true(all(dt$x > 0))
    expect_true(all(dt_us$x > 0))
})

#test varying size regions
vgr = qgr
end(vgr) = end(vgr) + 0:4*100

test_that("viewGrangesWin* sizes vary", {
    viewGRangesWin(bam_gr, qgr, n_tiles = 10, x0 = "left")
    viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, x0 = "left")
})

# fetchBam and fetchBw params
# win_method
# summary_FUN
