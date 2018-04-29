library(seqsetvis)
library(testthat)

qgr = CTCF_in_10a_overlaps_gr[1:5]
qgr = centerFixedSizeGRanges(qgr, 500)
#bed used to intersect bam
# rtracklayer::export.bed(qgr, con = "ctcf_5.bed")
bam_file = system.file("extdata/test.bam", package = "seqsetvis")
xls_file = system.file("extdata/test_peaks.xls", package = "seqsetvis")

# ssvFetchBam(bam_file, qgr)

test_that("fragLen_fromMacs2Xls", {
    fl = fragLen_fromMacs2Xls(xls_file)
    expect_equal(fl, 187)
})

test_that("fetchBam auto fragLen", {
    expect_message(fetchBam(bam_file, qgr = qgr), "fragLen was calculated as:")
})

test_that("fragLen_calcStranded parameters matter", {
    res2 = fragLen_calcStranded(bam_file, qgr = qgr[1:2], n_regions = 2, include_plot_in_output = TRUE)
    res1 = fragLen_calcStranded(bam_file, qgr = qgr[1:2], n_regions = 1, include_plot_in_output = TRUE)

    f1 = res2[[1]]
    f2 = res1[[1]]
    expect_failure(expect_equal(f1, f2))

    f3 = fragLen_calcStranded(bam_file, qgr = qgr[1], n_regions = 1, include_plot_in_output = TRUE, ma_distance = 1)[[1]]
    f4 = fragLen_calcStranded(bam_file, qgr = qgr[1], n_regions = 1, include_plot_in_output = TRUE, ma_distance = 33)[[1]]
    expect_failure(expect_equal(f3, f4))

    res5 = fragLen_calcStranded(bam_file, qgr = qgr, max_fragLen = 500, include_plot_in_output = TRUE)
    res6 = fragLen_calcStranded(bam_file, qgr = qgr, max_fragLen = 100, include_plot_in_output = TRUE)

    f5 = res5[[1]]
    f6 = res6[[1]]
    expect_failure(expect_equal(f5, f6))
})

test_that("fragLen_calcStranded can force no which", {
    expect_error({fragLen_calcStranded(bam_file, qgr = NULL)},
                 regexp = "This will probably be very slow and uneccessary")
    expect_failure(expect_error({fragLen_calcStranded(bam_file, qgr = NULL, force_no_which = TRUE)}))
})

test_that("viewGRangesWinSample_dt strand and position functions", {


    bam_score = fetchBam(bam_file, qgr = qgr)
    score_gr = bam_score
    window_size = 50
    anchor = "center"
    qgr_stranded = qgr
    GenomicRanges::strand(qgr_stranded) = c(rep("+", 2), rep("-", 3))

    b_dt_center = viewGRangesWinSample_dt(bam_score, qgr_stranded, 50, anchor = "center")

    b_dt_center_uns = viewGRangesWinSample_dt(bam_score, qgr_stranded, 50, anchor = "center_unstranded")

    b_dt_left = viewGRangesWinSample_dt(bam_score, qgr_stranded, 50, anchor = "left")

    b_dt_left_uns = viewGRangesWinSample_dt(bam_score, qgr_stranded, 50, anchor = "left_unstranded")

    # b_dt = rbindlist(list(center = b_dt_center,
    #                       center_unstranded = b_dt_center_uns,
    #                       left = b_dt_left,
    #                       left_unstranded = b_dt_left_uns), use.names = T, idcol = "group")
    # b_dt[, facet_label := paste(id, strand)]
    # ggplot(b_dt[grepl("center", group)]) + geom_path(aes(x = x, y = y, color = group)) + facet_wrap("facet_label")
    # ggplot(b_dt[!grepl("center", group)]) + geom_path(aes(x = x, y = y, color = group)) + facet_wrap("facet_label")


    #verify stranded equal for + and not equal for -
    expect_true(all(b_dt_center[strand == "+"]$x == b_dt_center_uns[strand == "+"]$x))
    expect_true(all(!b_dt_center[strand == "-"]$x == b_dt_center_uns[strand == "-"]$x))
    expect_true(all(b_dt_left[strand == "+"]$x == b_dt_left_uns[strand == "+"]$x))
    expect_true(all(!b_dt_left[strand == "-"]$x == b_dt_left_uns[strand == "-"]$x))
    #verify center and left not equal
    expect_true(all(!b_dt_center$x == b_dt_left$x))
    expect_true(all(!b_dt_center_uns$x == b_dt_left_uns$x))
    #verify center and left ARE equal if minus min
    expect_true(all(b_dt_center$x - min(b_dt_center$x) == b_dt_left$x - min(b_dt_left$x)))
    expect_true(all(b_dt_center_uns$x - min(b_dt_center_uns$x) == b_dt_left_uns$x - min(b_dt_left_uns$x)))
})
