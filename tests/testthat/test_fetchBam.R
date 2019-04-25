testthat::context("FetchBam")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
library(data.table)

qgr = CTCF_in_10a_overlaps_gr[1:5]
qgr = centerFixedSizeGRanges(qgr, 500)
#bed used to intersect bam
# rtracklayer::export.bed(qgr, con = "ctcf_5.bed")
bam_file = system.file("extdata/test.bam", package = "seqsetvis", mustWork = TRUE)
xls_file = system.file("extdata/test_peaks.xls", package = "seqsetvis", mustWork = TRUE)

# ssvFetchBam(bam_file, qgr)

test_that("fragLen_fromMacs2Xls", {
    fl = fragLen_fromMacs2Xls(xls_file)
    expect_equal(fl, 187)
})

test_that("fetchBam auto fragLen", {
    expect_message(seqsetvis:::fetchBam(bam_file, qgr = qgr), "fragLen was calculated as:")
})

test_that("fragLen_calcStranded parameters matter", {
    # res2 = fragLen_calcStranded(bam_file, qgr = qgr, n_regions = 5, include_plot_in_output = TRUE)
    # res1 = fragLen_calcStranded(bam_file, qgr = qgr, n_regions = 1, include_plot_in_output = TRUE)
    #
    # f1 = res2[[1]]
    # f2 = res1[[1]]
    # expect_failure(expect_equal(f1, f2))

    res5 = fragLen_calcStranded(bam_file, qgr = qgr, test_fragLen = 180, include_plot_in_output = TRUE)
    res6 = fragLen_calcStranded(bam_file, qgr = qgr, test_fragLen = 200, include_plot_in_output = TRUE)

    f5 = res5[[1]]
    f6 = res6[[1]]
    expect_failure(expect_equal(f5, f6))
})

test_that("viewGRangesWinSample_dt strand and position functions", {
    bam_score = seqsetvis:::fetchBam(bam_file, qgr = qgr)
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

# #strand of qgr should be ignored for pileup in favor of target strand
# #strand of qgr should be used for feature orientation
# test_that("viewGRangesWinSample_dt strand and position functions", {
#     bam_score_unstranded = fetchBam(bam_file, qgr = qgr)
#     score_gr = bam_score
#     window_size = 50
#
#     qgr_stranded = qgr
#     GenomicRanges::strand(qgr_stranded) = c(rep("+", 2), rep("-", 3))
#
#     bam_score_stranded = fetchBam(bam_file, qgr = qgr_stranded)
#     plot(bam_score_unstranded$score, bam_score_stranded$score)
# })


test_that("ssvFetchBam query GRanges output id set", {
    skip_on_os("windows")
    test_gr = qgr
    gr_sample = ssvFetchBam(bam_file, win_size = 5, qgr = test_gr, return_data.table = TRUE)
    expect_true("id" %in% colnames(gr_sample))

    names(test_gr) = NULL
    gr_sample = ssvFetchBam(bam_file, win_size = 5, qgr = test_gr, return_data.table = TRUE)
    expect_true("id" %in% colnames(gr_sample))

    test_gr$id = seq_along(test_gr) #when id set, id in output is missing
    gr_sample = ssvFetchBam(bam_file, win_size = 5, qgr = test_gr, return_data.table = TRUE)
    expect_true("id" %in% colnames(gr_sample))

    test_gr = qgr
    test_gr$id = seq_along(test_gr) #when id set, id in output is missing
    gr_sample = ssvFetchBam(bam_file, win_size = 5, qgr = test_gr, return_data.table = TRUE)
    expect_true("id" %in% colnames(gr_sample))

})

test_that("ssvFetchBam query GRanges name gets used", {
    test_gr = qgr
    names(test_gr) = paste("test_region", seq_along(test_gr))
    gr_sample = ssvFetchBam(bam_file, win_size = 5, qgr = test_gr, return_data.table = TRUE)
    expect_true(all(unique(gr_sample$id) == unique(names(test_gr))))
})

test_that("ssvFetchBam bam_file as data.frame/table", {
    test_gr = CTCF_in_10a_narrowPeak_grs$MCF10A_CTCF
    gr_sample = ssvFetchBam(
        data.frame(rep(bam_file,3),
                   colors = c("red", "green", "blue"),
                   sample = c("10a", "10b", "10c")),
        win_size = 5, qgr = test_gr,
        return_data.table = TRUE)
    expect_true(all(levels(gr_sample$colors) %in%
                        c("red", "green", "blue")))
    expect_true(all(levels(gr_sample$sample) %in%
                        c("10a", "10b", "10c")))
})

test_that("ssvFetchBam sample method correct bins", {
    skip_on_os("windows")
    test_qgr2 = qgr
    gr_sample = ssvFetchBam(bam_file, win_size = 5, qgr = test_qgr2)
    #for non-overlapping, all reduce granges width should be equal to input qgr
    expect_true(all(width(reduce(gr_sample)) == width(test_qgr2)))
    #for non-overlapping, all granges width should be equal to win_size
    expect_true(all(width(gr_sample) == 5))

    #same as above but width are auto adjusted to be divisible by win_size
    gr_sample = ssvFetchBam(bam_file, win_size = 3, qgr = test_qgr2)
    expect_true(all(width(reduce(gr_sample)) == 501))
    expect_true(all(width(gr_sample) == 3))

    ###varying widths
    end(test_qgr2) = end(test_qgr2) + seq_along(test_qgr2)*2
    gr_sample = ssvFetchBam(bam_file, win_size = 5, qgr = test_qgr2, win_method = "sample")
    #for non-overlapping, all reduce granges width should be equal to input qgr
    expect_true(all(width(reduce(gr_sample)) == 510))
    #for non-overlapping, all granges width should be equal to win_size
    expect_true(all(width(gr_sample) == 5))

    #same as above but individual widths now vary while total is unchanged
    gr_sample = ssvFetchBam(bam_file, win_size = 4, qgr = test_qgr2, win_method = "sample")
    expect_true(all(width(reduce(gr_sample)) == 508))
    expect_true(all(width(gr_sample) == 4))
})

test_that("ssvFetchBam summary method correct bins", {
    skip_on_os("windows")
    test_qgr2 = qgr

    ###invariant widths
    gr_sample = ssvFetchBam(bam_file, win_size = 5, qgr = test_qgr2, win_method = "summary")
    #for non-overlapping, all reduce granges width should be equal to input qgr
    expect_true(all(width(reduce(gr_sample)) == width(test_qgr2)))
    #for non-overlapping, all granges width should be equal to win_size
    expect_true(all(width(gr_sample) == 500 / 5))

    #same as above but individual widths now vary while total is unchanged
    gr_sample = ssvFetchBam(bam_file, win_size = 3, qgr = test_qgr2, win_method = "summary")
    expect_true(all(width(reduce(gr_sample)) == width(test_qgr2)))
    expect_true(all(width(gr_sample) %in% 166:167))

    ###varying widths
    end(test_qgr2) = end(test_qgr2) + seq_along(test_qgr2)*2
    gr_sample = ssvFetchBam(bam_file, win_size = 5, qgr = test_qgr2, win_method = "summary")
    #for non-overlapping, all reduce granges width should be equal to input qgr
    expect_true(all(width(reduce(gr_sample)) == width(test_qgr2)))
    #for non-overlapping, all granges width should be equal to win_size
    expect_true(all(width(gr_sample) %in% 100:102))

    #same as above but individual widths now vary while total is unchanged
    gr_sample = ssvFetchBam(bam_file, win_size = 3, qgr = test_qgr2, win_method = "summary")
    expect_true(all(width(reduce(gr_sample)) == width(test_qgr2)))
    expect_true(all(width(gr_sample) %in% 167:170))
})

test_that(".rm_dupes removes duplicates", {
    dt = data.table(which_label = 1:10, seqnames = "chr1", strand = c("+", "-"), start = 1:10, end = 1:10+10)
    dt[, width := end - start + 1]
    make_dupes = seq_len(nrow(dt))
    make_dupes = rep(make_dupes, make_dupes)
    dt = dt[make_dupes]
    dtl = lapply(1:10, function(x)seqsetvis:::.rm_dupes(dt, max_dupes = x))
    # max_dupes 1 should yield 10 unique entries
    expect_true(all(dtl[[1]]$which_label == 1:10))
    # max_dupes 10 should perform no dupe removal in this case
    expect_true(all(dtl[[10]]$which_label == make_dupes))
    lens = sapply(dtl, nrow)
    # as max_dupes increases, returned entries increases
    expect_true(all(lens[-length(lens)] - lens[-1] < 0))
})

test_that("ssvFetchBam removes duplicates if max_dupes set", {
    rFull = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, win_method = "summary")$y
    r1 = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, win_method = "summary", max_dupes = 1)$y
    r2 = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, win_method = "summary", max_dupes = 2)$y
    expect_true(all(r2 >= r1))
    expect_true(all(!r2 < r1))
    expect_true(all(rFull >= r1))
    expect_true(all(!rFull < r1))
})

test_that("ssvFetchBam strand of qgr doesn't matter", {
    qres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr[1], fragLens = 200)
    pgr = qgr
    strand(pgr) = "+"
    pres = ssvFetchBam(bam_file, win_size = 5, qgr = pgr[1], fragLens = 200)
    ngr = qgr
    strand(ngr) = "-"
    nres = ssvFetchBam(bam_file, win_size = 5, qgr = ngr[1], fragLens = 200)
    expect_equal(pres$y,  qres$y)
    expect_equal(nres$y,  rev(qres$y))
})

test_that("ssvFetchBam strand of targt_strand does matter", {
    strand(qgr) = "*"
    qres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, fragLens = 200, target_strand = "*")
    pres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, fragLens = 200, target_strand = "+")
    nres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, fragLens = 200, target_strand = "-")

    expect_true(!all(pres$y == qres$y))
    expect_true(!all(nres$y == qres$y))
    expect_true(all(nres$y + pres$y == qres$y))

    strand(qgr) = "+"
    qres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, fragLens = 200, target_strand = "*")
    pres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, fragLens = 200, target_strand = "+")
    nres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, fragLens = 200, target_strand = "-")

    expect_true(!all(pres$y == qres$y))
    expect_true(!all(nres$y == qres$y))
    expect_true(all(nres$y + pres$y == qres$y))

    strand(qgr) = "-"
    qres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, fragLens = 200, target_strand = "*")
    pres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, fragLens = 200, target_strand = "+")
    nres = ssvFetchBam(bam_file, win_size = 5, qgr = qgr, fragLens = 200, target_strand = "-")

    expect_true(!all(pres$y == qres$y))
    expect_true(!all(nres$y == qres$y))
    expect_true(all(nres$y + pres$y == qres$y))
})

test_that("ssvFetchBam target_strand of both - sample", {
    strand(qgr) = "*"
    res_both = ssvFetchBam(bam_file, win_size = 5, qgr = qgr[1], fragLens = 200,
                           target_strand = "both", return_data.table = TRUE)
    strand(qgr) = "+"
    res_pos = ssvFetchBam(bam_file, win_size = 5, qgr = qgr[1], fragLens = 200,
                          target_strand = "both", return_data.table = TRUE)
    strand(qgr) = "-"
    res_neg = ssvFetchBam(bam_file, win_size = 5, qgr = qgr[1], fragLens = 200,
                          target_strand = "both", return_data.table = TRUE)
    res_neg[, x := -x]
    res_neg = res_neg[order(x)][order(id)][order(strand)]

    expect_true(nrow(res_both[strand == "+"]) == nrow(res_both[strand == "-"]))
    expect_true(!all(res_both[strand == "+"]$y == res_both[strand == "-"]))
    expect_true(all(res_both[strand == "+"] == res_pos[strand == "+"]))
    expect_true(all(res_both[strand == "-"] == res_pos[strand == "-"]))
    expect_true(all(res_both[strand == "+"]$y == res_neg[strand == "-"]$y))
    expect_true(all(res_both[strand == "-"]$y == res_neg[strand == "+"]$y))

})

test_that("ssvFetchBam target_strand of both - summary", {
    strand(qgr) = "*"
    res_both = ssvFetchBam(bam_file, win_size = 5, qgr = qgr[1], fragLens = 200,
                           target_strand = "both", return_data.table = TRUE,
                           win_method = "summary")
    strand(qgr) = "+"
    res_pos = ssvFetchBam(bam_file, win_size = 5, qgr = qgr[1], fragLens = 200,
                          target_strand = "both", return_data.table = TRUE,
                          win_method = "summary")
    strand(qgr) = "-"
    res_neg = ssvFetchBam(bam_file, win_size = 5, qgr = qgr[1], fragLens = 200,
                          target_strand = "both", return_data.table = TRUE,
                          win_method = "summary")
    res_neg[, x := -x]
    res_neg = res_neg[order(x)][order(id)][order(strand)]

    expect_true(nrow(res_both[strand == "+"]) == nrow(res_both[strand == "-"]))
    expect_true(!all(res_both[strand == "+"]$y == res_both[strand == "-"]))
    expect_true(all(res_both[strand == "+"] == res_pos[strand == "+"]))
    expect_true(all(res_both[strand == "-"] == res_pos[strand == "-"]))
    expect_true(all(res_both[strand == "+"]$y == res_neg[strand == "-"]$y))
    expect_true(all(res_both[strand == "-"]$y == res_neg[strand == "+"]$y))

})

test_that("ssvFetchBam single column data.table", {
    qdt = data.table(file = bam_file)
    res = ssvFetchBam(qdt,
                      qgr = qgr[1],
                      win_size = 5, fragLens = 200,
                           target_strand = "both", return_data.table = TRUE,
                           win_method = "summary")
    expect_equal(res$sample[1], bam_file)
})
