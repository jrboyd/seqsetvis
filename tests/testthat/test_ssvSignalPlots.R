library(seqsetvis)
library(testthat)
library(GenomicRanges)
# library(rtracklayer)
test_bw = system.file("extdata/test_bigwigs/test_loading.bw", package = "seqsetvis", mustWork = T)
pos = c(20, 180, 210, 440, 520, 521)
region_size = 30
test_qgr = GRanges("chrTest", IRanges(pos+1, pos + region_size))
exp_colnames = c("seqnames", "start", "end", "width", "strand", "id", "y", "x")

test_that("fetchWindowedBigwig return expected for valid even win_size", {
    #these should all work cleanly
    for(win in c(2, 6, 10, 30)){
        bw_dt = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = test_qgr)
        expect_is(bw_dt, "data.table")
        expect_equal(colnames(bw_dt), exp_colnames)
        for(tid in unique(bw_dt$id)){
            test_dt = bw_dt[id == tid]
            expect_equal(nrow(test_dt), region_size / win) #expected number of regions
            bw_gr = GRanges(test_dt)
            expect_equal(length(findOverlaps(bw_gr, bw_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("fetchWindowedBigwig return expected for valid odd win_size", {
    for(win in c(1, 3, 5, 15)){
        bw_dt = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = test_qgr)
        expect_is(bw_dt, "data.table")
        expect_equal(colnames(bw_dt), exp_colnames)
        for(tid in unique(bw_dt$id)){
            test_dt = bw_dt[id == tid]
            expect_equal(nrow(test_dt), region_size / win) #expected number of regions
            bw_gr = GRanges(test_dt)
            expect_equal(length(findOverlaps(bw_gr, bw_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("fetchWindowedBigwig use GRanges names as id", {
    for(win in c(1, 3, 5, 15)){
        qgr = test_qgr
        names(qgr) = paste0("myNames_", seq_along(qgr))
        bw_dt = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = qgr)
        expect_is(bw_dt, "data.table")
        expect_equal(colnames(bw_dt), exp_colnames)
        expect_true(all(grepl("myNames", bw_dt$id)))
        for(tid in unique(bw_dt$id)){
            test_dt = bw_dt[id == tid]
            expect_equal(nrow(test_dt), region_size / win) #expected number of regions
            bw_gr = GRanges(test_dt)
            expect_equal(length(findOverlaps(bw_gr, bw_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("fetchWindowedBigwig patches missing values", {
    for(win in c(1, 5, 20)){
        qgr = GRanges("chrTest", IRanges(1, 2000))
        bw_dt = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = qgr)
        expect_is(bw_dt, "data.table")
        expect_equal(colnames(bw_dt), exp_colnames)
    }
})

test_that("fetchWindowedBigwig throws error if widths aren't divisble by win_size", {
    mix_width_gr = test_qgr
    end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)
    for(win in c(7, 11, 4, 8, 20)){
        expect_error({
            fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = mix_width_gr)
        })
    }
})



test_that("fetchWindowedBigwig throws warning if widths vary", {
    mix_width_gr = test_qgr
    end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)*3
    for(win in c(1, 3)){
        expect_warning({
            fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = mix_width_gr)
        })
    }
})

test_that("fetchWindowedBigwigList works with proper inputs", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr)})
    expect_s3_class(res, "data.table")
    expect_equal(colnames(res), c(exp_colnames, "sample"))
})

test_that("fetchWindowedBigwigList can set variable name", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr,
                                                           bw_variable_name = "group")})
    expect_s3_class(res, "data.table")
    expect_equal(colnames(res), c(exp_colnames, "group"))
})

test_that("fetchWindowedBigwigList duplicate names throws error", {
    bw_files = rep(test_bw, 3)
    expect_error(
        fetchWindowedBigwigList(bw_files = bw_files,
                                win_size = 3,
                                qgr = test_qgr)
    )
})

test_that("ssvSignalBandedQuantiles works with proper inputs", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    for(win in c(1, 3)){
        hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                               win_size = win,
                                                               qgr = test_qgr,
                                                               bw_variable_name = "group")})
        p_noBy = ssvSignalBandedQuantiles(res)
        p_wBy = ssvSignalBandedQuantiles(res, by_ = "group")
        expect_s3_class(p_noBy, "ggplot")
        expect_s3_class(p_wBy, "ggplot")
    }
})

test_that("ssvSignalBandedQuantiles different hsv setting", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    for(win in c(1, 3)){
        hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                               win_size = win,
                                                               qgr = test_qgr,
                                                               bw_variable_name = "group")})

        res$group = factor(res$group)

        p_rev = ssvSignalBandedQuantiles(res, by_ = "group", hsv_reverse = T)
        expect_s3_class(p_rev, "ggplot")

        p_symm = ssvSignalBandedQuantiles(res, by_ = "group", hsv_symmetric = T)
        expect_s3_class(p_symm, "ggplot")

        p_symmgray = ssvSignalBandedQuantiles(res, by_ = "group", hsv_symmetric = T, hsv_grayscale = T)
        expect_s3_class(p_symmgray, "ggplot")

        p_unsymmgray = ssvSignalBandedQuantiles(res, by_ = "group", hsv_symmetric = F, hsv_grayscale = T)
        expect_s3_class(p_unsymmgray, "ggplot")

        p_symmrev = ssvSignalBandedQuantiles(res, by_ = "group", hsv_symmetric = T, hsv_reverse = T)
        expect_s3_class(p_symmrev, "ggplot")

        p_symmrev = ssvSignalBandedQuantiles(res, by_ = "group", hsv_symmetric = T, hsv_reverse = T)
        expect_s3_class(p_symmrev, "ggplot")
    }
})


test_that("ssvSignalScatterplot works with basic inputs", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr)})
    p1 = ssvSignalScatterplot(res, x_name = "bw_1", y_name = "bw_2")
    expect_s3_class(p1, "ggplot")
    p2 = ssvSignalScatterplot(res, x_name = "bw_1", y_name = "bw_2", plot_type = "volcano")
    expect_s3_class(p2, "ggplot")
})

test_that("ssvSignalScatterplot works with other inputs", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr)})
    p1 = ssvSignalScatterplot(res, x_name = "bw_1", y_name = "bw_2",
                              value_variable = "x")
    expect_s3_class(p1, "ggplot")
    p2 = ssvSignalScatterplot(res, x_name = "bw_1", y_name = "bw_3",
                              value_function = median)
    expect_s3_class(p2, "ggplot")
    p3 = ssvSignalScatterplot(res, x_name = "region_1", y_name = "region_2",
                              value_function = median,
                              by_ = "sample", xy_variable = "id")
    memb = data.frame(id = unique(res$id), plotting_group = letters[1:6])
    p4 = ssvSignalScatterplot(bw_dt = res, x_name = "bw_1", y_name = "bw_1",
                              plotting_group = memb)
    expect_s3_class(p4, "ggplot")

    p5 = ssvSignalScatterplot(bw_dt = res, x_name = "bw_1", y_name = "bw_1", fixed_coords = F,
                              plotting_group = memb)
    expect_s3_class(p5, "ggplot")
})

test_that("ssvSignalScatterplot works with help enabled inputs", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr)})
    hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr)})
    p1 = ssvSignalScatterplot(res, x_name = "bw_1", y_name = "bw_2",
                              value_variable = "x", show_help = T)
    expect_s3_class(p1, "ggplot")

    p1v = ssvSignalScatterplot(res, x_name = "bw_1", y_name = "bw_2",
                               value_variable = "x", show_help = T, plot_type = "volcano")
    expect_s3_class(p1v, "ggplot")

    p2 = ssvSignalScatterplot(res, x_name = "bw_1", y_name = "bw_3",
                              value_function = median, show_help = T)
    expect_s3_class(p2, "ggplot")

    p3 = ssvSignalScatterplot(res, x_name = "region_1", y_name = "region_2",
                              value_function = median,
                              by_ = "sample", xy_variable = "id", show_help = T)

    memb = data.frame(id = unique(res$id), plotting_group = letters[1:6])
    p4 = ssvSignalScatterplot(bw_dt = res, x_name = "bw_1", y_name = "bw_1",
                              plotting_group = memb, show_help = T)
    expect_s3_class(p4, "ggplot")
})

a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1:7 + 10))
b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1:6 + 8))

test_that("centerFixedSizeGRanges final width equals fixed_size", {
    expect_equal(width(centerFixedSizeGRanges(a, 2)), rep(2, length(a)))
    expect_equal(width(centerFixedSizeGRanges(a, 5)), rep(5, length(a)))
    expect_equal(width(centerFixedSizeGRanges(a, 80)), rep(80, length(a)))
    expect_equal(width(centerFixedSizeGRanges(b, 1)), rep(1, length(b)))
    expect_equal(width(centerFixedSizeGRanges(b, 7)), rep(7, length(b)))
    expect_equal(width(centerFixedSizeGRanges(b, 50)), rep(50, length(b)))
})
test_that("centerFixedSizeGRanges size shifts are centered", {
    a_larger = centerFixedSizeGRanges(grs = a, fixed_size = max(width(a)) + 2)
    expect_true(all(start(a_larger) < start(a)))
    expect_true(all(end(a_larger) > end(a)))
    a_smaller = centerFixedSizeGRanges(grs = a_larger, fixed_size = min(width(a)) - 2)
    expect_true(all(start(a_larger) < start(a_smaller)))
    expect_true(all(end(a_larger) > end(a_smaller)))
})
test_that("centerFixedSizeGRanges size sifts are reversible", {
    #primarily concerened about impacts of rounding
    a4 = centerFixedSizeGRanges(grs = a, fixed_size = 4)
    a5 = centerFixedSizeGRanges(grs = a, fixed_size = 5)
    a6 = centerFixedSizeGRanges(grs = a, fixed_size = 6)
    a7 = centerFixedSizeGRanges(grs = a, fixed_size = 7)
    a9 = centerFixedSizeGRanges(grs = a, fixed_size = 9)
    #derivations
    a7_from5 = centerFixedSizeGRanges(grs = a5, fixed_size = 7)
    a4_from5 = centerFixedSizeGRanges(grs = a5, fixed_size = 4)
    a6_from4 = centerFixedSizeGRanges(grs = a4, fixed_size = 6)
    a9_from4 = centerFixedSizeGRanges(grs = a4, fixed_size = 9)
    a7_from4 = centerFixedSizeGRanges(grs = a4, fixed_size = 7)
    a5_from4 = centerFixedSizeGRanges(grs = a4, fixed_size = 5)
    #reversals
    a5_from7rev = centerFixedSizeGRanges(grs = a7_from5, fixed_size = 5)
    a5_from4rev = centerFixedSizeGRanges(grs = a4_from5, fixed_size = 5)
    a4_from6rev = centerFixedSizeGRanges(grs = a6_from4, fixed_size = 4)
    a4_from9rev = centerFixedSizeGRanges(grs = a9_from4, fixed_size = 4)
    #recover even from even
    expect_true(all(a4 == a4_from6rev))
    expect_true(all(a6 == a6_from4))

    #recover odd from even
    expect_true(all(a5 == a5_from4rev))
    expect_true(all(a9 == a9_from4))#fail
    expect_true(all(a7 == a7_from4))
    expect_true(all(a5 == a5_from4))

    #recover odd from odd
    expect_true(all(a5 == a5_from7rev))
    expect_true(all(a7 == a7_from5))

    #recover even from odd
    expect_true(all(a4 == a4_from5))
    expect_true(all(a4 == a4_from9rev))
})

test_that("ssvSignalHeatmap works with minimal inputs", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({
        res = fetchWindowedBigwigList(bw_files = bw_files,
                                      win_size = 3,
                                      qgr = test_qgr)
    })
    hidden = capture_output({
        p = ssvSignalHeatmap(res, nclust = 5)
    })
    expect_s3_class(p, "ggplot")
})

test_that("ssvSignalHeatmap works with maxCols and maxRows", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({
        res = fetchWindowedBigwigList(bw_files = bw_files,
                                      win_size = 3,
                                      qgr = test_qgr)
    })
    hidden = capture_output({
        p_max_rows = ssvSignalHeatmap(res, nclust = 2, max_rows = 4)
    })
    expect_s3_class(p_max_rows, "ggplot")

    hidden = capture_output({
        p_max_cols = ssvSignalHeatmap(res, nclust = 2, max_cols = 2)
    })
    expect_s3_class(p_max_cols, "ggplot")
})

test_that("ssvSignalHeatmap works with manual clustering", {
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({
        res = fetchWindowedBigwigList(bw_files = bw_files,
                                      win_size = 3,
                                      qgr = test_qgr)
    })
    hidden = capture_output({
        clust_dt = ssvSignalClustering(res, nclust = 2)
        p1 = ssvSignalHeatmap(clust_dt)
        p2 = ssvSignalHeatmap(clust_dt, perform_clustering = "no")
    })
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})

test_that("ssvSignalTrackplot", {
    #from examples
    test_plots = list(
        ssvSignalTrackplot(CTCF_in_10a_profiles_dt[id %in% 1:3], facet_ = "sample"),
        ssvSignalTrackplot(CTCF_in_10a_profiles_dt[id %in% 1:3],
                           facet_ = "sample~.",
                           facet_method = facet_grid),
        ssvSignalTrackplot(CTCF_in_10a_profiles_dt[id %in% 1:3],
                           facet_ = paste("sample", "~", "id"), facet_method = facet_grid),
        ssvSignalTrackplot(CTCF_in_10a_profiles_dt[id %in% 1:3]),
        ssvSignalTrackplot(CTCF_in_10a_profiles_dt[id %in% 1:3], facet_ = "id"),
        ssvSignalTrackplot(CTCF_in_10a_profiles_dt[id %in% 1:3],
                           facet_ = "id", spline_n = 10)
    )
    lapply(test_plots, function(p1){
        expect_s3_class(p1, "ggplot")
    })
})

test_that("ssvSignalTrackplotAgg", {
    #from examples
    test_plots = list(
        ssvSignalTrackplotAgg(CTCF_in_10a_profiles_dt) +
            labs(title = "agg regions by sample."),
        ssvSignalTrackplotAgg(CTCF_in_10a_profiles_dt, spline_n = 10) +
            labs(title = "agg regions by sample, with spline smoothing."),
        ssvSignalTrackplotAgg(CTCF_in_10a_profiles_dt[id %in% 1:10],
                              sample_ = "id", color_ = "id") +
            labs(title = "agg samples by region id (weird)"),
        ssvSignalTrackplotAgg(CTCF_in_10a_profiles_dt[id %in% 1:10], sample_ = "id",
                              color_ = "id", spline_n = 10) +
            labs(title = "agg samples by region id (weird), with spline smoothing")
    )
    lapply(test_plots, function(p1){
        expect_s3_class(p1, "ggplot")
    })
})