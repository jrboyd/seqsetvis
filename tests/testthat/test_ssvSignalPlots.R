library(seqsetvis)
library(testthat)
library(GenomicRanges)

res = CTCF_in_10a_profiles_dt
res$sample = factor(res$sample)

test_that("ssvSignalBandedQuantiles works with proper inputs", {
    p_noBy = ssvSignalBandedQuantiles(res)
    p_wBy = ssvSignalBandedQuantiles(res, by_ = "sample")
    expect_s3_class(p_noBy, "ggplot")
    expect_s3_class(p_wBy, "ggplot")
})

test_that("ssvSignalBandedQuantiles different hsv setting", {
    p_rev = ssvSignalBandedQuantiles(res, by_ = "sample", hsv_reverse = T)
    expect_s3_class(p_rev, "ggplot")

    p_symm = ssvSignalBandedQuantiles(res, by_ = "sample", hsv_symmetric = T)
    expect_s3_class(p_symm, "ggplot")

    p_symmgray = ssvSignalBandedQuantiles(res, by_ = "sample", hsv_symmetric = T, hsv_grayscale = T)
    expect_s3_class(p_symmgray, "ggplot")

    p_unsymmgray = ssvSignalBandedQuantiles(res, by_ = "sample", hsv_symmetric = F, hsv_grayscale = T)
    expect_s3_class(p_unsymmgray, "ggplot")

    p_symmrev = ssvSignalBandedQuantiles(res, by_ = "sample", hsv_symmetric = T, hsv_reverse = T)
    expect_s3_class(p_symmrev, "ggplot")

    p_symmrev = ssvSignalBandedQuantiles(res, by_ = "sample", hsv_symmetric = T, hsv_reverse = T)
    expect_s3_class(p_symmrev, "ggplot")
})


test_that("ssvSignalScatterplot works with basic inputs", {
    p1 = ssvSignalScatterplot(res, x_name = "MCF10A", y_name = "MCF10AT1")
    expect_s3_class(p1, "ggplot")
    p2 = ssvSignalScatterplot(res, x_name = "MCF10A", y_name = "MCF10AT1", plot_type = "volcano")
    expect_s3_class(p2, "ggplot")
})

test_that("ssvSignalScatterplot throws correct error for bad x_name or y_name", {
expect_error(ssvSignalScatterplot(res, x_name = "badx", y_name = "MCF10AT1",
                                  value_variable = "x"), "badx")
expect_error(ssvSignalScatterplot(res, x_name = "MCF10A", y_name = "bady",
                                  value_variable = "x"))
})

test_that("ssvSignalScatterplot works with other inputs", {
    p1 = ssvSignalScatterplot(res, x_name = "MCF10A", y_name = "MCF10AT1",
                              value_variable = "x")
    expect_s3_class(p1, "ggplot")
    p2 = ssvSignalScatterplot(res, x_name = "MCF10A", y_name = "MCF10CA1",
                              value_function = median)
    expect_s3_class(p2, "ggplot")
    p3 = ssvSignalScatterplot(res, x_name = "1", y_name = "2",
                              value_function = median,
                              by_ = "sample", xy_variable = "id")
    memb = data.frame(id = unique(res$id), plotting_group = letters[1:5])
    p4 = ssvSignalScatterplot(bw_dt = res, x_name = "MCF10A", y_name = "MCF10A",
                              plotting_group = memb)
    expect_s3_class(p4, "ggplot")

    p5 = ssvSignalScatterplot(bw_dt = res, x_name = "MCF10A", y_name = "MCF10A", fixed_coords = F,
                              plotting_group = memb)
    expect_s3_class(p5, "ggplot")
})

test_that("ssvSignalScatterplot works with help enabled inputs", {
    p1 = ssvSignalScatterplot(res, x_name = "MCF10A", y_name = "MCF10AT1",
                              value_variable = "x", show_help = T)
    expect_s3_class(p1, "ggplot")

    p1v = ssvSignalScatterplot(res, x_name = "MCF10A", y_name = "MCF10AT1",
                               value_variable = "x", show_help = T, plot_type = "volcano")
    expect_s3_class(p1v, "ggplot")

    p2 = ssvSignalScatterplot(res, x_name = "MCF10A", y_name = "MCF10CA1",
                              value_function = median, show_help = T)
    expect_s3_class(p2, "ggplot")

    p3 = ssvSignalScatterplot(res, x_name = "1", y_name = "2",
                              value_function = median,
                              by_ = "sample", xy_variable = "id", show_help = T)

    memb = data.frame(id = unique(res$id), plotting_group = letters[1:5])
    p4 = ssvSignalScatterplot(bw_dt = res, x_name = "MCF10A", y_name = "MCF10A",
                              plotting_group = memb, show_help = T)
    expect_s3_class(p4, "ggplot")
})

test_that("ssvSignalHeatmap works with minimal inputs", {
    hidden = capture_output({
        p = ssvSignalHeatmap(res, nclust = 5)
    })
    expect_s3_class(p, "ggplot")
})

test_that("ssvSignalHeatmap works with maxCols and maxRows", {
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
        ssvSignalTrackplot(res[id %in% 1:3], facet_ = "sample"),
        ssvSignalTrackplot(res[id %in% 1:3],
                           facet_ = "sample~.",
                           facet_method = facet_grid),
        ssvSignalTrackplot(res[id %in% 1:3],
                           facet_ = paste("sample", "~", "id"), facet_method = facet_grid),
        ssvSignalTrackplot(res[id %in% 1:3]),
        ssvSignalTrackplot(res[id %in% 1:3], facet_ = "id"),
        ssvSignalTrackplot(res,
                           facet_ = "id", spline_n = 10)
    )
    lapply(test_plots, function(p1){
        expect_s3_class(p1, "ggplot")
    })
})

test_that("ssvSignalTrackplotAgg", {
    #from examples
    test_plots = list(
        ssvSignalTrackplotAgg(res) +
            labs(title = "agg regions by sample."),
        ssvSignalTrackplotAgg(res, spline_n = 10) +
            labs(title = "agg regions by sample, with spline smoothing."),
        ssvSignalTrackplotAgg(res[id %in% 1:10],
                              sample_ = "id", color_ = "id") +
            labs(title = "agg samples by region id (weird)"),
        ssvSignalTrackplotAgg(res[id %in% 1:10], sample_ = "id",
                              color_ = "id", spline_n = 10) +
            labs(title = "agg samples by region id (weird), with spline smoothing")
    )
    lapply(test_plots, function(p1){
        expect_s3_class(p1, "ggplot")
    })
})

