testthat::context("SignalPlots")
#all ssvSignalPlots functions should accept GRanges in addition to data.table
# as originally implemented
library(seqsetvis)
library(testthat)
library(GenomicRanges)

doTest_ssvSignalPlots = function(test_object){

    test_that("ssvSignalBandedQuantiles works with proper inputs", {
        p_noBy = ssvSignalBandedQuantiles(test_object)
        p_wBy = ssvSignalBandedQuantiles(test_object, by_ = "sample")
        expect_s3_class(p_noBy, "ggplot")
        expect_s3_class(p_wBy, "ggplot")
    })

    test_that("ssvSignalBandedQuantiles different hsv setting", {
        p_rev = ssvSignalBandedQuantiles(test_object,
                                         by_ = "sample", hsv_reverse = TRUE)
        expect_s3_class(p_rev, "ggplot")

        p_symm = ssvSignalBandedQuantiles(test_object,
                                          by_ = "sample", hsv_symmetric = TRUE)
        expect_s3_class(p_symm, "ggplot")

        p_symmgray = ssvSignalBandedQuantiles(
            test_object,
            by_ = "sample",
            hsv_symmetric = TRUE,
            hsv_grayscale = TRUE
        )
        expect_s3_class(p_symmgray, "ggplot")

        p_unsymmgray = ssvSignalBandedQuantiles(
            test_object,
            by_ = "sample",
            hsv_symmetric = FALSE,
            hsv_grayscale = TRUE
        )
        expect_s3_class(p_unsymmgray, "ggplot")

        p_symmrev = ssvSignalBandedQuantiles(
            test_object,
            by_ = "sample",
            hsv_symmetric = TRUE,
            hsv_reverse = TRUE
        )
        expect_s3_class(p_symmrev, "ggplot")

        p_symmrev = ssvSignalBandedQuantiles(
            test_object,
            by_ = "sample",
            hsv_symmetric = TRUE,
            hsv_reverse = TRUE
        )
        expect_s3_class(p_symmrev, "ggplot")
    })


    test_that("ssvSignalScatterplot works with basic inputs", {
        p1 = ssvSignalScatterplot(test_object,
                                  x_name = "MCF10A_CTCF",
                                  y_name = "MCF10AT1_CTCF")
        expect_s3_class(p1, "ggplot")
        p2 = ssvSignalScatterplot(
            test_object,
            x_name = "MCF10A_CTCF",
            y_name = "MCF10AT1_CTCF",
            plot_type = "volcano"
        )
        expect_s3_class(p2, "ggplot")
    })

    test_that("ssvSignalScatterplot throws correct error for bad x_name or y_name",
              {
                  expect_error(
                      ssvSignalScatterplot(
                          test_object,
                          x_name = "badx",
                          y_name = "MCF10AT1_CTCF",
                          value_variable = "x"
                      ),
                      "badx"
                  )
                  expect_error(
                      ssvSignalScatterplot(
                          test_object,
                          x_name = "MCF10A_CTCF",
                          y_name = "bady",
                          value_variable = "x"
                      )
                  )
              })

    test_that("ssvSignalScatterplot works with other inputs", {
        p1 = ssvSignalScatterplot(
            test_object,
            x_name = "MCF10A_CTCF",
            y_name = "MCF10AT1_CTCF",
            value_variable = "x"
        )
        expect_s3_class(p1, "ggplot")
        p2 = ssvSignalScatterplot(
            test_object,
            x_name = "MCF10A_CTCF",
            y_name = "MCF10CA1_CTCF",
            value_function = median
        )
        expect_s3_class(p2, "ggplot")
        p3 = ssvSignalScatterplot(
            test_object,
            x_name = "1",
            y_name = "2",
            value_function = median,
            by_ = "sample",
            xy_variable = "id"
        )
        memb = data.frame(id = unique(test_object$id),  group = letters[1:5])
        p4 = ssvSignalScatterplot(
            bw_data = test_object,
            x_name = "MCF10A_CTCF",
            y_name = "MCF10A_CTCF",
            color_table = memb
        )
        expect_s3_class(p4, "ggplot")

        p5 = ssvSignalScatterplot(
            bw_data = test_object,
            x_name = "MCF10A_CTCF",
            y_name = "MCF10A_CTCF",
            fixed_coords = FALSE,
            color_table = memb
        )
        expect_s3_class(p5, "ggplot")
    })

    test_that("ssvSignalScatterplot works with help enabled inputs", {
        p1 = ssvSignalScatterplot(
            test_object,
            x_name = "MCF10A_CTCF",
            y_name = "MCF10AT1_CTCF",
            value_variable = "x",
            show_help = TRUE
        )
        expect_s3_class(p1, "ggplot")

        p1v = ssvSignalScatterplot(
            test_object,
            x_name = "MCF10A_CTCF",
            y_name = "MCF10AT1_CTCF",
            value_variable = "x",
            show_help = TRUE,
            plot_type = "volcano"
        )
        expect_s3_class(p1v, "ggplot")

        p2 = ssvSignalScatterplot(
            test_object,
            x_name = "MCF10A_CTCF",
            y_name = "MCF10CA1_CTCF",
            value_function = median,
            show_help = TRUE
        )
        expect_s3_class(p2, "ggplot")

        p3 = ssvSignalScatterplot(
            test_object,
            x_name = "1",
            y_name = "2",
            value_function = median,
            by_ = "sample",
            xy_variable = "id",
            show_help = TRUE
        )

        memb = data.frame(id = unique(test_object$id),  group = letters[1:5])
        p4 = ssvSignalScatterplot(
            bw_data = test_object,
            x_name = "MCF10A_CTCF",
            y_name = "MCF10A_CTCF",
            color_table = memb,
            show_help = TRUE
        )
        expect_s3_class(p4, "ggplot")
    })

    test_that("ssvSignalHeatmap works with minimal inputs", {
        hidden = capture_output({
            p = ssvSignalHeatmap(test_object, nclust = 5)
        })
        expect_s3_class(p, "ggplot")
    })

    test_that("ssvSignalHeatmap works with maxCols and maxRows", {
        expect_message({
            p_max_rows = ssvSignalHeatmap(test_object, nclust = 2, max_rows = 8)
        }, regexp = "92 rows were discarded")
        expect_s3_class(p_max_rows, "ggplot")

        expect_message({
            p_max_cols = ssvSignalHeatmap(test_object, nclust = 2, max_cols = 2)
        }, regexp = "12 columns were discarded")
        expect_s3_class(p_max_cols, "ggplot")
    })

    test_that("ssvSignalHeatmap works with manual clustering", {
        hidden = capture_output({
            clust_dt = ssvSignalClustering(test_object, nclust = 2)
            p1 = ssvSignalHeatmap(clust_dt)
            p2 = ssvSignalHeatmap(clust_dt, perform_clustering = "no")
        })
        expect_s3_class(p1, "ggplot")
        expect_s3_class(p2, "ggplot")
    })

    test_that("ssvSignalLineplot", {
        #from examples
        sub_to = test_object[test_object$id %in% 1:3]
        test_plots = list(
            ssvSignalLineplot(sub_to, facet_ = "sample"),
            ssvSignalLineplot(sub_to,
                              facet_ = "sample~.",
                              facet_method = facet_grid),
            ssvSignalLineplot(
                sub_to,
                facet_ = paste("sample", "~", "id"),
                facet_method = facet_grid
            ),
            ssvSignalLineplot(sub_to),
            ssvSignalLineplot(sub_to, facet_ = "id"),
            ssvSignalLineplot(sub_to,
                              facet_ = "id", spline_n = 10)
        )
        lapply(test_plots, function(p1) {
            expect_s3_class(p1, "ggplot")
        })
    })

    test_that("ssvSignalLineplotAgg", {
        #from examples
        test_plots = list(
            ssvSignalLineplotAgg(test_object) +
                labs(title = "agg regions by sample."),
            ssvSignalLineplotAgg(test_object, spline_n = 10) +
                labs(title = "agg regions by sample, with spline smoothing."),
            ssvSignalLineplotAgg(test_object[test_object$id %in% 1:10],
                                 sample_ = "id", color_ = "id") +
                labs(title = "agg samples by region id (weird)"),
            ssvSignalLineplotAgg(
                test_object[test_object$id %in% 1:10],
                sample_ = "id",
                color_ = "id",
                spline_n = 10
            ) +
                labs(title = paste("agg samples by region id",
                                   "(weird), with spline smoothing"))
        )
        lapply(test_plots, function(p1) {
            expect_s3_class(p1, "ggplot")
        })
    })
}

test_object = CTCF_in_10a_profiles_dt
test_object$sample = factor(test_object$sample)
res_dt = test_object
res_gr = GRanges(res_dt)

doTest_ssvSignalPlots(test_object = res_dt)
doTest_ssvSignalPlots(test_object = res_gr)

test_that("ssvSignalHeatmap 1 clust", {
    expect_failure(expect_error(ssvSignalHeatmap(res_dt, nclust = 1)))
})

test_that("ssvSignalHeatmap too many clust", {
    t_dt = res_dt[id %in% unique(id)[1:5]]
    expect_failure(expect_error(ssvSignalHeatmap(t_dt, nclust = 5)))
    expect_failure(expect_error(ssvSignalHeatmap(t_dt, nclust = 6)))
})
