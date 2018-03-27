library(seqsetvis)
library(testthat)
library(GenomicRanges)
gr_a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1))
gr_b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1))
gr_c = GRanges("chr1", IRanges(8:10*10+5, 8:10*10+6))
setL = list("gr_a" = gr_a, "gr_b" = gr_b, "gr_c" = gr_c)

test_that("ssvFeatureVenn various set sizes, no outside counts. 1.", {
    p1 = ssvFeatureVenn(setL[1])
    p1
    expect_s3_class(p1, class = "ggplot")
})

test_that("ssvFeatureVenn various set sizes, no outside counts. 2.", {
    p2 = ssvFeatureVenn(setL[1:2])
    p2
    expect_s3_class(p2, class = "ggplot")
})

test_that("ssvFeatureVenn various set sizes, no outside counts. 3.", {
    p3 = ssvFeatureVenn(setL[1:3])
    p3
    expect_s3_class(p3, class = "ggplot")
})

test_that("ssvFeatureVenn various set sizes, with outside counts. 1.", {
    p1 = ssvFeatureVenn(setL[1], show_outside_count = T)
    p1
    expect_s3_class(p1, class = "ggplot")
})

test_that("ssvFeatureVenn various set sizes, with outside counts. 2.", {
    p2 = ssvFeatureVenn(setL[1:2], show_outside_count = T)
    p2
    expect_s3_class(p2, class = "ggplot")
})

test_that("ssvFeatureVenn various set sizes, with outside counts. 3.", {
    p3 = ssvFeatureVenn(setL[1:3], show_outside_count = T)
    p3
    expect_s3_class(p3, class = "ggplot")
})

test_that("ssvFeatureVenn can't plot more than 3 sets.", {
    setL = list("gr_a" = gr_a, "gr_b" = gr_b, "gr_c" = gr_c, "gr_d" = gr_c, "gr_e" = gr_a)
    expect_error(ssvFeatureVenn(setL[1:4]))
    expect_error(ssvFeatureVenn(setL[1:5]))
})

test_that("ssvFeatureVenn circle colors", {
    p1 = ssvFeatureVenn(setL, circle_colors = "blue")
    expect_s3_class(p1, class = "ggplot")
})

test_that("ssvFeatureEuler fill, no fill", {
    p1 = ssvFeatureEuler(setL)
    expect_s3_class(p1, class = "ggplot")

    p2 = ssvFeatureEuler(setL, fill_circles = F)
    expect_s3_class(p2, class = "ggplot")
})

test_that("ssvFeatureBinaryHeatmap no raster option", {
    p1 = ssvFeatureBinaryHeatmap(setL, raster_approximation = F)
    expect_s3_class(p1, class = "ggplot")
})

test_that("ssvFeatureBinaryHeatmap downsample", {
    p1 = ssvFeatureBinaryHeatmap(setL, raster_approximation = T, raster_width_min = 2)
    expect_s3_class(p1, class = "ggplot")

    p2 = ssvFeatureBinaryHeatmap(setL, raster_approximation = T, raster_height_min = 3)
    expect_s3_class(p2, class = "ggplot")

    p3 = ssvFeatureBinaryHeatmap(setL, raster_approximation = T, raster_height_min = 3, raster_width_min = 2)
    expect_s3_class(p3, class = "ggplot")

})

test_that("ssvFeatureVenn various other paramters don't throw error", {
    olap = ssvOverlapIntervalSets(list("gr_a" = gr_a, "gr_b" = gr_b, "gr_c" = gr_c))
    p = ssvFeatureVenn(olap, circle_colors = c("red", "blue", "green"), fill_alpha = .1,
                       counts_txt_size = 10, show_outside_count = T,
                       counts_as_labels = T)
    expect_s3_class(p, class = "ggplot")
    p = ssvFeatureVenn(olap, circle_colors = c("red", "blue", "green"), fill_circles = F,
                       counts_txt_size = 10, show_outside_count = T,
                       counts_as_labels = T)
    expect_s3_class(p, class = "ggplot")
})

test_that("col2hex", {
    string_colors = c("red", "blue", "green")
    rgb_colors = col2rgb(string_colors)
    hex_colors = c("#FF0000", "#0000FF", "#00FF00")
    expect_equal(col2hex(string_colors), hex_colors)
    expect_error(col2hex(rgb_colors))
    expect_error(col2hex("asdf"))
})

