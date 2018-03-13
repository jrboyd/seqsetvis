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

test_that("ssvFeatureBinaryHeatmap raster options", {
    p1 = ssvFeatureBinaryHeatmap(setL)
    expect_s3_class(p1, class = "ggplot")

    p2 = ssvFeatureEuler(setL, fill_circles = F)
    expect_s3_class(p2, class = "ggplot")
})
