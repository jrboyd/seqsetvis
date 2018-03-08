library(seqsetvis)
library(testthat)
library(GenomicRanges)
gr_a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1))
gr_b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1))
gr_c = GRanges("chr1", IRanges(8:10*10+5, 8:10*10+6))
setL = list("gr_a" = gr_a, "gr_b" = gr_b, "gr_c" = gr_c)

test_that("setPlotVenn various set sizes, no outside counts. 1.", {
    p1 = setPlotVenn(setL[1])
    p1
    expect_s3_class(p1, class = "ggplot")
})

test_that("setPlotVenn various set sizes, no outside counts. 2.", {
    p2 = setPlotVenn(setL[1:2])
    p2
    expect_s3_class(p2, class = "ggplot")
})

test_that("setPlotVenn various set sizes, no outside counts. 3.", {
    p3 = setPlotVenn(setL[1:3])
    p3
    expect_s3_class(p3, class = "ggplot")
})

test_that("setPlotVenn various set sizes, with outside counts. 1.", {
    p1 = setPlotVenn(setL[1], show_outside_count = T)
    p1
    expect_s3_class(p1, class = "ggplot")
})

test_that("setPlotVenn various set sizes, with outside counts. 2.", {
    p2 = setPlotVenn(setL[1:2], show_outside_count = T)
    p2
    expect_s3_class(p2, class = "ggplot")
})

test_that("setPlotVenn various set sizes, with outside counts. 3.", {
    p3 = setPlotVenn(setL[1:3], show_outside_count = T)
    p3
    expect_s3_class(p3, class = "ggplot")
})

test_that("setPlotVenn can't plot more than 3 sets.", {
    setL = list("gr_a" = gr_a, "gr_b" = gr_b, "gr_c" = gr_c, "gr_d" = gr_c, "gr_e" = gr_a)
    expect_error(setPlotVenn(setL[1:4]))
    expect_error(setPlotVenn(setL[1:5]))
})

test_that("setPlotVenn circle colors", {
    p1 = setPlotVenn(setL, circle_color = "blue")
    expect_s3_class(p1, class = "ggplot")
})
