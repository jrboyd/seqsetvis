library(seqsetvis)
library(testthat)
library(GenomicRanges)
a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1))
b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1))
c = GRanges("chr1", IRanges(8:10*10+5, 8:10*10+6))

test_that("setPlotVenn various set sizes, no outside counts.", {
    setL = list("a" = a, "b" = b, "c" = c)
    p1 = setPlotVenn(setL[1])
    p1
    expect_s3_class(p1, class = "ggplot")

    p2 = setPlotVenn(setL[1:2])
    p2
    expect_s3_class(p2, class = "ggplot")

    p3 = setPlotVenn(setL[1:3])
    p3
    expect_s3_class(p3, class = "ggplot")
})

test_that("setPlotVenn various set sizes, with outside counts.", {
    setL = list("a" = a, "b" = b, "c" = c)
    p1 = setPlotVenn(setL[1], show_outside_count = T)
    p1
    expect_s3_class(p1, class = "ggplot")

    p2 = setPlotVenn(setL[1:2], show_outside_count = T)
    p2
    expect_s3_class(p2, class = "ggplot")

    p3 = setPlotVenn(setL[1:3], show_outside_count = T)
    p3
    expect_s3_class(p3, class = "ggplot")
})

test_that("setPlotVenn can't plot more than 3 sets.", {
    setL = list("a" = a, "b" = b, "c" = c, "d" = c, "e" = a)
    expect_error(setPlotVenn(setL[1:4]))
    expect_error(setPlotVenn(setL[1:5]))
})

test_that("setPlotVenn circle colors", {
    setL = list("a" = a, "b" = b, "c" = c)
    p1 = setPlotVenn(setL, circle_color = "blue")
    expect_s3_class(p1, class = "ggplot")
})
