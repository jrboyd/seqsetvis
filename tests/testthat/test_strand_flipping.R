testthat::context("StrandFlipping")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
library(data.table)

qgr = CTCF_in_10a_overlaps_gr[1:5]
qgr = centerFixedSizeGRanges(qgr, 500)

#bed used to intersect bam
# rtracklayer::export.bed(qgr, con = "ctcf_5.bed")
bam_file = system.file("extdata/test.bam", package = "seqsetvis")
xls_file = system.file("extdata/test_peaks.xls", package = "seqsetvis")

test_that("ssvFetchBam strand flipping both, + and -", {
    # skip_on_os("windows")
    test_gr = qgr
    strand(test_gr) = rep("+", 5)
    test_res = ssvFetchBam(bam_file, win_size = 20,
                           fragLens = NA, target_strand = "both",
                           qgr = test_gr, return_data.table = TRUE)

    flip_gr = qgr
    strand(flip_gr) = rep("-", 5)
    flip_res = ssvFetchBam(bam_file, win_size = 20,
                           fragLens = NA, target_strand = "both",
                           qgr = flip_gr, return_data.table = TRUE)

    # ggplot(test_res, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id")
    # ggplot(flip_res, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id")

    expect_equal(test_res[id == 1 & strand == "+"]$y,
                 rev(flip_res[id == 1 & strand == "-"]$y))
    expect_equal(test_res[id == 2 & strand == "+"]$y,
                 rev(flip_res[id == 2 & strand == "-"]$y))
    expect_equal(test_res[id == 3 & strand == "+"]$y,
                 rev(flip_res[id == 3 & strand == "-"]$y))

    expect_equal(test_res[id == 1 & strand == "-"]$y,
                 rev(flip_res[id == 1 & strand == "+"]$y))
    expect_equal(test_res[id == 2 & strand == "-"]$y,
                 rev(flip_res[id == 2 & strand == "+"]$y))
    expect_equal(test_res[id == 3 & strand == "-"]$y,
                 rev(flip_res[id == 3 & strand == "+"]$y))
})

test_that("ssvFetchBam strand flipping both, * and -", {
    # skip_on_os("windows")
    test_gr = qgr
    strand(test_gr) = rep("*", 5)
    test_res = ssvFetchBam(bam_file, win_size = 20,
                           fragLens = NA, target_strand = "both",
                           qgr = test_gr, return_data.table = TRUE)

    flip_gr = qgr
    strand(flip_gr) = rep("-", 5)
    flip_res = ssvFetchBam(bam_file, win_size = 20,
                           fragLens = NA, target_strand = "both",
                           qgr = flip_gr, return_data.table = TRUE)

    expect_equal(test_res[id == 1 & strand == "+"]$y,
                 rev(flip_res[id == 1 & strand == "-"]$y))
    expect_equal(test_res[id == 2 & strand == "+"]$y,
                 rev(flip_res[id == 2 & strand == "-"]$y))
    expect_equal(test_res[id == 3 & strand == "+"]$y,
                 rev(flip_res[id == 3 & strand == "-"]$y))

    expect_equal(test_res[id == 1 & strand == "-"]$y,
                 rev(flip_res[id == 1 & strand == "+"]$y))
    expect_equal(test_res[id == 2 & strand == "-"]$y,
                 rev(flip_res[id == 2 & strand == "+"]$y))
    expect_equal(test_res[id == 3 & strand == "-"]$y,
                 rev(flip_res[id == 3 & strand == "+"]$y))
})

test_that("ssvFetchBam strand flipping +, + and -", {
    # skip_on_os("windows")
    test_gr = qgr
    strand(test_gr) = rep("+", 5)
    test_res = ssvFetchBam(bam_file, win_size = 20,
                           fragLens = NA, target_strand = "+",
                           qgr = test_gr[1], return_data.table = TRUE)

    flip_gr = qgr
    strand(flip_gr) = rep("-", 5)
    flip_res = ssvFetchBam(bam_file, win_size = 20,
                           fragLens = NA, target_strand = "-",
                           qgr = flip_gr[1], return_data.table = TRUE)

    expect_equal(test_res$y, rev(flip_res$y))
})

test_that("ssvFetchBam strand flipping fragLen", {
    #flip_strand precedes fragLen extenstion
    test_gr = qgr
    strand(test_gr) = rep("+", 5)
    test_res = ssvFetchBam(bam_file, win_size = 20,
                           fragLens = 200, target_strand = "both",
                           qgr = test_gr[1], return_data.table = TRUE)

    flip_gr = qgr
    strand(flip_gr) = rep("+", 5)
    flip_res = ssvFetchBam(bam_file, win_size = 20,
                           fragLens = 200, target_strand = "both",
                           qgr = flip_gr[1], return_data.table = TRUE, flip_strand = TRUE)
    test_res$group = "normal"
    flip_res$group = "flip"
    #this is much more obvious in a plot
    # ggplot(rbind(test_res, flip_res), aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("group", ncol = 1)
    test_cor = cor(dcast(test_res, "x~strand", value.var = "y")[,-1])[1,2]
    flip_cor = cor(dcast(flip_res, "x~strand", value.var = "y")[,-1])[1,2]
    expect_true(test_cor - flip_cor > 1)

    expect_false(all(test_res$y == rev(flip_res$y)))
})
