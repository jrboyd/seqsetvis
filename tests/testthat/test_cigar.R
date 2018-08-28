library(data.table)
library(testthat)
library(GenomicRanges)



bam_dt = fread(system.file("extdata/test_cigar_dt.csv", package = "seqsetvis"))
bam_dt[, start := start - min(start) + 1]
bam_dt = bam_dt[order(start)]
expand_cigar_dt = seqsetvis:::.expand_cigar_dt
test_that("expand_cigar_dt simple alignment", {
    bdt = bam_dt[cigar == "82M"]
    exp_dt = expand_cigar_dt(bdt)
    expect_equal(nrow(exp_dt), 1)
    expect_equal(exp_dt$end, end(IRanges(start = bdt$start, width = bdt$width)))
    expect_equal(exp_dt$start, bdt$start)
    expect_equal(exp_dt$cigar_type, "M")
})

test_that("expand_cigar_dt deletion spanning alignment", {
    bdt = bam_dt[cigar == "46M2D45M",]
    exp_dt = expand_cigar_dt(bdt)
    expect_equal(nrow(exp_dt), 3)
    expect_equal(exp_dt$end, c(46, 46+2, 46+2+45))
    expect_equal(exp_dt$start, c(1, 1+46, 1+46+2))
    expect_equal(exp_dt$cigar_type, c("M", "D", "M"))
})

test_that("expand_cigar_dt skip region alignment", {
    bdt = bam_dt[cigar == "63M5005N28M"]
    exp_dt = expand_cigar_dt(bdt)
    expect_equal(nrow(exp_dt), 2)
    exp_dt = expand_cigar_dt(bdt, op_2count = c("M", "D", "=", "X", "N"))
    expect_equal(nrow(exp_dt), 3)
    expect_equal(exp_dt$end, c(21179, 21179 + 5005, 21179 + 5005+28))
    expect_equal(exp_dt$start, c(21117, 21117+63, 21117+63+5005))
    expect_equal(exp_dt$cigar_type, c("M", "N", "M"))
})

test_that("expand_cigar_dt skip region multiple", {
    #there happen to be 6 reads spanning the same gap
    bdt = bam_dt[grepl("6397N", cigar)]
    exp_dt = expand_cigar_dt(bdt)
    expect_equal(nrow(exp_dt), 12)
    exp_dt = expand_cigar_dt(bdt, op_2count = c("M", "D", "=", "X", "N"))
    expect_equal(nrow(exp_dt), 18)
    expect_equal(sum(exp_dt$end == 723), 6)
    expect_equal(sum(exp_dt$end == 7120), 6)
    expect_equal(sum(exp_dt$start == 724), 6)
    expect_equal(sum(exp_dt$start == 7121), 6)
    expect_equal(sum(exp_dt$cigar_type == "M"), 12)
    expect_equal(sum(exp_dt$cigar_type == "N"), 6)
})

