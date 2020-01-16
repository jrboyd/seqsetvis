testthat::context("force_skip_centerFix")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
library(data.table)

qgr = CTCF_in_10a_overlaps_gr[1:5]
k = unique(as.character(seqnames(qgr)))
seqlevels(qgr) = k
qgr = centerFixedSizeGRanges(qgr, 500)

exp_factor = 100
qgr = resize(qgr, seq_along(qgr)*exp_factor, fix = "center")
#bed used to intersect bam
# rtracklayer::export.bed(qgr, con = "ctcf_5.bed")
bam_file = system.file("extdata/test.bam", package = "seqsetvis", mustWork = TRUE)
bw_file = system.file("extdata/test_loading.bw", package = "seqsetvis", mustWork = TRUE)
win_size = 50


test_that("ssvFetchBam - force_skip_centerFix", {
    bam_dt = ssvFetchBam(bam_file, qgr,
                         win_size = win_size,
                         force_skip_centerFix = TRUE,
                         return_data.table = TRUE, fragLens = 200)
    size_dt = bam_dt[, .(region_width = diff(range(x)) + win_size), .(id)]

    is_exp = size_dt[, region_width / as.numeric(id)] == exp_factor
    expect_true(all(is_exp))
})

test_that("ssvFetchBigwig - force_skip_centerFix", {
    skip_on_os("windows")
    pos = c(20, 180, 210, 440, 520, 521)
    region_size = 30
    test_qgr = GRanges("chrTest", IRanges(pos+1, pos + region_size))
    test_qgr = resize(test_qgr, seq_along(test_qgr)*region_size, fix = "center")
    names(test_qgr) = seq_along(test_qgr)
    bw_dt = ssvFetchBigwig(bw_file, test_qgr,
                         win_size = 30,
                         force_skip_centerFix = TRUE,
                         return_data.table = TRUE)
    size_dt = bw_dt[, .(region_width = diff(range(x)) + 30), .(id)]

    is_exp = size_dt[, region_width / as.numeric(id)] == region_size
    expect_true(all(is_exp))
})

test_that("ssvFetchBamPE - force_skip_centerFix", {
    pe_file = system.file("extdata/Bcell_PE.mm10.bam", package = "seqsetvis", mustWork = TRUE)
    data("Bcell_peaks")
    qgrPE = Bcell_peaks
    qgrPE$name = seq_along(qgrPE)
    qgrPE = resize(qgrPE, seq_along(qgrPE)*exp_factor, fix = "center")
    bam_dt = ssvFetchBamPE(pe_file, qgrPE,
                         win_size = win_size,
                         force_skip_centerFix = TRUE,
                         return_data.table = TRUE)
    size_dt = bam_dt[, .(region_width = diff(range(x)) + win_size), .(id)]
    is_exp = size_dt[, region_width / as.numeric(id)] == exp_factor
    expect_true(all(is_exp))

})

test_that("ssvFetchGranges - force_skip_centerFix", {
    qgrGR = CTCF_in_10a_overlaps_gr[1:5]
    qgrGR = resize(qgrGR, seq_along(qgrGR)*exp_factor, fix = "center")
    gr_dt = ssvFetchGRanges(CTCF_in_10a_overlaps_gr[1:5], qgrGR,
                           win_size = win_size,
                           force_skip_centerFix = TRUE,
                           return_data.table = TRUE)
    size_dt = gr_dt[, .(region_width = diff(range(x)) + win_size), .(id)]
    is_exp = size_dt[, region_width / as.numeric(id)] == exp_factor
    expect_true(all(is_exp))
})
