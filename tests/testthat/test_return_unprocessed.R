testthat::context("return_unprocessed")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
library(data.table)

qgr = CTCF_in_10a_overlaps_gr#[1:5]
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


test_that("ssvFetchBamPE - return_unprocessed", {
    pe_file = system.file("extdata/Bcell_PE.mm10.bam", package = "seqsetvis", mustWork = TRUE)
    data("Bcell_peaks")
    qgrPE = Bcell_peaks
    qgrPE$name = seq_along(qgrPE)
    qgrPE = resize(qgrPE, seq_along(qgrPE)*exp_factor, fix = "center")
    raw_dt = ssvFetchBamPE(pe_file, qgrPE, return_unprocessed = TRUE)
    expect_s3_class(raw_dt, "data.table")
    expect_setequal(colnames(raw_dt),
                    c("which_label", "seqnames", "qname", "strand", "start", "width", "cigar", "isize", "id", "sample"))
})

test_that("ssvFetchBam - return_unprocessed", {
    raw_dt = ssvFetchBam(bam_file, qgr,
                         win_size = win_size,
                         force_skip_centerFix = TRUE,
                         return_data.table = TRUE,
                         return_unprocessed = TRUE)
    expect_s3_class(raw_dt, "data.table")
    expect_setequal(colnames(raw_dt),
                    c("which_label", 'seqnames', "strand", "start", "width",
                    "cigar", "read_id", "flag", "mapq", "mrnm", "mpos", "isize",
                    "seq", "qual", "id", "sample"))
})

test_that("ssvFetchBamPE - expandCigar", {
    pe_file = system.file("extdata/Bcell_PE.mm10.bam", package = "seqsetvis", mustWork = TRUE)
    data("Bcell_peaks")
    qgrPE = Bcell_peaks
    qgrPE$name = seq_along(qgrPE)
    qgrPE = resize(qgrPE, seq_along(qgrPE)*exp_factor, fix = "center")
    raw_dt = ssvFetchBamPE(pe_file, qgrPE, return_unprocessed = TRUE)
    exp_gr = expandCigar(raw_dt)
    expect_s4_class(exp_gr, "GRanges")
    expect_setequal(colnames(mcols(exp_gr)), c("which_label", "cigar_type"))
    exp_dt = expandCigar(raw_dt, return_data.table = TRUE)
    expect_s3_class(exp_dt, "data.table")
    expect_setequal(colnames(exp_dt), c("which_label", "seqnames", "strand", "start", "end", "cigar_type"))
})

test_that("ssvFetchBam - expandCigar", {
    raw_dt = ssvFetchBam(bam_file, qgr,
                         win_size = win_size,
                         force_skip_centerFix = TRUE,
                         return_data.table = TRUE,
                         return_unprocessed = TRUE)
    exp_gr = expandCigar(raw_dt)
    expect_s4_class(exp_gr, "GRanges")
    expect_setequal(colnames(mcols(exp_gr)), c("which_label", "cigar_type"))
    exp_dt = expandCigar(raw_dt, return_data.table = TRUE)
    expect_s3_class(exp_dt, "data.table")
    expect_setequal(colnames(exp_dt), c("which_label", "seqnames", "strand", "start", "end", "cigar_type"))
})
