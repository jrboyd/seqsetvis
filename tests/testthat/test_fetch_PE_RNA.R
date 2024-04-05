testthat::context("ssvFetchBamPE.RNA")
# this test both weak and overly complicated
library(seqsetvis)
library(GenomicRanges)
library(data.table)
library(testthat)

library(GenomicRanges)
pkg_dir = system.file(package = "seqsetvis", "extdata", mustWork = TRUE)
bam_files_esr1 = dir(pkg_dir, pattern = "H1.+R1.ESR1_RNA.+bam$", full.names = TRUE)
names(bam_files_esr1) = sub("_R.+", "", basename(bam_files_esr1))
query_gr = GenomicRanges::GRanges("chr6:151656691-152129619:+")
query_gr = GenomicRanges::GRanges("chr6:152121500-152123500:+")

strand(query_gr) = "+"

prof_dt.pe1 = ssvFetchBamPE.RNA(bam_files_esr1, query_gr, unique_names = "PE_noFlip", return_data.table = TRUE, win_size = 1)
prof_dt.pe2 = ssvFetchBamPE.RNA(bam_files_esr1, query_gr, unique_names = "PE_flip", return_data.table = TRUE, win_size = 1, flip_strand = TRUE)

prof_dt1 = ssvFetchBam(
    bam_files_esr1,
    unique_names = "pos_noFlip",
    fragLens = NA,
    query_gr,
    return_data.table = TRUE,
    win_size = 1,
    splice_strategy = "ignore",
    target_strand = "both")
prof_dt2 = ssvFetchBam(
    bam_files_esr1,
    unique_names = "pos_flip",
    fragLens = NA,
    query_gr,
    return_data.table = TRUE,
    win_size = 1,
    splice_strategy = "ignore",
    target_strand = "both",
    flip_strand = TRUE)

strand(query_gr) = "-"

prof_dt3 = ssvFetchBam(
    bam_files_esr1,
    unique_names = "neg_noFlip",
    fragLens = NA,
    query_gr,
    return_data.table = TRUE,
    win_size = 1,
    splice_strategy = "ignore",
    target_strand = "both")
prof_dt4 = ssvFetchBam(
    bam_files_esr1,
    unique_names = "neg_flip",
    fragLens = NA,
    query_gr,
    return_data.table = TRUE,
    win_size = 1,
    splice_strategy = "ignore",
    target_strand = "both",
    flip_strand = TRUE)

prof_dt = rbind(prof_dt1, prof_dt2, prof_dt3, prof_dt4, prof_dt.pe1, prof_dt.pe2)

prof_dt[, gen_x := (start + end) / 2]

ggplot(prof_dt, aes(x = gen_x, y = y, color = strand)) + geom_path() +
    facet_grid(sample~strand)

prof_dt.agg = prof_dt[, .(y = sum(y)), .(x, gen_x, sample)]
ggplot(prof_dt.agg, aes(x = gen_x, y = y)) + geom_path() +
    facet_grid(sample~.)
prof_dt.agg.w = dcast(prof_dt.agg, gen_x~sample, value.var = "y")
prof_dt.agg.w[pos_noFlip > 1]
prof_dt.agg.w[pos_noFlip == PE_noFlip]
test_that("all are equivalent", {
    expect_equal(prof_dt.agg.w$pos_noFlip, prof_dt.agg.w$pos_flip)
    expect_equal(prof_dt.agg.w$pos_noFlip, prof_dt.agg.w$neg_noFlip)
    expect_equal(prof_dt.agg.w$pos_noFlip, prof_dt.agg.w$neg_flip)
    expect_equal(prof_dt.agg.w$pos_noFlip, prof_dt.agg.w$PE_noFlip)
    expect_equal(prof_dt.agg.w$pos_noFlip, prof_dt.agg.w$PE_flip)
})
