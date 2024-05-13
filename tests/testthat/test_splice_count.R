testthat::context("splice_count")

library(seqsetvis)
library(testthat)

#### positive strand query_gr ####
query_gr = GenomicRanges::GRanges("chr6:151656691-152129619:+")
# se_rna_bam_file = system.file(package = "seqsetvis", "extdata/MCF7_R1.ESR1_RNA.bam", mustWork = TRUE)
se_rna_bam_file = system.file(package = "seqsetvis", "extdata/H1_ESC_PE_R1.ESR1_RNA.bam", mustWork = TRUE)


prof_dt.no_flip = ssvFetchBam(
    se_rna_bam_file,
    query_gr,
    return_data.table = TRUE,
    splice_strategy = "splice_count",
    flag = Rsamtools::scanBamFlag(isFirstMateRead = TRUE),
    target_strand = "both"
)



prof_dt.flip = ssvFetchBam(
    se_rna_bam_file,
    query_gr,
    return_data.table = TRUE,
    splice_strategy = "splice_count",
    flag = Rsamtools::scanBamFlag(isFirstMateRead = TRUE),
    target_strand = "both",
    flip_strand = TRUE
)

prof_dt.no_flip
prof_dt.flip

# prof_dt = ssvFetchBam(
#     se_rna_bam_file, query_gr, return_data.table = TRUE, fragLens = NA)
# ggplot(prof_dt, aes(x = (start + end)/2 , y = y)) +
#     geom_path() +
#     ggbio::geom_arch(data = prof_dt.no_flip, aes(x = start, xend = end, height = N, color = strand))


prof_dt.flip.neg = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = TRUE,
    target_strand = "-")
prof_dt.flip.pos = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = TRUE,
    target_strand = "+")
prof_dt.flip.both = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = TRUE,
    target_strand = "both")

prof_dt.no_flip.pos = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = FALSE,
    target_strand = "+")
prof_dt.no_flip.neg = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = FALSE,
    target_strand = "-")

prof_dt.no_flip.both = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = FALSE,
    target_strand = "both")


my_test_strand = function(dt, pos_exp, neg_exp, ns_exp){
    testthat::expect_equal(sum(dt$strand == "+"), pos_exp)
    testthat::expect_equal(sum(dt$strand == "-"), neg_exp)
    testthat::expect_equal(sum(dt$strand == "*"), ns_exp)
}

test_that("splice_count strand - default", {
    my_test_strand(prof_dt.no_flip, 1, 3, 0)
})

test_that("splice_count strand - flip", {
    my_test_strand(prof_dt.flip, 3, 1, 0)
})

test_that("splice_count strand - no flip and target_strand", {
my_test_strand(prof_dt.no_flip.both, 3, 4, 0)
    my_test_strand(prof_dt.no_flip.pos, 3, 0, 0)
    my_test_strand(prof_dt.no_flip.neg, 0, 4, 0)
})

test_that("splice_count strand - flip and target_strand", {
    my_test_strand(prof_dt.flip.both, 4, 3, 0)
    my_test_strand(prof_dt.flip.pos, 4, 0, 0)
    my_test_strand(prof_dt.flip.neg, 0, 3, 0)
})
#### star strand query_gr ####
query_gr = GenomicRanges::GRanges("chr6:151656691-152129619:*")

prof_dt.no_flip.pos = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = FALSE,
    target_strand = "+")
prof_dt.no_flip.neg = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = FALSE,
    target_strand = "-")

prof_dt.no_flip.both = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = FALSE,
    target_strand = "both")


test_that("splice_count strand - no flip and target_strand", {
    my_test_strand(prof_dt.no_flip.both, 3, 4, 0)
    my_test_strand(prof_dt.no_flip.pos, 3, 0, 0)
    my_test_strand(prof_dt.no_flip.neg, 0, 4, 0)
})



#### negative strand query_gr ####
query_gr = GenomicRanges::GRanges("chr6:151656691-152129619:-")

prof_dt.no_flip.pos = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = FALSE,
    target_strand = "+")
prof_dt.no_flip.neg = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = FALSE,
    target_strand = "-")

prof_dt.no_flip.both = ssvFetchBam(
    se_rna_bam_file, query_gr, return_data.table = TRUE, splice_strategy = "splice_count",
    flip_strand = FALSE,
    target_strand = "both")


test_that("splice_count strand - no flip and target_strand", {
    my_test_strand(prof_dt.no_flip.both, 3, 4, 0)
    my_test_strand(prof_dt.no_flip.pos, 3, 0, 0)
    my_test_strand(prof_dt.no_flip.neg, 0, 4, 0)
})


