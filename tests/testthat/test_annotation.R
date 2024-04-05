testthat::context("ssvAnnotateSubjectGRanges")
library(seqsetvis)
library(GenomicRanges)
library(testthat)

np_grs = CTCF_in_10a_narrowPeak_grs
olap_gr = ssvOverlapIntervalSets(np_grs)

ssvAnnotateSubjectGRanges(np_grs$MCF10A_CTCF, olap_gr, annotation_name = "MCF10A")
ssvAnnotateSubjectGRanges(np_grs$MCF10A_CTCF, olap_gr)
ssvAnnotateSubjectGRanges(np_grs, olap_gr)
ssvAnnotateSubjectGRanges(np_grs, olap_gr, LETTERS[1:3])
ssvAnnotateSubjectGRanges(GRangesList(np_grs), olap_gr, LETTERS[1:3])

test_that("ssvAnnotateSubjectGRanges", {
})
