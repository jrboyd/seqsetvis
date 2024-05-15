testthat::context("ssvAnnotateSubjectGRanges")
library(seqsetvis)
library(GenomicRanges)
library(testthat)

data(CTCF_in_10a_narrowPeak_grs)
np_grs = CTCF_in_10a_narrowPeak_grs
olap_gr = ssvOverlapIntervalSets(np_grs)

# get_exp = function(x){
#     message('c(\n"', paste(x, collapse = '", \n"'), '"\n)'    )
# }
# get_exp(colnames(mcols(test_gr)))

test_that("ssvAnnotateSubjectGRanges - single GR no names", {
    test_gr = ssvAnnotateSubjectGRanges(np_grs$MCF10A_CTCF, olap_gr)
    expect_equal(colnames(mcols(test_gr)),
                 c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1_CTCF",
                   "name", "score", "signalValue", "pValue", "qValue", "peak"))
})

test_that("ssvAnnotateSubjectGRanges - single GR with name", {
    test_gr = ssvAnnotateSubjectGRanges(np_grs$MCF10A_CTCF, olap_gr, annotation_name = "MCF10A")
    expect_equal(colnames(mcols(test_gr)),
                 c("MCF10A_CTCF", "MCF10AT1_CTCF",
                   "MCF10CA1_CTCF", "name.MCF10A",
                   "score.MCF10A", "signalValue.MCF10A",
                   "pValue.MCF10A", "qValue.MCF10A", "peak.MCF10A"))
})

test_that("ssvAnnotateSubjectGRanges - list of GRs", {
    test_gr = ssvAnnotateSubjectGRanges(np_grs, olap_gr)
    expect_equal(colnames(mcols(test_gr)),
                 c(
                     "MCF10A_CTCF",
                     "MCF10AT1_CTCF",
                     "MCF10CA1_CTCF",
                     "name.MCF10A_CTCF",
                     "score.MCF10A_CTCF",
                     "signalValue.MCF10A_CTCF",
                     "pValue.MCF10A_CTCF",
                     "qValue.MCF10A_CTCF",
                     "peak.MCF10A_CTCF",
                     "name.MCF10AT1_CTCF",
                     "score.MCF10AT1_CTCF",
                     "signalValue.MCF10AT1_CTCF",
                     "pValue.MCF10AT1_CTCF",
                     "qValue.MCF10AT1_CTCF",
                     "peak.MCF10AT1_CTCF",
                     "name.MCF10CA1_CTCF",
                     "score.MCF10CA1_CTCF",
                     "signalValue.MCF10CA1_CTCF",
                     "pValue.MCF10CA1_CTCF",
                     "qValue.MCF10CA1_CTCF",
                     "peak.MCF10CA1_CTCF"
                 ))
})

test_that("ssvAnnotateSubjectGRanges - list of GRs - override names", {
    test_gr = ssvAnnotateSubjectGRanges(np_grs, olap_gr, LETTERS[1:3])
    expect_equal(colnames(mcols(test_gr)),
                 c(
                     "MCF10A_CTCF",
                     "MCF10AT1_CTCF",
                     "MCF10CA1_CTCF",
                     "name.A",
                     "score.A",
                     "signalValue.A",
                     "pValue.A",
                     "qValue.A",
                     "peak.A",
                     "name.B",
                     "score.B",
                     "signalValue.B",
                     "pValue.B",
                     "qValue.B",
                     "peak.B",
                     "name.C",
                     "score.C",
                     "signalValue.C",
                     "pValue.C",
                     "qValue.C",
                     "peak.C"
                 ))
})

test_that("ssvAnnotateSubjectGRanges - GRangesList - override names", {
    test_gr = ssvAnnotateSubjectGRanges(GRangesList(np_grs), olap_gr, LETTERS[1:3])
    expect_equal(colnames(mcols(test_gr)),
                 c(
                     "MCF10A_CTCF",
                     "MCF10AT1_CTCF",
                     "MCF10CA1_CTCF",
                     "name.A",
                     "score.A",
                     "signalValue.A",
                     "pValue.A",
                     "qValue.A",
                     "peak.A",
                     "name.B",
                     "score.B",
                     "signalValue.B",
                     "pValue.B",
                     "qValue.B",
                     "peak.B",
                     "name.C",
                     "score.C",
                     "signalValue.C",
                     "pValue.C",
                     "qValue.C",
                     "peak.C"
                 ))
})

test_that("ssvAnnotateSubjectGRanges - list of GRs - error no names", {
    bad_grs = np_grs
    names(bad_grs) = NULL
    expect_error(ssvAnnotateSubjectGRanges(bad_grs, olap_gr), "annotation_source must have names set.")
})

test_that("ssvOverlapIntervalSets and ssvConsensusIntervalSets - message", {
    #when np_grs overlap multi-multi it's not clear which metadata should be preserved.
    expect_message({
        olap_anno_gr = ssvOverlapIntervalSets(np_grs, preserve_mcols = TRUE, ext = 50000000)
    },
    "Resolving multiple overlapping annotatation_source items with default function"
    )

    expect_message({
        cons_anno_gr = ssvConsensusIntervalSets(np_grs, preserve_mcols = TRUE, ext = 500000)
    },
    "Resolving multiple overlapping annotatation_source items with default function"
    )

    expect_message({
        anno_gr = ssvAnnotateSubjectGRanges(annotation_source = np_grs$MCF10A_CTCF, subject_gr = olap_anno_gr[1])
    },
    "Resolving multiple overlapping annotatation_source items with default function"
    )


})

test_that("ssvAnnotateSubjectGranges - multiple overlapping annotatation_source", {
    test_multi_resolver_FUN = function(x, variable.name = NULL){
        if(is.character(x)){
            x = "a"
        }else if(is.numeric(x)){
            x = 1
        }else if(is.factor(x)){
            x = factor("a")
        }else if(is.logical(x)){
            x = TRUE
        }
        x
    }

    olap_anno_gr = ssvOverlapIntervalSets(np_grs, preserve_mcols = TRUE, ext = 50000000)
    mcols(olap_anno_gr) = NULL

    anno_multi = ssvAnnotateSubjectGRanges(annotation_source = np_grs, subject_gr = olap_anno_gr[1])
    expect_equal(length(anno_multi), 1)
    expect_equal(anno_multi$score.MCF10A_CTCF, 3683)
    expect_equal(anno_multi$name.MCF10A_CTCF, "MCF10A_CTCF_pooled_peak_1190-from_11")

    anno_test = ssvAnnotateSubjectGRanges(annotation_source = np_grs, subject_gr = olap_anno_gr[1], multi_resolver_FUN = test_multi_resolver_FUN)
    expect_equal(length(anno_test), 1)
    expect_equal(anno_test$score.MCF10A_CTCF, 1)
    expect_equal(anno_test$name.MCF10A_CTCF, "a")

})



