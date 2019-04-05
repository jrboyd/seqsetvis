testthat::context("MakeMembTable")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
set_a = 1:8
set_b = 5:9
set_c = 8:10

gr_a = GRanges("chr1", IRanges(set_a*10, set_a*10+3))
gr_b = GRanges("chr1", IRanges(set_b*10-2, set_b*10+3))
gr_c = GRanges("chr1", IRanges(set_c*10+-4, set_c*10+1))

expected_table = as.matrix(data.frame("set_A" = 1:10 %in% set_a,
                                      "set_B" = 1:10 %in% set_b,
                                      "set_C" = 1:10 %in% set_c, row.names = 1:10))

expected_table_named = as.matrix(data.frame("named_set_A" = 1:10 %in% set_a,
                                            "named_set_B" = 1:10 %in% set_b,
                                            "named_set_C" = 1:10 %in% set_c, row.names = 1:10))

test_that("ssvMakeMembTable unnamed numeric set list", {
    sets_list = list(set_a, set_b, set_c)
    mt = ssvMakeMembTable(sets_list)
    expect_is(mt, "matrix")
    expect_equal(as.logical(mt), as.logical(expected_table))
    expect_equal(rownames(mt), rownames(expected_table))
    expect_equal(colnames(mt), colnames(expected_table))
})

test_that("ssvMakeMembTable unnamed character set list", {
    sets_list = list(set_a, set_b, set_c)
    sets_list = lapply(sets_list, function(x)letters[x])
    mt = ssvMakeMembTable(sets_list)
    expect_is(mt, "matrix")
    expect_equal(as.logical(mt), as.logical(expected_table))
    expect_equal(rownames(mt), letters[1:nrow(expected_table)])
    expect_equal(colnames(mt), colnames(expected_table))
})

test_that("ssvMakeMembTable unnamed matrix", {
    sets_list = list(set_a, set_b, set_c)
    memb_table = ssvMakeMembTable(sets_list)
    colnames(memb_table) = NULL
    mt = ssvMakeMembTable(memb_table)
    expect_is(mt, "matrix")
    expect_equal(as.logical(mt), as.logical(expected_table))
    expect_equal(rownames(mt), rownames(expected_table))
    expect_equal(colnames(mt), colnames(expected_table))
})

test_that("ssvMakeMembTable unnamed data.frame", {
    sets_list = list(set_a, set_b, set_c)
    memb_table = ssvMakeMembTable(sets_list)
    colnames(memb_table) = NULL
    memb_table = as.data.frame(memb_table)
    mt = ssvMakeMembTable(memb_table)
    expect_is(mt, "matrix")
    expect_equal(as.logical(mt), as.logical(expected_table))
    expect_equal(rownames(mt), rownames(expected_table))
    expect_equal(colnames(mt), colnames(expected_table))
})

test_that("ssvMakeMembTable unnamed list of GRanges", {
    sets_gr_list = list(gr_a, gr_b, gr_c)

    suppressMessages({mt = ssvMakeMembTable(sets_gr_list)})
    expect_is(mt, "matrix")
    expect_equal(as.logical(mt), as.logical(expected_table))
    expect_equal(rownames(mt), rownames(expected_table))
    expect_equal(colnames(mt), colnames(expected_table))
})

test_that("ssvMakeMembTable NAMED list of GRanges", {
    sets_gr_list = list("named_set_A" = gr_a, "named_set_B" = gr_b, "named_set_C" = gr_c)

    suppressMessages({mt = ssvMakeMembTable(sets_gr_list)})
    expect_is(mt, "matrix")
    expect_equal(as.logical(mt), as.logical(expected_table_named))
    expect_equal(rownames(mt), rownames(expected_table_named))
    expect_equal(colnames(mt), colnames(expected_table_named))
})

test_that("ssvMakeMembTable on mixed type list causes error", {
    expect_error(ssvMakeMembTable(list(1, "a")))
})
