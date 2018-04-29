#tests seqsetvis::fetchWindowedBigwig and seqsetvis::fetchWindowedBigwigList
#originally data.table version of each which may explain some strangeness
library(seqsetvis)
library(testthat)
library(GenomicRanges)
# library(rtracklayer)
test_bw = system.file("extdata/test_loading.bw", package = "seqsetvis", mustWork = TRUE)
pos = c(20, 180, 210, 440, 520, 521)
region_size = 30
test_qgr = GRanges("chrTest", IRanges(pos+1, pos + region_size))
exp_colnames = c("seqnames", "start", "end", "width", "strand", "id", "y", "x")[6:8]


test_that("fetchWindowedBigwig return expected for valid even win_size", {
    skip_on_os("windows")
    #these should all work cleanly
    for(win in c(2, 6, 10, 30)){
        bw_gr = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = test_qgr)
        expect_is(bw_gr, "GRanges")
        expect_equal(colnames(mcols(bw_gr)), exp_colnames)
        for(tid in unique(bw_gr$id)){
            test_gr = bw_gr[bw_gr$id == tid]
            expect_equal(length(test_gr), region_size / win) #expected number of regions
            expect_equal(length(findOverlaps(test_gr, test_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("fetchWindowedBigwig return expected for valid odd win_size", {
    skip_on_os("windows")
    for(win in c(1, 3, 5, 15)){
        bw_gr = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = test_qgr)
        expect_is(bw_gr, "GRanges")
        expect_equal(colnames(mcols(bw_gr)), exp_colnames)
        for(tid in unique(bw_gr$id)){
            test_gr = bw_gr[bw_gr$id == tid]
            expect_equal(length(test_gr), region_size / win) #expected number of regions
            expect_equal(length(findOverlaps(test_gr, test_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("fetchWindowedBigwig use GRanges names as id", {
    skip_on_os("windows")
    for(win in c(1, 3, 5, 15)){
        qgr = test_qgr
        names(qgr) = paste0("myNames_", seq_along(qgr))
        bw_gr = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = qgr)
        expect_is(bw_gr, "GRanges")
        expect_equal(colnames(mcols(bw_gr)), exp_colnames)
        expect_true(all(grepl("myNames", bw_gr$id)))
        for(tid in unique(bw_gr$id)){
            test_gr = bw_gr[bw_gr$id == tid]
            expect_equal(length(test_gr), region_size / win) #expected number of regions
            expect_equal(length(findOverlaps(test_gr, test_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("fetchWindowedBigwig patches missing values", {
    skip_on_os("windows")
    for(win in c(1, 5, 20)){
        qgr = GRanges("chrTest", IRanges(1, 2000))
        bw_gr = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = qgr)
        expect_is(bw_gr, "GRanges")
        expect_equal(colnames(mcols(bw_gr)), exp_colnames)
    }
})

test_that("fetchWindowedBigwig throws message if widths aren't divisble by win_size", {
    skip_on_os("windows")
    mix_width_gr = test_qgr
    end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)
    for(win in c(7, 31)){
        expect_message(regexp = "widths of qgr were not identical and evenly divisible by win_size", {
            fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = mix_width_gr)
        })
    }
})



test_that("fetchWindowedBigwig throws warning if widths vary", {
    skip_on_os("windows")
    mix_width_gr = test_qgr
    end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)*3
    for(win in c(1, 3)){
        expect_message(regexp = "widths of qgr were not identical and evenly divisible by win_size", {
            fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = mix_width_gr)
        })
    }
})

test_that("fetchWindowedBigwigList works with character vector bw_files", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = fetchWindowedBigwigList(file_paths = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr)})
    expect_is(res, "GRanges")
    expect_equal(colnames(mcols(res)), c(exp_colnames, "sample"))
})

test_that("fetchWindowedBigwigList works with list bw_files", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    bw_files = as.list(bw_files)
    hidden = capture_output({res = fetchWindowedBigwigList(file_paths = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr)})
    expect_is(res, "GRanges")
    expect_equal(colnames(mcols(res)), c(exp_colnames, "sample"))
})

test_that("fetchWindowedBigwigList can set variable name", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = fetchWindowedBigwigList(file_paths = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr,
                                                           names_variable = "group")})
    expect_is(res, "GRanges")
    expect_equal(colnames(mcols(res)), c(exp_colnames, "group"))
})

test_that("fetchWindowedBigwigList duplicate names throws error", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    expect_error(
        fetchWindowedBigwigList(file_paths = bw_files,
                                win_size = 3,
                                qgr = test_qgr)
    )
})

