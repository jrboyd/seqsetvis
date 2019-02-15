#tests seqsetvis::ssvFetchBigwig.single and seqsetvis::ssvFetchBigwig
#originally data.table version of each which may explain some strangeness
library(seqsetvis)
library(testthat)
library(GenomicRanges)
# library(rtracklayer)
test_bw = system.file("extdata/test_loading.bw", package = "seqsetvis", mustWork = TRUE)
pos = c(20, 180, 210, 440, 520, 521)
region_size = 30
test_gr = GRanges("chrTest", IRanges(pos+1, pos + region_size))
exp_colnames = c("seqnames", "start", "end", "width", "strand", "id", "y", "x")[6:8]


test_that("ssvFetchBigwig.single return expected for valid even win_size", {
    skip_on_os("windows")
    #these should all work cleanly
    for(win in c(2, 6, 10, 30)){
        bw_gr = ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = test_gr)
        expect_is(bw_gr, "GRanges")
        expect_equal(colnames(mcols(bw_gr)), exp_colnames)
        for(tid in unique(bw_gr$id)){
            tgr = bw_gr[bw_gr$id == tid]
            expect_equal(length(tgr), region_size / win) #expected number of regions
            expect_equal(length(findOverlaps(tgr, tgr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("ssvFetchBigwig.single return expected for valid odd win_size", {
    skip_on_os("windows")
    for(win in c(1, 3, 5, 15)){
        bw_gr = ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = test_gr)
        expect_is(bw_gr, "GRanges")
        expect_equal(colnames(mcols(bw_gr)), exp_colnames)
        for(tid in unique(bw_gr$id)){
            tgr = bw_gr[bw_gr$id == tid]
            expect_equal(length(tgr), region_size / win) #expected number of regions
            expect_equal(length(findOverlaps(tgr, tgr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("ssvFetchBigwig.single use GRanges names as id", {
    skip_on_os("windows")
    for(win in c(1, 3, 5, 15)){
        qgr = test_gr
        names(qgr) = paste0("myNames_", seq_along(qgr))
        bw_gr = ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = qgr)
        expect_is(bw_gr, "GRanges")
        expect_equal(colnames(mcols(bw_gr)), exp_colnames)
        expect_true(all(grepl("myNames", bw_gr$id)))
        for(tid in unique(bw_gr$id)){
            tgr = bw_gr[bw_gr$id == tid]
            expect_equal(length(tgr), region_size / win) #expected number of regions
            expect_equal(length(findOverlaps(tgr, tgr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("ssvFetchBigwig.single patches missing values", {
    skip_on_os("windows")
    for(win in c(1, 5, 20)){
        qgr = GRanges("chrTest", IRanges(1, 2000))
        bw_gr = ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = qgr)
        expect_is(bw_gr, "GRanges")
        expect_equal(colnames(mcols(bw_gr)), exp_colnames)
    }
})

test_that("ssvFetchBigwig.single throws message if widths aren't divisble by win_size", {
    skip_on_os("windows")
    mix_width_gr = test_gr
    end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)
    for(win in c(7, 31)){
        expect_message(regexp = "widths of qgr were not identical and evenly divisible by win_size", {
            ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = mix_width_gr)
        })
    }
})



test_that("ssvFetchBigwig.single throws warning if widths vary", {
    skip_on_os("windows")
    mix_width_gr = test_gr
    end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)*3
    for(win in c(1, 3)){
        expect_message(regexp = "widths of qgr were not identical and evenly divisible by win_size", {
            ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = mix_width_gr)
        })
    }
})

test_that("ssvFetchBigwig query GRanges output id set", {
    skip_on_os("windows")
    tgr = test_gr
    bw_files = test_bw
    gr_sample = ssvFetchBigwig(bw_files, win_size = 5, qgr = test_gr, return_data.table = TRUE)
    expect_true("id" %in% colnames(gr_sample))

    names(tgr) = NULL
    gr_sample = ssvFetchBigwig(bw_files, win_size = 5, qgr = test_gr, return_data.table = TRUE)
    expect_true("id" %in% colnames(gr_sample))

    tgr$id = seq_along(tgr) #when id set, id in output is missing
    gr_sample = ssvFetchBigwig(bw_files, win_size = 5, qgr = test_gr, return_data.table = TRUE)
    expect_true("id" %in% colnames(gr_sample))

    tgr = test_gr
    tgr$id = seq_along(tgr) #when id set, id in output is missing
    gr_sample = ssvFetchBigwig(bw_files, win_size = 5, qgr = test_gr, return_data.table = TRUE)
    expect_true("id" %in% colnames(gr_sample))

})

test_that("ssvFetchBigwig works with character vector bw_files", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = ssvFetchBigwig(file_paths = bw_files,
                                                           win_size = 3,
                                                           qgr = test_gr)})
    expect_is(res, "GRanges")
    expect_equal(colnames(mcols(res)), c(exp_colnames, "sample"))
})

test_that("ssvFetchBigwig works with list bw_files", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    bw_files = as.list(bw_files)
    hidden = capture_output({res = ssvFetchBigwig(file_paths = bw_files,
                                                           win_size = 3,
                                                           qgr = test_gr)})
    expect_is(res, "GRanges")
    expect_equal(colnames(mcols(res)), c(exp_colnames, "sample"))
})

test_that("ssvFetchBigwig can set variable name", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = ssvFetchBigwig(file_paths = bw_files,
                                                           win_size = 3,
                                                           qgr = test_gr,
                                                           names_variable = "group")})
    expect_is(res, "GRanges")
    expect_equal(colnames(mcols(res)), c(exp_colnames, "group"))
})

test_that("ssvFetchBigwig duplicate names throws error", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    expect_error(
        ssvFetchBigwig(file_paths = bw_files,
                                win_size = 3,
                                qgr = test_gr)
    )
})

test_that("ssvFetchBigwig sample method correct bins", {
    skip_on_os("windows")
    test_gr2 = test_gr[c(1,2,5)]
    gr_sample = ssvFetchBigwig(test_bw, win_size = 3, qgr = test_gr2)
    #for non-overlapping, all reduce granges width should be equal to input qgr
    expect_true(all(width(reduce(gr_sample)) == width(test_gr2)))
    #for non-overlapping, all granges width should be equal to win_size
    expect_true(all(width(gr_sample) == 3))

    #same as above but width are auto adjusted to be divisible by win_size
    gr_sample = ssvFetchBigwig(test_bw, win_size = 4, qgr = test_gr2)
    expect_true(all(width(reduce(gr_sample)) == 32))
    expect_true(all(width(gr_sample) == 4))
})

test_that("ssvFetchBigwig summary method correct bins", {
    skip_on_os("windows")
    test_gr2 = test_gr[c(1,2,5)]
    gr_sample = ssvFetchBigwig(test_bw, win_size = 3, qgr = test_gr2, win_method = "summary")
    #for non-overlapping, all reduce granges width should be equal to input qgr
    expect_true(all(width(reduce(gr_sample)) == width(test_gr2)))
    #for non-overlapping, all granges width should be equal to win_size
    expect_true(all(width(gr_sample) == 30 / 3))

    #same as above but width are auto adjusted to be divisible by win_size
    gr_sample = ssvFetchBigwig(test_bw, win_size = 4, qgr = test_gr2, win_method = "summary")
    expect_true(all(width(reduce(gr_sample)) == 30))
    expect_true(all(width(gr_sample) %in% 7:8))
})
