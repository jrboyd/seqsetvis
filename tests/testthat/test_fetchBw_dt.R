testthat::context("FetchBw_dt")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
# library(rtracklayer)
test_bw = system.file("extdata/test_loading.bw", package = "seqsetvis", mustWork = TRUE)
pos = c(20, 180, 210, 440, 520, 521)
region_size = 30
test_qgr = GRanges("chrTest", IRanges(pos+1, pos + region_size))
exp_colnames = c("seqnames", "start", "end", "width", "strand", "id", "y", "x")

ssvFetchBigwig.single = seqsetvis:::ssvFetchBigwig.single

test_that("ssvFetchBigwig.single return expected for valid even win_size", {
    skip_on_os("windows")
    #these should all work cleanly
    for(win in c(2, 6, 10, 30)){
        bw_dt = ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = test_qgr, return_data.table = TRUE)
        expect_is(bw_dt, "data.table")
        expect_equal(colnames(bw_dt), exp_colnames)
        for(tid in unique(bw_dt$id)){
            test_dt = bw_dt[id == tid]
            expect_equal(nrow(test_dt), region_size / win) #expected number of regions
            bw_gr = GRanges(test_dt)
            expect_equal(length(findOverlaps(bw_gr, bw_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

test_that("ssvFetchBigwig.single return expected for valid odd win_size", {
    skip_on_os("windows")
    for(win in c(1, 3, 5, 15)){
        bw_dt = ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = test_qgr, return_data.table = TRUE)
        expect_is(bw_dt, "data.table")
        expect_equal(colnames(bw_dt), exp_colnames)
        for(tid in unique(bw_dt$id)){
            test_dt = bw_dt[id == tid]
            expect_equal(nrow(test_dt), region_size / win) #expected number of regions
            bw_gr = GRanges(test_dt)
            expect_equal(length(findOverlaps(bw_gr, bw_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
        }
    }
})

# test_that("ssvFetchBigwig.single use GRanges names as id", {
#     skip_on_os("windows")
#     for(win in c(1, 3, 5, 15)){
#         qgr = test_qgr
#         names(qgr) = paste0("myNames_", seq_along(qgr))
#         bw_dt = ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = qgr, return_data.table = TRUE)
#         expect_is(bw_dt, "data.table")
#         expect_equal(colnames(bw_dt), exp_colnames)
#         expect_true(all(grepl("myNames", bw_dt$id)))
#         for(tid in unique(bw_dt$id)){
#             test_dt = bw_dt[id == tid]
#             expect_equal(nrow(test_dt), region_size / win) #expected number of regions
#             bw_gr = GRanges(test_dt)
#             expect_equal(length(findOverlaps(bw_gr, bw_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
#         }
#     }
# })

test_that("ssvFetchBigwig.single patches missing values", {
    skip_on_os("windows")
    for(win in c(1, 5, 20)){
        qgr = GRanges("chrTest", IRanges(1, 2000))
        bw_dt = ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = qgr, return_data.table = TRUE)
        expect_is(bw_dt, "data.table")
        expect_equal(colnames(bw_dt), exp_colnames)
    }
})

test_that("ssvFetchBigwig.single throws message if widths aren't divisble by win_size", {
    skip_on_os("windows")
    mix_width_gr = test_qgr
    end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)
    for(win in c(7, 31)){
        expect_message(regexp = "widths of qgr were not identical and evenly divisible by win_size", {
            ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = mix_width_gr, return_data.table = TRUE)
        })
    }
})



test_that("ssvFetchBigwig.single throws message if widths vary", {
    skip_on_os("windows")
    mix_width_gr = test_qgr
    end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)*3
    for(win in c(1, 3)){
        expect_message(regexp = "widths of qgr were not identical and evenly divisible by win_size", {
            ssvFetchBigwig.single(bw_file = test_bw, win_size = win, qgr = mix_width_gr, return_data.table = TRUE)
        })
    }
})

test_that("ssvFetchBigwig works with character vector bw_files", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = ssvFetchBigwig(file_paths = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr, return_data.table = TRUE)})
    expect_s3_class(res, "data.table")
    expect_equal(colnames(res), c(exp_colnames, "sample"))
})

test_that("ssvFetchBigwig works with list bw_files", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    bw_files = as.list(bw_files)
    hidden = capture_output({res = ssvFetchBigwig(file_paths = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr, return_data.table = TRUE)})
    expect_s3_class(res, "data.table")
    expect_equal(colnames(res), c(exp_colnames, "sample"))
})

test_that("ssvFetchBigwig can set variable name", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)
    hidden = capture_output({res = ssvFetchBigwig(file_paths = bw_files,
                                                           win_size = 3,
                                                           qgr = test_qgr,
                                                           names_variable = "group", return_data.table = TRUE)})
    expect_s3_class(res, "data.table")
    expect_equal(colnames(res), c(exp_colnames, "group"))
})

test_that("ssvFetchBigwig duplicate names throws error", {
    skip_on_os("windows")
    bw_files = rep(test_bw, 3)
    expect_error(
        ssvFetchBigwig(file_paths = bw_files,
                                win_size = 3,
                                qgr = test_qgr, return_data.table = TRUE)
    )
})

