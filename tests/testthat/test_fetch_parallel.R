testthat::context("parallel")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
library(data.table)

qgr = CTCF_in_10a_overlaps_gr[1:5]
qgr = centerFixedSizeGRanges(qgr, 500)
names(qgr) = NULL
#bed used to intersect bam
# rtracklayer::export.bed(qgr, con = "ctcf_5.bed")
bam_file = system.file("extdata/test.bam", package = "seqsetvis", mustWork = TRUE)
test_bw = system.file("extdata/test_loading.bw", package = "seqsetvis", mustWork = TRUE)

test_that("ssvFetchBam summary method correct bins", {
    skip_on_os("windows")
    test_qgr2 = qgr

    gr_sample_split1 = ssvFetchBam(bam_file, win_size = 5,
                                   fragLens = NA,
                                   qgr = test_qgr2, win_method = "summary",
                                   n_region_splits = 1)
    gr_sample_split2 = ssvFetchBam(bam_file, win_size = 5,
                                   fragLens = NA,
                                   qgr = test_qgr2, win_method = "summary",
                                   n_region_splits = 2)
    gr_sample_split3 = ssvFetchBam(bam_file, win_size = 5,
                                   fragLens = NA,
                                   qgr = test_qgr2, win_method = "summary",
                                   n_region_splits = 5)

    expect_equal(gr_sample_split1$id, gr_sample_split2$id)
    expect_equal(gr_sample_split1$x, gr_sample_split2$x)
    expect_equal(gr_sample_split1$y, gr_sample_split2$y)
    expect_equal(gr_sample_split1$sample, gr_sample_split2$sample)

    expect_equal(gr_sample_split1$id, gr_sample_split3$id)
    expect_equal(gr_sample_split1$x, gr_sample_split3$x)
    expect_equal(gr_sample_split1$y, gr_sample_split3$y)
    expect_equal(gr_sample_split1$sample, gr_sample_split3$sample)
})

test_that("ssvFetchBigwig works with character vector bw_files", {
    skip_on_os("windows")
    pos = c(20, 180, 210, 440, 520, 521)
    region_size = 30
    test_qgr2 = GRanges("chrTest", IRanges(pos+1, pos + region_size))
    bw_files = rep(test_bw, 3)
    names(bw_files) = paste0("bw_", 1:3)

    gr_sample_split1 = ssvFetchBigwig(bw_files, win_size = 5,
                                      qgr = test_qgr2, win_method = "summary",
                                      n_region_splits = 1, return_data.table = TRUE)
    gr_sample_split2 = ssvFetchBigwig(bw_files, win_size = 5,
                                      qgr = test_qgr2, win_method = "summary",
                                      n_region_splits = 2, return_data.table = TRUE)
    gr_sample_split3 = ssvFetchBigwig(bw_files, win_size = 5,
                                      qgr = test_qgr2, win_method = "summary",
                                      n_region_splits = 5, return_data.table = TRUE)

    expect_equal(gr_sample_split1$id, gr_sample_split2$id)
    expect_equal(gr_sample_split1$x, gr_sample_split2$x)
    expect_equal(gr_sample_split1$y, gr_sample_split2$y)
    expect_equal(gr_sample_split1$sample, gr_sample_split2$sample)

    expect_equal(gr_sample_split1$id, gr_sample_split3$id)
    expect_equal(gr_sample_split1$x, gr_sample_split3$x)
    expect_equal(gr_sample_split1$y, gr_sample_split3$y)
    expect_equal(gr_sample_split1$sample, gr_sample_split3$sample)
})
