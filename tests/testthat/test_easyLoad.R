testthat::context("EasyLoad")
library(seqsetvis)
library(testthat)
library(GenomicRanges)


fail_tests = function(met_nam){
    met = get(met_nam)

    test_that(paste(met_nam, "error on non-character inputs"), {
        exp = "is.character\\(file_paths\\) is not TRUE"
        expect_error(met(1), exp)
        expect_error(met(1L), exp)
        expect_error(met(mean), exp)
        expect_error(met(list()), exp)
    })

    test_that(paste(met_nam, "error on non-existant file"), {
        exp = "all\\(file.exists.+is not TRUE"
        bad_file = "not a file.anywhere.ever.beddybye"
        expect_error(met(bad_file), exp)
        expect_error(met(bad_file), exp)
        expect_error(met(bad_file), exp)
        expect_error(met(bad_file), exp)
    })
}

for(met_nam in c("easyLoad_bed",
                 "easyLoad_broadPeak",
                 "easyLoad_narrowPeak")){
    fail_tests(met_nam)
}

test_dir = system.file("extdata/", package = "seqsetvis")
test_files = c("test_loading.narrowPeak", "test_loading.broadPeak", "test_loading.bed")
test_files = paste0(test_dir, "/", test_files)
stopifnot(all(file.exists(test_files)))

test_that("easyLoad_narrowPeak works, 1 file", {
    tgr = easyLoad_narrowPeak(test_files[1])
    expect_is(tgr, "list")
    expect_is(tgr[[1]], "GRanges")
    expect_equal(length(tgr), 1)
    expect_equal(names(tgr), "test_loading.narrowPeak")
    expect_equal(ncol(mcols(tgr[[1]])), 6)
})

test_that("easyLoad_narrowPeak works, 3 files and names", {
    tgr = easyLoad_narrowPeak(rep(test_files[1],3), file_names = LETTERS[1:3])
    expect_is(tgr, "list")
    expect_is(tgr[[1]], "GRanges")
    expect_equal(length(tgr), 3)
    expect_equal(names(tgr), LETTERS[1:3])
    expect_equal(ncol(mcols(tgr[[1]])), 6)
})

test_that("easyLoad_narrowPeak works, 1 file", {
    tgr = easyLoad_broadPeak(test_files[2])
    expect_is(tgr, "list")
    expect_is(tgr[[1]], "GRanges")
    expect_equal(length(tgr), 1)
    expect_equal(names(tgr), "test_loading.broadPeak")
    expect_equal(ncol(mcols(tgr[[1]])), 5)
})

test_that("easyLoad_narrowPeak works, 3 files and names", {
    tgr = easyLoad_broadPeak(rep(test_files[2],3), file_names = LETTERS[1:3])
    expect_is(tgr, "list")
    expect_is(tgr[[1]], "GRanges")
    expect_equal(length(tgr), 3)
    expect_equal(names(tgr), LETTERS[1:3])
    expect_equal(ncol(mcols(tgr[[1]])), 5)
})

test_that("easyLoad_bed works, 1 file", {
    tgr = easyLoad_bed(test_files[3])
    expect_is(tgr, "list")
    expect_is(tgr[[1]], "GRanges")
    expect_equal(length(tgr), 1)
    expect_equal(names(tgr), "test_loading.bed")
    expect_equal(ncol(mcols(tgr[[1]])), 2)
})

test_that("easyLoad_bed works, 3 files and names", {
    tgr = easyLoad_bed(rep(test_files[3],3), file_names = LETTERS[1:3])
    expect_is(tgr, "list")
    expect_is(tgr[[1]], "GRanges")
    expect_equal(length(tgr), 3)
    expect_equal(names(tgr), LETTERS[1:3])
    expect_equal(ncol(mcols(tgr[[1]])), 2)
})

test_that("easyLoad_bed empty file", {
    empty_bed_file = dir(test_dir, pattern = "empty.bed", full.names = TRUE)
    bed_empty = easyLoad_bed(empty_bed_file)[[1]]
    bed_normal = easyLoad_bed(test_files[3])[[1]]
    expect_is(bed_empty, "GRanges")
    expect_equal(length(bed_empty), 0)
    expect_equal(
        colnames(GenomicRanges::mcols(bed_empty)),
        colnames(GenomicRanges::mcols(bed_normal))
    )

})

test_that("easyLoad_bed empty file", {
    empty_np_file = dir(test_dir, pattern = "empty.narrowPeak", full.names = TRUE)
    np_empty = easyLoad_narrowPeak(empty_np_file)[[1]]
    np_normal = easyLoad_narrowPeak(test_files[1])[[1]]
    expect_is(np_empty, "GRanges")
    expect_equal(length(np_empty), 0)
    expect_equal(
        colnames(GenomicRanges::mcols(np_empty)),
        colnames(GenomicRanges::mcols(np_normal))
    )
})
