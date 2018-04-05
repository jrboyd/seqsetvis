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
        expect_error(met("not a file.anywhere.ever.beddybye"), exp)
        expect_error(met(1L), exp)
        expect_error(met(mean), exp)
        expect_error(met(list()), exp)
    })
}

for(met_nam in c("easyLoad_bed",
             "easyLoad_broadPeak",
             "easyLoad_narrowPeak")){
    fail_tests(met_nam)
}
