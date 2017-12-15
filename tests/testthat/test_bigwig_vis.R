# library(peakvisr)
library(testthat)
# library(rtracklayer)
test_bw = system.file("extdata/test_bigwigs/test_loading.bw", package = "peakvisr", mustWork = T)
pos = c(20, 180, 210, 440, 520, 521)
region_size = 30
test_qgr = GRanges("chrTest", IRanges(pos+1, pos + region_size))

test_that("fetchWindowedBigwig return expected for valid even win_size", {
  #these should all work cleanly
  for(win in c(2, 6, 10, 30)){
    bw_dt = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = test_qgr)
    expect_is(bw_dt, "data.table")
    for(tid in unique(bw_dt$id)){
      test_dt = bw_dt[id == tid]
      expect_equal(nrow(test_dt), region_size / win) #expected number of regions
      bw_gr = GRanges(test_dt)
      expect_equal(length(findOverlaps(bw_gr, bw_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
    }
  }
})

test_that("fetchWindowedBigwig return expected for valid odd win_size", {
  for(win in c(1, 3, 5, 15)){
    bw_dt = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = test_qgr)
    expect_is(bw_dt, "data.table")
    for(tid in unique(bw_dt$id)){
      test_dt = bw_dt[id == tid]
      expect_equal(nrow(test_dt), region_size / win) #expected number of regions
      bw_gr = GRanges(test_dt)
      expect_equal(length(findOverlaps(bw_gr, bw_gr)), region_size / win, info = "each returned GRange should only intersect itself.")
    }
  }
})

test_that("fetchWindowedBigwig throws error if widths aren't divisble by win_size", {
  mix_width_gr = test_qgr
  end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)
  for(win in c(7, 11, 4, 8, 20)){
    expect_error({
      fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = mix_width_gr)
    })
  }
})

test_that("fetchWindowedBigwig throws warning if widths vary", {
  mix_width_gr = test_qgr
  end(mix_width_gr) =  end(mix_width_gr) + seq_along(mix_width_gr)*3
  for(win in c(1, 3)){
    expect_warning({
      fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = mix_width_gr)
    })
  }
})


