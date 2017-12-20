# library(peakvisr)
library(testthat)
# library(rtracklayer)
test_bw = system.file("extdata/test_bigwigs/test_loading.bw", package = "peakvisr", mustWork = T)
pos = c(20, 180, 210, 440, 520, 521)
region_size = 30
test_qgr = GRanges("chrTest", IRanges(pos+1, pos + region_size))
exp_colnames = c("seqnames", "start", "end", "width", "strand", "id", "FE", "x")

test_that("fetchWindowedBigwig return expected for valid even win_size", {
  #these should all work cleanly
  for(win in c(2, 6, 10, 30)){
    bw_dt = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = test_qgr)
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

test_that("fetchWindowedBigwig return expected for valid odd win_size", {
  for(win in c(1, 3, 5, 15)){
    bw_dt = fetchWindowedBigwig(bw_file = test_bw, win_size = win, qgr = test_qgr)
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

test_that("fetchWindowedBigwigList works with proper inputs", {
  bw_files = rep(test_bw, 3)
  names(bw_files) = paste0("bw_", 1:3)
  hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                         win_size = 3,
                                                         qgr = test_qgr)})
  expect_s3_class(res, "data.table")
  expect_equal(colnames(res), c(exp_colnames, "sample"))
})

test_that("fetchWindowedBigwigList can set variable name", {
  bw_files = rep(test_bw, 3)
  names(bw_files) = paste0("bw_", 1:3)
  hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                         win_size = 3,
                                                         qgr = test_qgr,
                                                         bw_variable_name = "group")})
  expect_s3_class(res, "data.table")
  expect_equal(colnames(res), c(exp_colnames, "group"))
})

test_that("fetchWindowedBigwigList duplicate names throws error", {
  bw_files = rep(test_bw, 3)
  expect_error(
    fetchWindowedBigwigList(bw_files = bw_files,
                            win_size = 3,
                            qgr = test_qgr)
  )
})

test_that("fetchWindowedBigwigList works with proper inputs", {
  bw_files = rep(test_bw, 3)
  names(bw_files) = paste0("bw_", 1:3)
  for(win in c(1, 3)){
    hidden = capture_output({res = fetchWindowedBigwigList(bw_files = bw_files,
                                                           win_size = win,
                                                           qgr = test_qgr,
                                                           bw_variable_name = "group")})
    expect_s3_class(res, "data.table")
    expect_equal(colnames(res), c(exp_colnames, "sample"))
  }
})

a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1:7 + 10))
b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1:6 + 8))

test_that("centerFixedSizeGRanges final width equals fixed_size", {
  expect_equal(width(centerFixedSizeGRanges(a, 2)), rep(2, length(a)))
  expect_equal(width(centerFixedSizeGRanges(a, 5)), rep(5, length(a)))
  expect_equal(width(centerFixedSizeGRanges(a, 80)), rep(80, length(a)))
  expect_equal(width(centerFixedSizeGRanges(b, 1)), rep(1, length(b)))
  expect_equal(width(centerFixedSizeGRanges(b, 7)), rep(7, length(b)))
  expect_equal(width(centerFixedSizeGRanges(b, 50)), rep(50, length(b)))
})
test_that("centerFixedSizeGRanges size shifts are centered", {
  a_larger = centerFixedSizeGRanges(grs = a, fixed_size = max(width(a)) + 2)
  expect_true(all(start(a_larger) < start(a)))
  expect_true(all(end(a_larger) > end(a)))
  a_smaller = centerFixedSizeGRanges(grs = a_larger, fixed_size = min(width(a)) - 2)
  expect_true(all(start(a_larger) < start(a_smaller)))
  expect_true(all(end(a_larger) > end(a_smaller)))
})
test_that("centerFixedSizeGRanges size sifts are reversible", {
  #primarily concerened about impacts of rounding
  a4 = centerFixedSizeGRanges(grs = a, fixed_size = 4)
  a5 = centerFixedSizeGRanges(grs = a, fixed_size = 5)
  a6 = centerFixedSizeGRanges(grs = a, fixed_size = 6)
  a7 = centerFixedSizeGRanges(grs = a, fixed_size = 7)
  a9 = centerFixedSizeGRanges(grs = a, fixed_size = 9)
  #derivations
  a7_from5 = centerFixedSizeGRanges(grs = a5, fixed_size = 7)
  a4_from5 = centerFixedSizeGRanges(grs = a5, fixed_size = 4)
  a6_from4 = centerFixedSizeGRanges(grs = a4, fixed_size = 6)
  a9_from4 = centerFixedSizeGRanges(grs = a4, fixed_size = 9)
  a7_from4 = centerFixedSizeGRanges(grs = a4, fixed_size = 7)
  a5_from4 = centerFixedSizeGRanges(grs = a4, fixed_size = 5)
  #reversals
  a5_from7rev = centerFixedSizeGRanges(grs = a7_from5, fixed_size = 5)
  a5_from4rev = centerFixedSizeGRanges(grs = a4_from5, fixed_size = 5)
  a4_from6rev = centerFixedSizeGRanges(grs = a6_from4, fixed_size = 4)
  a4_from9rev = centerFixedSizeGRanges(grs = a9_from4, fixed_size = 4)
  #recover even from even
  expect_true(all(a4 == a4_from6rev))
  expect_true(all(a6 == a6_from4))

  #recover odd from even
  expect_true(all(a5 == a5_from4rev))
  expect_true(all(a9 == a9_from4))#fail
  expect_true(all(a7 == a7_from4))
  expect_true(all(a5 == a5_from4))

  #recover odd from odd
  expect_true(all(a5 == a5_from7rev))
  expect_true(all(a7 == a7_from5))

  #recover even from odd
  expect_true(all(a4 == a4_from5))
  expect_true(all(a4 == a4_from9rev))
})

