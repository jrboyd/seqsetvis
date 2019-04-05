testthat::context("CenterAtMax")
library(seqsetvis)
library(testthat)
library(GenomicRanges)

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
