testthat::context("normalization")

library(seqsetvis)
library(testthat)

test_that("append_ynorm defaults", {
    #these should all work cleanly
    norm_dt = append_ynorm(CTCF_in_10a_profiles_dt)
    expect_equal(sort(colnames(norm_dt)), sort(c(colnames(CTCF_in_10a_profiles_dt), "y_cap_value", "y_norm")))
})

test_that("append_ynorm manual names", {
    #these should all work cleanly
    dt = CTCF_in_10a_profiles_dt
    dt[, ytest := y]
    dt$y = NULL
    norm_dt = append_ynorm(dt, value_ = "ytest", cap_value_ = "captest", norm_value_ = "normtest")
    expect_equal(sort(colnames(norm_dt)), sort(c(colnames(dt), "captest", "normtest")))
})

test_that("append_ynorm manual cap_dt", {
    norm_dt1 = append_ynorm(CTCF_in_10a_profiles_dt)

    cap_dt = calc_norm_factors(CTCF_in_10a_profiles_dt)
    cap_dt[, y_cap_value := y_cap_value * 10]

    norm_dt2 = append_ynorm(CTCF_in_10a_profiles_dt, cap_dt = cap_dt)

    expect_true(all(norm_dt2$y_norm < norm_dt1$y_norm | norm_dt2$y_norm == 0))
})


test_that("calc_norm_factors defaults", {
    cap_dt = calc_norm_factors(CTCF_in_10a_profiles_dt)
    expect_equal(sort(colnames(cap_dt)), sort(c("sample", "y_cap_value")))
})

test_that("calc_norm_factors manual names", {
    cap_dt = calc_norm_factors(CTCF_in_10a_profiles_dt, cap_value_ = "captest", by2 = "seqnames")
    expect_equal(sort(colnames(cap_dt)), sort(c("seqnames", "captest")))
})
