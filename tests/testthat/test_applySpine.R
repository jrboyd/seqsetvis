library(testthat)
library(seqsetvis)
xs = (0:5-1.5)*5
fun = function(x)x^2
test_dt2 = data.table::data.table(xvals = c(xs, xs-5, xs - 10), yvals = c(fun(xs) - 100, fun(xs-5)+60 , fun(xs - 10)+150 ), grp = rep(letters[1:3], each = length(xs)))
#
# ggplot(test_dt2, aes(x = xvals, y = yvals)) + geom_point() + stat_summary(fun.y = mean, geom="line")
# ggplot(test_dt2, aes(x = xvals, y = yvals, col = grp)) + geom_line() + geom_point()
#
# #should set by
# sp_dt_noBy =applySpline(test_dt2, x_ = "xvals", y_ = "yvals", n = 4)
# ggplot(sp_dt_noBy, aes(x = xvals, y = yvals)) + geom_point() + stat_summary(fun.y = mean, geom="line") +
#   geom_point(data = test_dt2, aes(col = grp)) + stat_summary(data = test_dt2, aes(col = grp), fun.y = mean, geom="line")
# #proper use of by
# sp_dt_wBy =applySpline(test_dt2, x_ = "xvals", y_ = "yvals", n = 4, by = "grp")
# ggplot(sp_dt_wBy, aes(x = xvals, y = yvals, col = grp)) + geom_point() + stat_summary(fun.y = mean, geom="line") +
#   geom_point(data = test_dt2, aes(col = grp), alpha = .4) + stat_summary(data = test_dt2, aes(col = grp), fun.y = mean, geom="line", alpha = .4)

test_that("applySpline warnings for by_ specification", {
  expect_warning(applySpline(test_dt2, x_ = "xvals", y_ = "yvals", n = 4), regexp = "applySpline")
  expect_failure(expect_warning(applySpline(test_dt2, x_ = "xvals", y_ = "yvals", n = 4, by_ = "grp")))
})

test_that("applySpline errors for unmatched variable name", {
  expect_error(
    applySpline(test_dt2),
    regexp = "applySpline")
  expect_error(
    applySpline(test_dt2, x_ = "xvals"),
    regexp = "applySpline")
  expect_error(
    applySpline(test_dt2, x_ = "xvals", y_ = "yvals", by_ = "blerp"),
    regexp = "applySpline")
})

test_that("applySpline n*length in result", {
  sp_dt4 = applySpline(test_dt2, x_ = "xvals", y_ = "yvals", by_ = "grp", n = 4)
  expect_equal(nrow(test_dt2)*4, nrow(sp_dt4))
  sp_dt10 = applySpline(test_dt2, x_ = "xvals", y_ = "yvals", by_ = "grp", n = 10)
  expect_equal(nrow(test_dt2)*10, nrow(sp_dt10))
})

test_that("applySpline no warning from multiple by_", {
    test_dt2[, grp2 := grp]
    expect_failure(expect_warning(applySpline(dt = test_dt2, x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2"), n = 4)))
})



