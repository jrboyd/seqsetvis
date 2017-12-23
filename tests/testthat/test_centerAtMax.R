library(seqsetvis)
library(testthat)

n = 8
xs = (1:n-5)*5
fun = function(x)(x^2)
xs = c(xs, xs, xs)
ys = -c(fun(xs + floor((seq_along(xs)-1) / n)*5)) + floor((seq_along(xs)-1) / n)*50
#setup data.table of parabolic curves, each transposed on the x-axis
test_dt = data.table::data.table(xvals = xs, yvals = ys, grp = rep(letters[1:3], each = n))
ggplot(test_dt, aes(x = xvals, y = yvals, col = grp)) + geom_line() + geom_point()

#should set by
cm_dt_noBy = centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by_ = "", check_by_dupes = F)
ggplot(cm_dt_noBy, aes(x = xvals, y = yvals, col = grp)) +
  geom_line() +
  geom_point() +
  geom_point(data = test_dt, aes(col = grp), alpha = .3) +
  geom_line(data = test_dt, aes(col = grp), alpha = .3) +
  labs(title = "With by_ unset, all profiles are shifted to the global maximum")
#proper use of by
cm_dt_wBy =centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by = "grp")
ggplot(cm_dt_wBy, aes(x = xvals, y = yvals, col = grp)) +
  geom_line() +
  geom_point() +
  geom_point(data = test_dt, aes(col = grp), alpha = .3) +
  geom_line(data = test_dt, aes(col = grp), alpha = .3) +
  labs(title = "With by_ set properly, all profiles are shifted to their individual maximum")
#effect of trimming
cm_dt_wBy_noTrim =centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by = "grp", trim_to_valid = F)
ggplot(cm_dt_wBy, aes(x = xvals, y = yvals, col = grp)) +
  geom_line() +
  geom_point() +
  geom_point(data = cm_dt_wBy_noTrim, aes(col = grp), alpha = .3) +
  geom_line(data = cm_dt_wBy_noTrim, aes(col = grp), alpha = .3) +
  labs(title = "When trim_to_valid is TRUE, non-universal values of x are removed.")

test_that("centerAtMax warnings for by_ specification", {
  expect_message(centerAtMax(test_dt, x_ = "xvals", y_ = "yvals"), regexp = "centerAtMax")
  expect_failure(expect_message(centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by_ = "grp")))
})

test_that("centerAtMax errors for unmatched variable name", {
  expect_error(
    centerAtMax(test_dt),
    regexp = "centerAtMax")
  expect_error(
    centerAtMax(test_dt, x_ = "xvals"),
    regexp = "centerAtMax")
  expect_error(
    centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by_ = "blerp"),
    regexp = "centerAtMax")
})

test_that("centerAtMax trimming reduces ranges", {
  cm_noTrim = centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by_ = "grp", trim_to_valid = F)
  expect_equal(nrow(cm_noTrim), nrow(test_dt))
  cm_wTrim = centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by_ = "grp", trim_to_valid = T)
  expect_lt(nrow(cm_wTrim), nrow(test_dt))
})



