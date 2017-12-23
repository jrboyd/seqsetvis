divideBy <- function(dividend, divisor) {
  if (divisor == 0)
    return(NA)
  dividend / divisor
}

test_that("divideBy works properly", {
  expect_equal(divideBy(4, 2), 2)
  expect_true(is.na(divideBy(4, 0)))
  expect_equal(divideBy(4, 1.2345), 3.24, tolerance = 1.0e-4)
})
