testthat::context("safeBrew")

library(seqsetvis)
library(testthat)

test_that("safeBrew default", {
    cols = safeBrew(21)
    expect_equal(length(cols), 21)
    expect_null(names(cols))
    expect_equal(cols[1], "#1B9E77")
    expect_equal(cols[1], cols[9])
    expect_equal(cols[1], cols[17])
})

test_that("safeBrew character n", {
    cols = safeBrew(LETTERS)
    expect_equal(length(cols), length(LETTERS))
    expect_equal(names(cols), LETTERS)
    expect_equal(unname(cols[1]), "#1B9E77")
    expect_equal(unname(cols[1]), unname(cols[9]))
    expect_equal(unname(cols[1]), unname(cols[17]))
})

test_that("safeBrew rev character n", {
    cols = safeBrew(rev(LETTERS))
    expect_equal(length(cols), length(LETTERS))
    expect_equal(names(cols), rev(LETTERS))
    expect_equal(unname(cols[1]), "#1B9E77")
    expect_equal(unname(cols[1]), unname(cols[9]))
    expect_equal(unname(cols[1]), unname(cols[17]))
})

test_that("safeBrew rev factor n", {
    fLET = factor(LETTERS, levels = rev(LETTERS))
    cols = safeBrew(fLET)
    expect_equal(length(cols), length(LETTERS))
    expect_equal(names(cols), rev(LETTERS))
    expect_equal(unname(cols[1]), "#1B9E77")
    expect_equal(unname(cols[1]), unname(cols[9]))
    expect_equal(unname(cols[1]), unname(cols[17]))
})

test_that("safeBrew rev factor n", {
    c.D = safeBrew(3, pal = "Dark2")
    c.d = safeBrew(3, pal = "dark2")
    c.s = safeBrew(3, pal = "set1")
    c.g1 = safeBrew(3, pal = "gg")
    c.g2 = safeBrew(3, pal = "ggplot")
    c.g3 = safeBrew(3, pal = "ggplot2")
    expect_equal(c.D[1], "#1B9E77")
    expect_equal(c.D, c.d)
    expect_equal(c.s[1], "#E41A1C")
    expect_equal(c.g1[1], "#F8766D")
    expect_equal(c.g1, c.g2)
    expect_equal(c.g1, c.g3)
})
