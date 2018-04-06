library(seqsetvis)

test_that("ssvFactorizeMembTable expected output format for GRanges", {
    gr = CTCF_in_10a_overlaps_gr
    names(gr) = paste("peak", seq_along(gr))
    fac = ssvFactorizeMembTable(gr)
    expect_true(is.data.frame(fac))
    expect_equal(colnames(fac), c("id", "group"))
    expect_equal(nrow(fac), length(CTCF_in_10a_overlaps_gr))
    expect_true(is.factor(fac$group))
    expect_true(is.factor(fac$id))
})

test_that("ssvFactorizeMembTable expected output values for GRanges", {
    gr = CTCF_in_10a_overlaps_gr
    fac = ssvFactorizeMembTable(gr)
    expect_equal(as.character(fac$group[2]), "MCF10A_CTCF")
    expect_equal(as.character(fac$group[5]), "MCF10A_CTCF & MCF10AT1_CTCF")
    expect_equal(as.character(fac$group[10]),  "MCF10A_CTCF & MCF10AT1_CTCF & MCF10CA1_CTCF")
})

test_that("ssvFactorizeMembTable expected output values for list", {
    setL = list(1:3, 2:3, c(3:5))
    fac = ssvFactorizeMembTable(setL)
    expect_equal(nrow(fac), 5)
    expect_equal(as.character(fac$group[1]), "set_A")
    expect_equal(as.character(fac$group[2]),  "set_A & set_B")
    expect_equal(as.character(fac$group[3]),  "set_A & set_B & set_C")
    expect_equal(as.character(fac$group[4]),  "set_C")
    expect_equal(as.character(fac$group[5]),  "set_C")
})

