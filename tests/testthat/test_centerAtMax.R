library(seqsetvis)
library(testthat)
library(data.table)
n = 8
xs = (1:n-5)*5
fun = function(x)(x^2)
xs = c(xs, xs, xs)
ys = -c(fun(xs + floor((seq_along(xs)-1) / n)*5)) + floor((seq_along(xs)-1) / n)*50
#setup data.table of parabolic curves, each transposed on the x-axis
test_dt = data.table::data.table(xvals = xs, yvals = ys, locus = rep(letters[1:3], each = n))
ggplot(test_dt, aes(x = xvals, y = yvals, col = locus)) + geom_line() + geom_point()

test_dt2 = rbind(test_dt, test_dt)
test_dt2$sample = rep(1:2, each = nrow(test_dt))
test_dt2[sample == 2, xvals := xvals + 10]
test_dt2[sample == 2, yvals := yvals + 300]
ggplot(test_dt2, aes(x = xvals, y = yvals, col = locus)) + geom_line() + geom_point() + facet_grid(sample ~ .)
#should set by
cm_dt_noBy = centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by_ = "", check_by_dupes = F)
ggplot(cm_dt_noBy, aes(x = xvals, y = yvals, col = locus)) +
    geom_line() +
    geom_point() +
    geom_point(data = test_dt, aes(col = locus), alpha = .3) +
    geom_line(data = test_dt, aes(col = locus), alpha = .3) +
    labs(title = "With by_ unset, all profiles are shifted to the global maximum")
#proper use of by
cm_dt_wBy = centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by_ = "locus")
ggplot(cm_dt_wBy, aes(x = xvals, y = yvals, col = locus)) +
    geom_line() +
    geom_point() +
    geom_point(data = test_dt, aes(col = locus), alpha = .3) +
    geom_line(data = test_dt, aes(col = locus), alpha = .3) +
    labs(title = "With by_ set properly, all profiles are shifted to their individual maximum")
#effect of trimming
cm_dt_wBy_noTrim =centerAtMax(test_dt, x_ = "xvals", y_ = "yvals", by = "locus", trim_to_valid = F)
ggplot(cm_dt_wBy, aes(x = xvals, y = yvals, col = locus)) +
    geom_line() +
    geom_point() +
    geom_point(data = cm_dt_wBy_noTrim, aes(col = locus), alpha = .3) +
    geom_line(data = cm_dt_wBy_noTrim, aes(col = locus), alpha = .3) +
    labs(title = "When trim_to_valid is TRUE, non-universal values of x are removed.")



doTests_centerAtMax = function(test_object, test_object2){
    test_that("centerAtMax warnings for by_ specification", {
        expect_message(centerAtMax(test_object, x_ = "xvals", y_ = "yvals"), regexp = "centerAtMax")
        expect_failure(expect_message(centerAtMax(test_object, x_ = "xvals", y_ = "yvals", by_ = "locus")))
    })

    test_that("centerAtMax errors for unmatched variable name", {
        expect_error(
            centerAtMax(test_object),
            regexp = "centerAtMax")
        expect_error(
            centerAtMax(test_object, x_ = "xvals"),
            regexp = "centerAtMax")
        expect_error(
            centerAtMax(test_object, x_ = "xvals", y_ = "yvals", by_ = "blerp"),
            regexp = "centerAtMax")
    })

    test_that("centerAtMax trimming reduces ranges", {
        cm_noTrim = centerAtMax(test_object, x_ = "xvals", y_ = "yvals", by_ = "locus", trim_to_valid = F)
        expect_equal(nrow(cm_noTrim), nrow(test_object))
        cm_wTrim = centerAtMax(test_object, x_ = "xvals", y_ = "yvals", by_ = "locus", trim_to_valid = T)
        nr_test = nrow(as.data.frame(test_object))
        nr_wTrim = nrow(as.data.frame(cm_wTrim))
        nr_noTrim = nrow(as.data.frame(test_object))
        expect_lt(nr_wTrim, nr_test)
    })

    test_that("centerAtMax input must be data.table", {
        expect_error(centerAtMax(dt = data.frame(1:3), x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2")), regexp = "must be of type data.table")
        expect_error(centerAtMax(dt = matrix(1:3), x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2")), regexp = "must be of type data.table")
        expect_error(centerAtMax(dt = (1:3), x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2")), regexp = "must be of type data.table")
        expect_error(centerAtMax(dt = "(1:3)", x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2")), regexp = "must be of type data.table")
    })

    test_that("centerAtMax viewSize", {
        dt1 = centerAtMax(dt = test_object2, x_ = "xvals", y_ = "yvals", by_ = c("locus"), view_size = 5, trim_to_valid = F, check_by_dupes = F)
        dt2 = centerAtMax(dt = test_object2, x_ = "xvals", y_ = "yvals", by_ = c("locus"), view_size = 10, trim_to_valid = F, check_by_dupes = F)
        nr1 = nrow(as.data.frame(dt1))
        nr2 = nrow(as.data.frame(dt2))
        expect_equal(nr1, nr2)
        expect_true(!all(dt2$xvals == dt1$xvals))
    })

    test_that("centerAtMax viewSize", {
        dt1 = centerAtMax(dt = test_object2, x_ = "xvals", y_ = "yvals", by_ = c("locus"), check_by_dupes = F, replace_x = F)
        dt2 = centerAtMax(dt = test_object2, x_ = "xvals", y_ = "yvals", by_ = c("locus"), check_by_dupes = F, replace_x = T)

        dt1 = as.data.frame(dt1)
        dt2 = as.data.frame(dt2)

        expect_equal(colnames(dt1)[ncol(dt1)], "xvals_summitPosition")
        expect_equal(ncol(dt1), ncol(dt2)+1)
    })
}

doTests_centerAtMax(test_dt, test_dt2)
test_gr = GRanges(seqnames = "chrTest",
                  IRanges(test_dt$xvals+100, test_dt$xvals+104),
                  xvals = test_dt$xvals,
                  yvals = test_dt$yvals,
                  locus = test_dt$locus)
test_gr2 = GRanges(seqnames = "chrTest",
                  IRanges(test_dt2$xvals+100, test_dt2$xvals+104),
                  xvals = test_dt2$xvals,
                  yvals = test_dt2$yvals,
                  locus = test_dt2$locus,
                  sample = test_dt2$sample)
doTests_centerAtMax(test_object = test_gr, test_object2 = test_gr2)
doTests_centerAtMax(as.data.table(test_gr), as.data.table(test_gr2))

test_that("centerAtMax multiple by_ arguments", {
    pdt = copy(test_dt2)
    pdt$centered = "centered: no"
    cm = centerAtMax(pdt, x_ = "xvals", y_ = "yvals", by_ = c("sample", "locus"), trim_to_valid = F)
    cm$centered = "centered: yes"
    pdt =  rbind(pdt, cm)
    p = ggplot(pdt, aes(x = xvals, y = yvals, col = locus)) +
        geom_line() +
        geom_point() +
        labs(title = "With compound by_, all profiles are shifted to their individual maximum") +
        facet_grid(sample ~ centered)
    maxpos = cm[, xvals[which.max(yvals)], by = c("sample", "locus")]$V1
    expect_true(all(maxpos == 0))
})

test_that("centerAtMax sample by_ argument", {
    pdt = copy(test_dt2)
    pdt$centered = "centered: no"
    cm = centerAtMax(pdt, x_ = "xvals", y_ = "yvals", by_ = c("sample"), trim_to_valid = F, check_by_dupes = F)
    cm$centered = "centered: yes"
    pdt =  rbind(pdt, cm)
    p = ggplot(pdt, aes(x = xvals, y = yvals, col = locus)) +
        geom_line() +
        geom_point() +
        labs(title = "With sample by_, profiles are shifted equally per locus",
             subtitle = "akin to per locus across samples, a normal use") +
        facet_grid(sample ~ centered)
    maxpos = cm[, xvals[which.max(yvals)], by = "sample"]$V1
    expect_true(all(maxpos == 0))
})

test_that("centerAtMax locus by_ argument", {
    pdt = copy(test_dt2)
    pdt$centered = "centered: no"
    cm = centerAtMax(pdt, x_ = "xvals", y_ = "yvals", by_ = c("locus"), trim_to_valid = F, check_by_dupes = F)
    cm$centered = "centered: yes"
    pdt =  rbind(pdt, cm)
    p = ggplot(pdt, aes(x = xvals, y = yvals, col = locus)) +
        geom_line() +
        geom_point() +
        labs(title = "With locus by_, all profiles are shifted equally per sample",
             subtitle = "akin to per sample across loci, a weird use") +
        facet_grid(sample ~ centered)
    maxpos = cm[, xvals[which.max(yvals)], by = "locus"]$V1
    expect_true(all(maxpos == 0))
})
