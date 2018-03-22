library(testthat)
library(seqsetvis)
library(GenomicRanges)
library(data.table)
xs = (0:5-1.5)*5
fun = function(x)x^2
test_dt2 = data.table::data.table(xvals = c(xs, xs-5, xs - 10), yvals = c(fun(xs) - 100, fun(xs-5)+60 , fun(xs - 10)+150 ), grp = rep(letters[1:3], each = length(xs)))
test_gr2 = GRanges(test_dt2[, .(seqnames = "chrTest", start = ceiling(xvals - min(xvals) - 2.5), end = ceiling(xvals - min(xvals) + 2), xvals, yvals), by = .(grp)])
test_dt3 = as.data.table(test_gr2)
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

doTest_applySpline = function(test_object){

    test_that("applySpline returned class same as input", {
        original_class = class(test_object)
        sp_dt = applySpline(test_object, x_ = "xvals", y_ = "yvals", by_ = "grp", n = 4)
        expect_equal(class(sp_dt), original_class)
    })

    test_that("applySpline returned colnames equal", {
        original_colnames = colnames(as.data.frame(test_object))
        sp_dt = applySpline(test_object, x_ = "xvals", y_ = "yvals", by_ = "grp", n = 4)
        res_colnames = colnames(as.data.frame(sp_dt))
        expect_equal(res_colnames, original_colnames)
    })

    test_that("applySpline warnings for by_ specification", {
        expect_warning(applySpline(test_object, x_ = "xvals", y_ = "yvals", n = 4), regexp = "applySpline")
        expect_failure(expect_warning(applySpline(test_object, x_ = "xvals", y_ = "yvals", n = 4, by_ = "grp")))
    })

    test_that("applySpline errors for unmatched variable name", {
        expect_error(
            applySpline(test_object),
            regexp = "applySpline")
        expect_error(
            applySpline(test_object, x_ = "xvals"),
            regexp = "applySpline")
        expect_error(
            applySpline(test_object, x_ = "xvals", y_ = "yvals", by_ = "blerp"),
            regexp = "applySpline")
    })

    test_that("applySpline n*length in result", {
        sp_dt4 = applySpline(test_object, x_ = "xvals", y_ = "yvals", by_ = "grp", n = 4)
        expect_equal(nrow(as.data.frame(test_object))*4, nrow(as.data.frame(sp_dt4)))
        sp_dt10 = applySpline(test_object, x_ = "xvals", y_ = "yvals", by_ = "grp", n = 10)
        expect_equal(nrow(as.data.frame(test_object))*10, nrow(as.data.frame(sp_dt10)))
    })

    test_that("applySpline no warning from multiple by_", {
        test_object$grp2  = test_object$grp
        expect_failure(expect_warning(applySpline(dt = test_object, x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2"), n = 4)))
    })

    test_that("applySpline error if not data.table or GRanges", {
        expect_error(applySpline(dt = data.frame(1:3), x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2"), n = 4), regexp = "data.table::is.data.table\\(dt\\) is not TRUE")
        expect_error(applySpline(dt = matrix(1:3), x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2"), n = 4), regexp = "data.table::is.data.table\\(dt\\) is not TRUE")
        expect_error(applySpline(dt = (1:3), x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2"), n = 4), regexp = "data.table::is.data.table\\(dt\\) is not TRUE")
        expect_error(applySpline(dt = "(1:3)", x_ = "xvals", y_ = "yvals", by_ = c("grp", "grp2"), n = 4), regexp = "data.table::is.data.table\\(dt\\) is not TRUE")
    })

    test_that("applySpline n < 1 throws caught error", {
        expect_error(applySpline(test_object, x_ = "xvals", y_ = "yvals", by_ = "grp", n = .9), regexp = "n >")
    })

    test_that("applySpline order doesn't matter", {
        r1 = sample(seq_len(nrow(as.data.frame(test_object))))
        r2 = sample(seq_len(nrow(as.data.frame(test_object))))
        t1 = test_object[r1]
        t2 = test_object[r2]
        # if(class(test_object)[1] == "GRanges"){
        #
        # }else{
        #
        # }
        sp_o = applySpline(test_object, x_ = "xvals", y_ = "yvals", by_ = "grp", n = 3)
        sp_1 = applySpline(t1, x_ = "xvals", y_ = "yvals", by_ = "grp", n = 3)
        sp_2 = applySpline(t2, x_ = "xvals", y_ = "yvals", by_ = "grp", n = 3)

        sp_o = sp_o[order(sp_o$xvals)]
        sp_o = sp_o[order(sp_o$grp)]

        sp_1 = sp_1[order(sp_1$xvals)]
        sp_1 = sp_1[order(sp_1$grp)]

        sp_2 = sp_2[order(sp_2$xvals)]
        sp_2 = sp_2[order(sp_2$grp)]

        expect_equal(sp_o$yvals, sp_1$yvals)
        expect_equal(as.data.frame(sp_o)$start, as.data.frame(sp_1)$start)

        expect_equal(sp_o$yvals, sp_2$yvals)
        expect_equal(as.data.frame(sp_o)$start, as.data.frame(sp_2)$start)
    })

}

#dt = copy(test_object); x_ = "xvals"; y_ = "yvals"; by_ = "grp"; n = 4; splineFun = stats::spline

doTest_applySpline(test_dt2)
doTest_applySpline(test_object = test_gr2)
doTest_applySpline(test_dt3)
