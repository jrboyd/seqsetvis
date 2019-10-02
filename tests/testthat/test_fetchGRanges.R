testthat::context("FetchGranges")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
library(data.table)

qgr = CTCF_in_10a_overlaps_gr[1:10]
qgr = centerFixedSizeGRanges(qgr, 500)
#bed used to intersect bam


test_that("ssvFetchGRanges win_size", {
    res50 = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 50)
    expect_equal(300, length(res50))
    res10 = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 50/5)
    expect_equal(300*5, length(res10))
})


test_that("ssvFetchGRanges attribs", {
    res10 = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 10, file_attribs = data.frame("a" = 1:3, "bc" = 4:6))
    expect_true(all(res10$a %in% 1:3))
    expect_true(all(res10$bc %in% 4:6))
})

test_that("ssvFetchGRanges win_method", {
    res5sum = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 5, win_method = "summary")
    expect_equal(length(qgr)*length(CTCF_in_10a_narrowPeak_grs)*5, length(res5sum))
    expect_true(!"tile_id" %in% colnames(mcols(res5sum)))
    expect_gte(min(res5sum$x), -.5)
    expect_lte(max(res5sum$x), .5)
})


test_that("ssvFetchGRanges anchor", {
    res5left = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr, win_size = 5,
                              win_method = "summary", anchor = "left")
    expect_equal(length(qgr)*length(CTCF_in_10a_narrowPeak_grs)*5, length(res5left))
    expect_gte(min(res5left$x), 0)
    expect_lte(max(res5left$x), 1)
})

test_that("ssvFetchGRanges target_strand", {
    res50star = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr,
                              target_strand = "*", return_data.table = TRUE)
    expect_equal(max(res50star$y), 1)
    res50minus = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr,
                               target_strand = "-", return_data.table = TRUE)
    expect_equal(max(res50minus$y), 0)
    res50plus = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr,
                               target_strand = "+", return_data.table = TRUE)
    expect_equal(max(res50plus$y), 0)
    res50both = ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, qgr,
                               target_strand = "both", return_data.table = TRUE)
    expect_equal(max(res50both$y), 0)

    np_gr = CTCF_in_10a_narrowPeak_grs
    strand(np_gr$MCF10A_CTCF[1:3]) = "+"
    strand(np_gr$MCF10AT1_CTCF[4:6]) = "-"

    res50minus = ssvFetchGRanges(np_gr, qgr,
                                 target_strand = "-", return_data.table = TRUE)
    expect_equal(max(res50minus$y), 1)
    res50plus = ssvFetchGRanges(np_gr, qgr,
                                target_strand = "+", return_data.table = TRUE)
    expect_equal(max(res50plus$y), 1)
    res50both = ssvFetchGRanges(np_gr, qgr,
                                target_strand = "both", return_data.table = TRUE)
    expect_equal(max(res50both$y), 1)

    expect_equal(nrow(res50both), nrow(res50minus)*2)
})

test_that("ssvFetchGRanges qualitative sample", {
    qual_gr = CTCF_in_10a_narrowPeak_grs
    qual_gr = lapply(qual_gr, function(x){
        x$quality = "low"
        x[x$qValue > 150]$quality = "high"
        x[x$qValue > 250]$quality = "very high"
        x
    })

    dt = ssvFetchGRanges(qual_gr, qgr = qgr,
                         attrib_var = "quality",
                         fill_value = "no peak",
                         return_data.table = TRUE)
    dt[, .N, sample]
    expect_equal(nrow(dt[is.na(y), .N, sample]), 0)
    expect_equal(length(table(dt$y)), 4)
})

test_that("ssvFetchGRanges qualitative summary", {
    qual_gr = CTCF_in_10a_narrowPeak_grs
    qual_gr = lapply(qual_gr, function(x){
        x$quality = "low"
        x[x$qValue > 150]$quality = "high"
        x[x$qValue > 250]$quality = "very high"
        x
    })

    dt = ssvFetchGRanges(qual_gr, qgr = qgr,
                         win_method = "summary",
                         attrib_var = "quality",
                         fill_value = "no peak",
                         return_data.table = TRUE)
    expect_true(!"tile_id" %in% colnames(dt))
    dt[, .N, sample]
    dt$quality = factor(dt$quality, levels = c("no peak", "low", "high", "very high"))
    # ggplot(dt, aes(x = x, y = y, group = id)) + geom_path() + facet_grid(quality~sample)
    ggplot(dt, aes(x = x, y = id, fill = y, group = id)) + geom_raster() + facet_grid(quality~sample)
    expect_equal(nrow(dt[is.na(y), .N, sample]), 0)

    expect_true(all(c("quality", "y") %in% colnames(dt)))
})

test_that("ssvFetchGRanges target_strand = 'both' - coverage", {
    qual_gr = CTCF_in_10a_narrowPeak_grs
    qual_gr = lapply(qual_gr, function(x){
        x$quality = "low"
        x[x$qValue > 150]$quality = "high"
        x[x$qValue > 250]$quality = "very high"
        strand(x) = "+"
        x
    })

    dt = ssvFetchGRanges(qual_gr, qgr = qgr,
                         win_method = "summary",
                         target_strand = "both",
                         fill_value = 0,
                         #attrib_var = "quality",
                         #fill_value = "no peak",
                         return_data.table = TRUE)
    expect_equal(max(dt$y), 1)
    expect_equal(min(dt$y), 0)

    dt = ssvFetchGRanges(qual_gr, qgr = qgr,
                         win_method = "sample",
                         target_strand = "both", fill_value = 0,
                         #attrib_var = "quality",
                         #fill_value = "no peak",
                         return_data.table = TRUE)
    expect_equal(max(dt$y), 1)
    expect_equal(min(dt$y), 0)
})

test_that("ssvFetchGRanges target_strand = 'both' - quality", {
    qual_gr = CTCF_in_10a_narrowPeak_grs
    qual_gr = lapply(qual_gr, function(x){
        x$quality = "low"
        x[x$qValue > 150]$quality = "high"
        x[x$qValue > 250]$quality = "very high"
        strand(x) = "+"
        x
    })

    dt = ssvFetchGRanges(qual_gr, qgr = qgr,
                         win_method = "summary",
                         target_strand = "both",
                         attrib_var = "quality",
                         fill_value = "no peak",
                         return_data.table = TRUE)

    expect_equal(table(dt[strand == "+"]$quality) %>% length, 4)
    expect_equal(table(dt[strand == "-"]$quality) %>% length, 1)

    dt = ssvFetchGRanges(qual_gr, qgr = qgr,
                         win_method = "sample",
                         target_strand = "both",
                         attrib_var = "quality",
                         fill_value = "no peak",
                         return_data.table = TRUE)
    expect_equal(table(dt[strand == "+"]$y) %>% length, 4)
    expect_equal(table(dt[strand == "-"]$y) %>% length, 1)
})
