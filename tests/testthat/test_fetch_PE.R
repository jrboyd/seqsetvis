testthat::context("ssvFetchBamPE")

library(seqsetvis)
library(GenomicRanges)
library(data.table)
library(testthat)
pe_file = system.file("extdata/Bcell_PE.mm10.bam", package = "seqsetvis", mustWork = TRUE)

data("Bcell_peaks")
qgr = Bcell_peaks

pe_raw = ssvFetchBamPE(
    pe_file,
    qgr,
    return_data.table = TRUE,
    win_size = 50,
    return_unprocessed = TRUE
)

exp_cn = c("which_label", "seqnames",
           "qname", "strand", "start",
           "width", "cigar", "isize", "id",
           "sample", "flag", "mapq",
           "mrnm", "mpos", "seq", "qual")

test_that("ssvFetchBamPE return_unprocessed", {
    expect_setequal(
        colnames(pe_raw),
        exp_cn
    )
})

pe_dt = ssvFetchBamPE(pe_file,
                      qgr,
                      return_data.table = TRUE,
                      win_size = 50)

test_that("ssvFetchBamPE normal", {
    expect_equal(
        colnames(pe_dt),
        c(
            "seqnames",
            "start",
            "end",
            "width",
            "strand",
            "id",
            "y",
            "x",
            "sample"
        )
    )
})

pe_raw = ssvFetchBamPE(
    pe_file,
    qgr,
    return_data.table = TRUE,
    win_size = 50,
    return_fragSizes = TRUE
)


test_that("viewGRangesWinSample_dt strand and position functions", {
    bam_score = seqsetvis:::fetchBamPE(pe_file, qgr = qgr)
    score_gr = bam_score
    window_size = 50
    anchor = "center"
    qgr_stranded = qgr
    GenomicRanges::strand(qgr_stranded) = c(rep("+", 2), rep("-", 2))

    b_dt_center = viewGRangesWinSample_dt(bam_score, qgr_stranded, 50, anchor = "center")

    b_dt_center_uns = viewGRangesWinSample_dt(bam_score, qgr_stranded, 50, anchor = "center_unstranded")

    b_dt_left = viewGRangesWinSample_dt(bam_score, qgr_stranded, 50, anchor = "left")

    b_dt_left_uns = viewGRangesWinSample_dt(bam_score, qgr_stranded, 50, anchor = "left_unstranded")

    b_dt = rbindlist(list(center = b_dt_center,
                          center_unstranded = b_dt_center_uns,
                          left = b_dt_left,
                          left_unstranded = b_dt_left_uns), use.names = TRUE, idcol = "group")
    b_dt[, facet_label := paste(id, strand)]
    cowplot::plot_grid(
        ggplot(b_dt[grepl("center", group)]) + geom_path(aes(x = x, y = y, color = group)) + facet_wrap("facet_label") + labs(title = "centered"),
        ggplot(b_dt[!grepl("center", group)]) + geom_path(aes(x = x, y = y, color = group)) + facet_wrap("facet_label") + labs(title = "left")

    )


    #verify stranded equal for + and not equal for -
    expect_true(all(b_dt_center[strand == "+"]$x == b_dt_center_uns[strand == "+"]$x))
    expect_true(all(!b_dt_center[strand == "-"]$x == b_dt_center_uns[strand == "-"]$x))
    expect_true(all(b_dt_left[strand == "+"]$x == b_dt_left_uns[strand == "+"]$x))
    expect_true(all(!b_dt_left[strand == "-"]$x == b_dt_left_uns[strand == "-"]$x))
    #verify center and left not equal
    expect_true(all(!b_dt_center$x == b_dt_left$x))
    expect_true(all(!b_dt_center_uns$x == b_dt_left_uns$x))
    #verify center and left ARE equal if minus min
    expect_true(all(b_dt_center$x - min(b_dt_center$x) == b_dt_left$x - min(b_dt_left$x)))
    expect_true(all(b_dt_center_uns$x - min(b_dt_center_uns$x) == b_dt_left_uns$x - min(b_dt_left_uns$x)))
})


# pe_raw[isize > 0, mean(isize), .(id)]
# pe_raw$id = factor(pe_raw$id, levels = qgr$name)
# ggplot(pe_raw[isize > 0, ], aes(x = isize)) + geom_histogram() + facet_wrap("id")
# plist_pe = lapply(levels(pe_raw$id), function(qid) {
#     ggplot(pe_raw[id == qid &
#                       !is.na(strand)], aes(x = pos)) +
#         geom_histogram(bins = 30) +
#         facet_wrap("strand", ncol = 1) +
#         labs(title = qid)
#
# })
# cowplot::plot_grid(plotlist = plist_pe)

#
# se_dt_base = ssvFetchBam(
#     pe_file,
#     qgr,
#     return_data.table = TRUE,
#     win_size = 50,
#     fragLens = NA
# )
# se_dt_base[, y := y / 2]
# se_dt = ssvFetchBam(
#     pe_file,
#     qgr,
#     return_data.table = TRUE,
#     target_strand = "both",
#     win_size = 50,
#     fragLens = NA,
#     flip_strand = TRUE
# )
# pe_dt = ssvFetchBamPE(pe_file, qgr, return_data.table = TRUE, win_size = 50)
# pe_dt$strand = "pe"
#
# dt = rbindlist(
#     list(
#         base = se_dt_base,
#         stranded = se_dt,
#         pe = pe_dt
#     ),
#     use.names = TRUE,
#     idcol = "group"
# )
# dt$id = factor(dt$id, levels = qgr$name)
# ggplot(dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id")#facet_grid("group~id")

# ssvFetchBam(pe_file, qgr, win_size = 50,
#             target_strand = "both", return_data.table = TRUE, flip_strand = TRUE)

# strand_dt1 = ssvFetchBam(pe_file, qgr, fragLens = NA, win_size = 50,
#                         target_strand = "both", return_data.table = TRUE, flip_strand = FALSE)
# p1 = ggplot(strand_dt1, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id")+
#     labs(title = "noflip NA") + scale_color_manual(values = c("+" = "red", "-" = "blue"))
# strand_dt2 = ssvFetchBam(pe_file, qgr, fragLens = 138, win_size = 50,
#                         target_strand = "both", return_data.table = TRUE, flip_strand = FALSE)
# p2 = ggplot(strand_dt2, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id")+
#     labs(title = "noflip 138") + scale_color_manual(values = c("+" = "red", "-" = "blue"))
#
# strand_dt3 = ssvFetchBam(pe_file, qgr, fragLens = NA, win_size = 50,
#                         target_strand = "both", return_data.table = TRUE, flip_strand = TRUE)
# p3 = ggplot(strand_dt3, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id")+
#     labs(title = "flip NA") + scale_color_manual(values = c("+" = "red", "-" = "blue"))
# strand_dt4 = ssvFetchBam(pe_file, qgr, fragLens = 138, win_size = 50,
#                         target_strand = "both", return_data.table = TRUE, flip_strand = TRUE)
#
# p4 = ggplot(strand_dt4, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id") +
#     labs(title = "flip 138") + scale_color_manual(values = c("+" = "red", "-" = "blue"))
# cowplot::plot_grid(p1, p2, p3, p4)
#
# plist_noflip = lapply(seq_along(qgr), function(i){
#     gr = qgr[i]
#     p = seqsetvis::fragLen_calcStranded(
#         pe_file,
#         gr,
#         flip_strand = FALSE,
#         test_fragLen = seq(0, 600, by = 50),
#         include_plot_in_output = TRUE
#     )[[2]]
#     p + labs(title = qgr$name[i])
# })
# plist_flip = lapply(seq_along(qgr), function(i){
#     gr = qgr[i]
#     res = seqsetvis::fragLen_calcStranded(
#         pe_file,
#         gr,
#         flip_strand = TRUE,
#         test_fragLen = seq(0, 600, by = 50),
#         include_plot_in_output = TRUE
#     )
#     p = res[[2]]
#     p + labs(title = qgr$name[i])
# })
#
# shift_dt = crossCorrByRle(pe_file, qgr, flip_strand = TRUE)
# ggplot(shift_dt, aes(x = shift, y = correlation)) + geom_path() + facet_wrap("id")

# cowplot::plot_grid(plotlist = plist_flip)
# cowplot::plot_grid(plotlist = plist_noflip)
