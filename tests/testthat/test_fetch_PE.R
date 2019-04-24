testthat::context("fetchPE")

library(seqsetvis)
library(GenomicRanges)
library(data.table)
library(testthat)
pe_file = system.file("extdata/Bcell_PE.mm10.bam", package = "seqsetvis")

data("Bcell_peaks")
qgr = Bcell_peaks

pe_raw = ssvFetchBamPE(
    pe_file,
    qgr,
    return_data.table = TRUE,
    win_size = 50,
    return_unprocessed = TRUE
)

test_that("ssvFetchBamPE return_unprocessed", {
    expect_equal(
        colnames(pe_raw),
        c(
            "which_label",
            "qname",
            "rname",
            'strand',
            "pos",
            "qwidth",
            "cigar",
            "isize",
            "id",
            "sample"
        )
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
