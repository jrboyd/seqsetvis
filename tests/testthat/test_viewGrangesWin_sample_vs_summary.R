# flipping viewGranges
library(seqsetvis)
library(testthat)

qgr = CTCF_in_10a_overlaps_gr[1:5]
strand(qgr) = c("+", "-", "-", "+", "-")
# qgr = centerFixedSizeGRanges(qgr, 500)
#bed used to intersect bam
# rtracklayer::export.bed(qgr, con = "ctcf_5.bed")
bam_file = system.file("extdata/test.bam", package = "seqsetvis")

bam_gr = fetchBam(bam_file, qgr)

#sampling
samp_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, x0 = "center")
samp_dt$group = "stranded"
samp_dt_us = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, x0 = "center_unstranded")
samp_dt_us$group = "unstranded"

ssvSignalLineplot(rbind(samp_dt, samp_dt_us), sample_ = "id", color_ = "strand", group_ = "group")

samp_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, x0 = "left")
samp_dt$group = "stranded"
samp_dt_us = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 50, x0 = "left_unstranded")
samp_dt_us$group = "unstranded"

ssvSignalLineplot(rbind(samp_dt, samp_dt_us), sample_ = "id", color_ = "strand", group_ = "group")

#summary
samp_dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, x0 = "center")
samp_dt$group = "stranded"
samp_dt_us = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, x0 = "center_unstranded")
samp_dt_us$group = "unstranded"

ssvSignalLineplot(rbind(samp_dt, samp_dt_us), sample_ = "id", color_ = "strand", group_ = "group")

samp_dt = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, x0 = "left")
samp_dt$group = "stranded"
samp_dt_us = viewGRangesWinSummary_dt(bam_gr, qgr, n_tiles = 10, x0 = "left_unstranded")
samp_dt_us$group = "unstranded"

ssvSignalLineplot(rbind(samp_dt, samp_dt_us), sample_ = "id", color_ = "strand", group_ = "group")


# fetchBam and fetchBw params
# win_method
# summary_FUN
