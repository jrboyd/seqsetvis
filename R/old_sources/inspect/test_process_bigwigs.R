library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(ggplot2)
source("functions_process_bw.R")
MCF7bza_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs", pattern = "MCF7_bza_.+_FE.bw", full.names = T)
names(MCF7bza_bws) = sub("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_bza_", "", MCF7bza_bws) %>%
  sub(pattern = "_FE.bw", replacement = "")

###load test data
np = read.table("/slipstream/galaxy/uploads/working/qc_framework/output_AF_RUNX1_ChIP/AF-MCF10AT1_RUNX1_pooled/AF-MCF10AT1_RUNX1_pooled_peaks_passIDR.05.narrowPeak")
np = np[,1:4]
colnames(np) = c("seqnames", "start", "end", "id")
np_gr = GRanges(np)
mids = floor(start(np_gr) + width(np_gr)/2)
ext = 2000
start(np_gr) = mids - ext + 1
end(np_gr) = mids + ext
set.seed(0)
rand_i = lapply(1:1000, function(x){
  sample(1:length(np_gr), 1000)
})
i = 1
test_gr = np_gr[rand_i[[i]]]

bw_file = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/10A_progression/AF-MCF10AT1_RUNX1_pooled_FE.bw"
###method 1 - import ranges from bigwig
qgr = np_gr
# qdt = as.data.table(findOverlaps(qgr, qgr))
# a = qdt[queryHits != subjectHits][1,]$queryHits
# b = qdt[queryHits != subjectHits][1,]$subjectHits
# qgr = qgr[c(a,b)]
bw_gr = fetch_windowed_bw(bw_file = bw_file, qgr = qgr)
bw_dt = bw_gr2dt(bw_gr, qgr)
gg_bw_banded_quantiles(bw_dt)

w = 25
bw_gr = fetch_windowed_bw(bw_file = bw_file, win_size = w, qgr = qgr)
bw_dt = bw_gr2dt(bw_gr, qgr, win_size = w)
gg_bw_banded_quantiles(bw_dt, win_size = w)

w = 10
bw_gr = fetch_windowed_bw(bw_file = bw_file, win_size = w, qgr = qgr)
bw_dt = bw_gr2dt(bw_gr, qgr, win_size = w)
gg_bw_banded_quantiles(bw_dt, win_size = w)
