library(GenomicRanges)
library(rtracklayer)
library(seqsetvis)
qgr = GRanges("chr19", IRanges(5*10^6, 10*10^6))
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE77772&format=file&file=GSE77772%5FMCF10A%5FH3K27ME3%5FlogFE%2Ebw"
url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69377/suppl/GSE69377_MCF10A_H3K4AC_logFE.bw"
url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69377/suppl/GSE69377_MCF7_H3K4ME3_logFE.bw"

bw = import.bw(url, which = qgr)
bw$score = round(bw$score, digits = 4)
bw$score = as.integer(1000* bw$score)
print(object.size(bw), units = "MB")

peak_files = dir("inst/extdata/peak_calls/full", pattern = "narrowPeak", full.names = T)
names(peak_files) = peak_files %>% basename %>%
    sub("_pooled_peak.+roadPeak", "", .) %>%
    sub("_pooled_peak.+arrowPeak", "", .)

peak_df = lapply(peak_files, read.table)
peak_gr = lapply(peak_df, function(x){
    colnames(x) = c("seqnames", "start", "end", "peak_id", "score", "dot", "FE", "pval", "qval", "rel_summit")
    GRanges(x)
})

memb_gr = overlapIntervalSets(peak_gr)
memb_gr = subset(memb_gr, seqnames == "chr21" & start > 5*10^6 & end < 10*10^6)
memb_gr = centerFixedSizeGRanges(memb_gr, fixed_size = 2000)
bw1 = import.bw(url, which = subset(memb_gr, MCF7_H3K4ME3))
print(object.size(bw1), units = "MB")
bw2 = fetchWindowedBigwig(url, subset(memb_gr, MCF7_H3K4ME3), win_size = 100)
print(object.size(bw2), units = "MB")
centerAtMax(bw2, y_ = "FE")
bw2$FE = 10^bw2$FE

bw2$sample = "s1"
bw2 = centerAtMax(bw2, view_size = 500, y_ = "FE", by_ = "id")
bw2$x = bw2$x - 50
regionSetPlotHeatmap(bw2, nclust = 3, max_cols = 12)
