###load peak sets
library(magrittr)
library(GenomicRanges)
# library(peakvisr)

peak_files = dir("~/ShinyData/peaks/", pattern = "IDR", full.names = T)
names(peak_files) = peak_files %>% basename %>%
  sub("_pooled_peaks_passIDR.05.narrowPeak", "", .) %>%
  sub("AF-", "", .)

peak_df = lapply(peak_files, read.table)
peak_gr = lapply(peak_df, function(x){
  colnames(x)[1:3] = c("seqnames", "start", "end")
  GRanges(x)
})

peak_tab = overlapIntervalSets(peak_gr)
peak_memb = elementMetadata(peak_tab)

###venn diagrams
ggVenn(peak_memb[,1:3])
ggVenn(peak_memb[,4:6])
ggVenn(peak_memb[,c(1,4)])
ggVenn(peak_memb[,c(1,4)+1])
ggVenn(peak_memb[,c(1,4)+2])

###barplots
ggBars(peak_memb)

###pie
ggPie(peak_memb)

###euler
ggEuler(peak_memb[,1:6], line_width = 2)
