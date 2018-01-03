library(seqsetvis)
library(rtracklayer)
library(magrittr)

### set file paths
np_files = c("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10A_CTCF_pooled_peaks_passIDR.05.narrowPeak.gz",
             "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10AT1_CTCF_pooled_peaks_passIDR.05.narrowPeak.gz",
             "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10CA1_CTCF_pooled_peaks_passIDR.05.narrowPeak.gz")
names(np_files) = np_files %>%
    basename() %>%
    sub("GSE98551_", "", x = .) %>%
    sub("_pooled_peaks_passIDR.05.narrowPeak.gz", "", x = .)


# bw_files = c("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10A_CTCF_pooled_FE.bw",
#              "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10AT1_CTCF_pooled_FE.bw",
#              "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10CA1_CTCF_pooled_FE.bw")

# bw_files = https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE98551&format=file&file=GSE98551%5FMCF10AT1%5FCTCF%5Fpooled%5FFE%2Ebw
names(bw_files) = bw_files %>%
    basename() %>%
    sub("GSE98551_", "", x = .) %>%
    sub("_pooled_FE.bw", "", x = .)

### load peak calls

#from: https://charlesjb.github.io/How_to_import_narrowPeak/
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
np_grs <- lapply(np_files, function(f){
    import(f, format = "BED",
           extraCols = extraCols_narrowPeak)
})

### perform 3-way intersection
np_overlap = overlapIntervalSets(np_grs)

### plots for 3-way intersection
setPlotBars(np_overlap)
setPlotPie(np_overlap)
setPlotVenn(np_overlap)
setPlotEuler(np_overlap)
setPlotHeatmap(np_overlap)

### prepare overlap sites for retrieving data from bigwigs
## determine desired width
boxplot(width(np_overlap))
width_q75 = np_overlap %>% width() %>% quantile(., .75)
width_q75 = ceiling(width_q75 / 100) * 100

## apply width + a little extra
fixed_overlap = centerFixedSizeGRanges(np_overlap, width_q75 + 400)
sample_rate = .01
qgr = sample(fixed_overlap, length(fixed_overlap) * sample_rate)


### fetch profile data from bigwigs for all sites
fetchWindowedBigwigList(bw_files = bw_files, qgr = qgr, win_size =  50, bw_variable_name = "cell_line")

import.bed("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10AT1_CTCF_pooled_FE.bw")
import.bed("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10AT1_CTCF_pooled_peaks_passIDR.05.narrowPeak.gz")
import.bed15("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10AT1_CTCF_pooled_peaks_passIDR.05.narrowPeak.gz")
curl::curl_download("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10AT1_CTCF_pooled_peaks_passIDR.05.narrowPeak.gz", destfile = "tmp.gz")
import.
