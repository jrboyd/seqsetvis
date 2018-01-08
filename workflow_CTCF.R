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


bw_urls = c("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10A_CTCF_pooled_FE.bw",
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10AT1_CTCF_pooled_FE.bw",
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98551/suppl/GSE98551_MCF10CA1_CTCF_pooled_FE.bw")

# bw_urls = https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE98551&format=file&file=GSE98551%5FMCF10AT1%5FCTCF%5Fpooled%5FFE%2Ebw
names(bw_urls) = bw_urls %>%
    basename() %>%
    sub("GSE98551_", "", x = .) %>%
    sub("_pooled_FE.bw", "", x = .)

bw_files = paste0(names(bw_urls), "_FE.bw")
names(bw_files) = names(bw_urls)

### load peak calls

#from: https://charlesjb.github.io/How_to_import_narrowPeak/
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
np_grs <- lapply(np_files, function(f){
    import(f, format = "BED",
           extraCols = extraCols_narrowPeak)
})

np_grs= lapply(np_grs, function(gr){
    subset(gr, signalValue > 15)
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
fixed_overlap = centerFixedSizeGRanges(np_overlap, width_q75 * 3)
sample_rate = 1
qgr = sample(fixed_overlap, length(fixed_overlap) * sample_rate)

### download bigwigs for local access
for(i in seq_along(bw_urls)){
    if(!file.exists(bw_files[i])){
        curl::curl_download(bw_urls[i], destfile = bw_files[i])
    }
}
### fetch profile data from bigwigs for all sites
bw_dt = fetchWindowedBigwigList(bw_files = bw_files, qgr = qgr, win_size =  50, bw_variable_name = "cell_line")

h1 = regionSetPlotHeatmap(bw_dt[abs(x) <= width_q75/2], facet_ = "cell_line", nclust = 12)

cbw_dt = centerAtMax(bw_dt, y_ = "FE", by = "id", view_size = width_q75, check_by_dupes = F)
h2 = regionSetPlotHeatmap(cbw_dt, facet_ = "cell_line", nclust = 12)

pdf("h.pdf", width = 4, height = 16)
cowplot::plot_grid(h1, h2, ncol = 1)
dev.off()
q98 = function(x)quantile(x, .98)
p = regionSetPlotScatter(bw_dt = bw_dt,
                         x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF",
                         xy_variable = "cell_line", value_function = q98) +
    labs(title = "98th quantile in regions")
# p
# p + scale_alpha_continuous(range = .03)
p1 = p + geom_density_2d() + scale_alpha_continuous(range = .03)
# p + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_alpha_continuous(range = .03)

p = regionSetPlotScatter(bw_dt = bw_dt, x_name = "MCF10A_CTCF", y_name = "MCF10CA1_CTCF", xy_variable = "cell_line", value_function = q98) + labs(title = "98th quantile in regions")
# p
# p + scale_alpha_continuous(range = .03)
p2 = p + geom_density_2d() + scale_alpha_continuous(range = .03)
# p + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_alpha_continuous(range = .03)

p = regionSetPlotScatter(bw_dt = bw_dt, x_name = "MCF10AT1_CTCF", y_name = "MCF10CA1_CTCF", xy_variable = "cell_line", value_function = q98) + labs(title = "98th quantile in regions")
# p
# p + scale_alpha_continuous(range = .03)
p3 = p + geom_density_2d() + scale_alpha_continuous(range = .03)
# p + stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_alpha_continuous(range = .03)
cowplot::plot_grid(p1, p2, p3, nrow = 1)
