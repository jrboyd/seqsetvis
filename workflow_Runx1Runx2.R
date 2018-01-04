library(seqsetvis)
library(rtracklayer)
library(magrittr)
library(data.table)

###load reference
source("../jrb_R_scripts/parse_gtf.dt.R")
ref_dt = parse_gtf.dt("/slipstream/home/joeboyd/gencode.v25.annotation.gtf", additional_attribs = "gene_type")
ref_dt = ref_dt[gene_type %in% c("lincRNA", "antisense_RNA", "sense_intronic") & seqnames != "chrY"]
ref_dt[strand == "+", end := start]
ref_dt[strand == "-", start := end]
ext = 30*10^3
ref_dt[,start := start - ext]
ref_dt[,end := end + ext - 1]
ref_gr = GRanges(ref_dt)

### set file paths
np_files = c("/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx1_pooled/MDA231_Runx1_pooled_peaks_passIDR.05_FEgt6.bed",
             "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled_peaks_passIDR.05_FEgt6.bed")
names(np_files) = np_files %>%
    basename() %>%
    sub("GSE98551_", "", x = .) %>%
    sub("_pooled_peaks_passIDR.05_FEgt6.bed", "", x = .)


bw_files = c("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MDA231_MK_runx/MDA231_Runx1_pooled_FE.bw",
             "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MDA231_MK_runx/MDA231_Runx2_pooled_FE.bw")

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
    import(f, format = "BED")
})

### perform 3-way intersection
# np_overlap = overlapIntervalSets(np_grs)
np_overlap = overlapIntervalSets(append(list("lnc" = ref_gr), np_grs), use_first = T)
np_overlap$no_hit = NULL

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
fixed_overlap = centerFixedSizeGRanges(np_overlap, width_q75 + 2000)
sample_rate = 1
qgr = sample(fixed_overlap, length(fixed_overlap) * sample_rate)


### fetch profile data from bigwigs for all sites
bw_dt = fetchWindowedBigwigList(bw_files = bw_files, qgr = qgr, win_size =  50, bw_variable_name = "sample")

cdt = centerAtMax(bw_dt, view_size = 28000, y_ = "FE", by_ = c("id"), check_by_dupes = F, trim_to_valid = T)
cdt[, zscore := (FE - mean(FE)) / sd(FE), by = sample]

ggplot(cdt[id %in% c("5588", "4619", "7828", "2447")], mapping = aes(x = x, y = FE, color = sample)) +
    geom_line() + facet_grid(id ~ .)

cap = 10
cdt[zscore < cap, zscore := cap]
ggplot(cdt[id %in% c("5588", "4619", "7828", "2447")], mapping = aes(x = x, y = zscore, color = sample)) +
    geom_line() + facet_grid(id ~ .)


plot_dt = copy(cdt)

regionSetPlotHeatmap(plot_dt)
regionSetPlotBandedQuantiles(plot_dt, by_ = "sample", hsv_symmetric = T, hsv_reverse = T)

memb = mcols(qgr)
memb$id = names(qgr)
memb$plotting_group = "none"
memb[memb$MDA231_Runx1,]$plotting_group = "MDA231_Runx1"
memb[memb$MDA231_Runx2,]$plotting_group = "MDA231_Runx2"
memb[memb$MDA231_Runx1 & memb$MDA231_Runx2,]$plotting_group = "both"
memb$plotting_group = factor(memb$plotting_group)

regionSetPlotScatter(plot_dt, x_name = "MDA231_Runx1", y_name = "MDA231_Runx2", plotting_group = as.data.table(memb[,3:4]))
