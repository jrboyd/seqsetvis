###load peak sets
library(magrittr)
library(GenomicRanges)
library(seqsetvis )
library(data.table)

# #Goal, find interesting clustering of histone modifications at gene promoters
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# # genes(TxDb.Hsapiens.UCSC.hg38.knownGene, filter = list(tx_chrom = "chr19"))
# tx = transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene, filter = list(tx_chrom = "chr19"))
# # exons(TxDb.Hsapiens.UCSC.hg38.knownGene, filter = list(tx_chrom = "chr19"))
# # cds(TxDb.Hsapiens.UCSC.hg38.knownGene, filter = list(tx_chrom = "chr19"))
#
# library(biomaRt)
# ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
# subset(listFilters(ensembl), grepl("UCSC", description))
# subset(listFilters(ensembl), grepl("ENSG", description))
# subset(listFilters(ensembl), grepl("ENST", description))
# subset(listFilters(ensembl), grepl("HGNC", description))
# subset(listFilters(ensembl), grepl("GO", description))
# go = "GO:0030154" #cell differentiation
# tx_meta =  getBM(mart = ensembl,
#                            attributes = c("ucsc", "ensembl_gene_id", "hgnc_symbol", "ensembl_transcript_id"),
#                            filters = c("ucsc"),
#                            values = list(tx$tx_name))
#
# tx_meta_cell_diff =  getBM(mart = ensembl,
#                  attributes = c("ucsc", "ensembl_gene_id", "hgnc_symbol", "ensembl_transcript_id"),
#                  filters = c("ucsc", "go"),
#                  values = list(tx$tx_name, go))
#
# tx_cell_diff = subset(tx, tx_name %in% tx_meta_cell_diff$ucsc)
# dt.tx_cell_diff = as.data.table(merge(tx_cell_diff, tx_meta_cell_diff, by.x = "tx_name", by.y = "ucsc"))
# dt.tx_cell_diff[strand == "+", end := start]
# dt.tx_cell_diff[strand == "-", start := end]
# ext = 2000L
# dt.tx_cell_diff[, c("start", "end") := .(start - ext, end + ext - 1L)]
# tx_cell_diff = GRanges(dt.tx_cell_diff)
#
#
#
# subset(listFilters(ensembl), grepl("UCSC", description))
# listDatasets(ensembl)
# mart = biomaRt::listMarts()

qgr = read.table("/slipstream/home/joeboyd/H7_bivalent_enst_promoter_ext2500_5573genes.bed")
colnames(qgr) = c("seqnames", "start", "end", "id", "zero", "strand")
qgr = GRanges(qgr)
end(qgr) = end(qgr) - 1
start(qgr) = start(qgr) - 10000
end(qgr) = end(qgr) + 10000

peak_files = dir("inst/extdata/peak_calls/full/", pattern = "narrow", full.names = T)
names(peak_files) = peak_files %>% basename %>%
  sub("_pooled_peak.+roadPeak", "", .) %>%
  sub("_pooled_peak.+arrowPeak", "", .)
peak_df = lapply(peak_files, read.table)
peak_gr = lapply(peak_df, function(x){
  colnames(x) = c("seqnames", "start", "end", "peak_id", "score", "dot", "FE", "pval", "qval", "rel_summit")
  GRanges(x)
})

# names(qgr) = qgr$id
peak_gr$base = qgr
peak_gr = rev(peak_gr)
peak_olap = overlapIntervalSets(peak_gr, use_first = T)
as.data.table(mcols(peak_olap))
wsize = 50
scale_floor = .3
scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))

bw_files = dir("inst/extdata/bigwigs/full", pattern = ".bw$", full.names = T)
names(bw_files) = bw_files %>% basename %>%
  sub("_pooled_F.+w", "", .)
bw_list = lapply(names(bw_files), function(nam){
  f = bw_files[nam]
  dt = fetchWindowedBigwig(bw_file = f, win_size = wsize, qgr = qgr)
  dt$sample = nam
  dt
})
all_bw_dt = rbindlist(bw_list)
all_bw_dt[strand == "-", x := -x]
all_bw_dt = all_bw_dt[x > -4000 & x < 12000]


all_bw_dt$id = factor(all_bw_dt$id)
# ggplot(all_bw_dt) + geom_raster(aes(x = x, y = id, fill = FE)) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
#   scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)
#
all_bw_dt[, capFE := FE / max(FE), by = sample]
# ggplot(all_bw_dt) + geom_raster(aes(x = x, y = id, fill = capFE)) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
#   scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)
#
all_bw_dt[, lgFE := log2(FE + 1)]
# ggplot(all_bw_dt) + geom_raster(aes(x = x, y = id, fill = lgFE)) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
#   scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)

# pdf('qbands.pdf')
regionSetPlotBandedQuantiles(all_bw_dt, y_ = "capFE", by_ = "sample",
                             symm_colors = T, hsv_hue_min = 0, hsv_hue_max = .7,
                             hsv_reverse = T)
regionSetPlotBandedQuantiles(all_bw_dt, y_ = "capFE", by_ = "sample", hsv_grayscale = T,
                             symm_colors = T, hsv_hue_min = 0, hsv_hue_max = .7,
                             hsv_reverse = T)
# dev.off()

# regionSetPlotBandedQuantiles(all_bw_dt, y_ = "capFE", by_ = "sample", symm_colors = T, hsv_hue_min = 0, hsv_hue_max = .7)
# regionSetPlotBandedQuantiles(all_bw_dt, y_ = "capFE", by_ = "sample", symm_colors = T, hsv_hue_min = 0, hsv_hue_max = .7, )
# regionSetPlotBandedQuantiles(all_bw_dt, y_ = "lgFE", by_ = "sample")
# regionSetPlotBandedQuantiles(all_bw_dt, y_ = "FE", by_ = "sample")


plot_dt = copy(all_bw_dt[sample == "H3K27ME3"])
plot_dt = copy(all_bw_dt)
# dc_dt = dcast(plot_dt[abs(x) < 200], id ~ sample + x, value.var = "capFE")
dc_dt = dcast(plot_dt[abs(x) < 1000], id ~ sample + x, value.var = "lgFE")
dc_mat = as.matrix(dc_dt[,-1])
rownames(dc_mat) = dc_dt$id
rclusters = clusteringKmeansNestedHclust(dc_mat, nclust = 6)
plot_dt$id = factor(plot_dt$id, levels = rclusters$id)

scale_floor = .1
scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))
ggplot(plot_dt) + geom_raster(aes(x = x, y = id, fill = capFE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)

layout(matrix(1:4, ncol = 2))
plot_dt[, hist(lgFE, main = sample), by = sample]

plot_dt[, zFE := (lgFE - mean(lgFE)) / sd(lgFE), by = sample]
plot_dt[, hist(zFE, main = sample), by = sample]

dc_dt = dcast(plot_dt[abs(x) < 1000], id ~ sample + x, value.var = "zFE")
dc_mat = as.matrix(dc_dt[,-1])
rownames(dc_mat) = dc_dt$id
rclusters = clustering_kmeans_nested_hclust(dc_mat, nclust = 6)
plot_dt$id = factor(plot_dt$id, levels = rclusters$id)
scale_floor = .5
scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))
cap = 3
plot_dt[zFE < -cap, zFE := -cap]
plot_dt[zFE > cap, zFE := cap]

regionSetPlotHeatmap(all_bw_dt)
regionSetPlotHeatmap(all_bw_dt, y_ = "zFE")
# ggplot(plot_dt) + geom_raster(aes(x = x, y = id, fill = zFE)) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
#   scale_fill_gradient2(low = "black", mid = "black", high = "yellow") + facet_grid(. ~ sample)
#
# pdf("tmp1.pdf", width = 4, height = 40)
# p = ggplot(plot_dt) + geom_raster(aes(x = x, y = id, fill = zFE)) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
#   scale_fill_gradient2(low = "black", mid = "black", high = "yellow") + facet_grid(. ~ sample)
# print(p)
# dev.off()
#
# pdf("tmp2.pdf")
# for(grp in unique(rclusters$group)){
#   tp = rclusters[group == grp, id]
#   p = ggplot(plot_dt[id %in% tp]) + geom_raster(aes(x = x, y = id, fill = zFE)) +
#     theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
#     scale_fill_gradient2(low = "black", mid = "black", high = "yellow") + facet_grid(. ~ sample)
#   print(p)
# }
# dev.off()
#
#
# grp = 3
# tp = rclusters[group == grp, id]
# dc_dt = dcast(plot_dt[abs(x) < 1000 & id %in% tp], id ~ sample + x, value.var = "zFE")
# dc_mat = as.matrix(dc_dt[,-1])
# rownames(dc_mat) = dc_dt$id
# rclusters = clustering_kmeans_nested_hclust(dc_mat, nclust = 6)
# plot_dt$id = factor(plot_dt$id, levels = rclusters$id)
# scale_floor = .5
# scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))
# cap = 3
# plot_dt[zFE < -cap, zFE := -cap]
# plot_dt[zFE > cap, zFE := cap]
# ggplot(plot_dt[id %in% tp]) + geom_raster(aes(x = x, y = id, fill = zFE)) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
#   scale_fill_gradient2(low = "black", mid = "black", high = "yellow") + facet_grid(. ~ sample)
#
