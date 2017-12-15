###load peak sets
library(magrittr)
library(GenomicRanges)
# library(peakvisr)

bw_urls = c("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77772/suppl/GSE77772_MCF7_H3K27ME3_logFE.bw",
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69377/suppl/GSE69377_MCF7_H3K4ME3_logFE.bw")


qgr = read.table("/slipstream/home/joeboyd/H7_bivalent_enst_promoter_ext2500_5573genes.bed")
colnames(qgr) = c("seqnames", "start", "end", "id", "zero", "strand")
qgr = GRanges(qgr)
end(qgr) = end(qgr) - 1
start(qgr) = start(qgr) - 10000
end(qgr) = end(qgr) + 10000


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
ggplot(all_bw_dt) + geom_raster(aes(x = x, y = id, fill = FE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)

all_bw_dt[, capFE := FE / max(FE), by = sample]
ggplot(all_bw_dt) + geom_raster(aes(x = x, y = id, fill = capFE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)

all_bw_dt[, lgFE := log2(FE + 1)]
ggplot(all_bw_dt) + geom_raster(aes(x = x, y = id, fill = lgFE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)



plot_dt = copy(all_bw_dt)

# dc_dt = dcast(plot_dt[abs(x) < 200], id ~ sample + x, value.var = "capFE")
dc_dt = dcast(plot_dt[abs(x) < 1000], id ~ sample + x, value.var = "lgFE")
dc_mat = as.matrix(dc_dt[,-1])
rownames(dc_mat) = dc_dt$id
rclusters = clustering_kmeans_nested_hclust(dc_mat, nclust = 6)
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
ggplot(plot_dt) + geom_raster(aes(x = x, y = id, fill = zFE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_gradient2(low = "black", mid = "black", high = "yellow") + facet_grid(. ~ sample)

pdf("tmp1.pdf", width = 4, height = 40)
p = ggplot(plot_dt) + geom_raster(aes(x = x, y = id, fill = zFE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_gradient2(low = "black", mid = "black", high = "yellow") + facet_grid(. ~ sample)
print(p)
dev.off()

pdf("tmp2.pdf")
for(grp in unique(rclusters$group)){
  tp = rclusters[group == grp, id]
  p = ggplot(plot_dt[id %in% tp]) + geom_raster(aes(x = x, y = id, fill = zFE)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
    scale_fill_gradient2(low = "black", mid = "black", high = "yellow") + facet_grid(. ~ sample)
  print(p)
}
dev.off()


grp = 3
tp = rclusters[group == grp, id]
dc_dt = dcast(plot_dt[abs(x) < 1000 & id %in% tp], id ~ sample + x, value.var = "zFE")
dc_mat = as.matrix(dc_dt[,-1])
rownames(dc_mat) = dc_dt$id
rclusters = clustering_kmeans_nested_hclust(dc_mat, nclust = 6)
plot_dt$id = factor(plot_dt$id, levels = rclusters$id)
scale_floor = .5
scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))
cap = 3
plot_dt[zFE < -cap, zFE := -cap]
plot_dt[zFE > cap, zFE := cap]
ggplot(plot_dt[id %in% tp]) + geom_raster(aes(x = x, y = id, fill = zFE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_gradient2(low = "black", mid = "black", high = "yellow") + facet_grid(. ~ sample)
