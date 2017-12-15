###load peak sets
library(magrittr)
library(GenomicRanges)
# library(peakvisr)

peak_files = dir("inst/extdata/peak_calls/full", pattern = "narrowPeak", full.names = T)
names(peak_files) = peak_files %>% basename %>%
  sub("_pooled_peak.+roadPeak", "", .) %>%
  sub("_pooled_peak.+arrowPeak", "", .)

peak_df = lapply(peak_files, read.table)
peak_gr = lapply(peak_df, function(x){
  colnames(x) = c("seqnames", "start", "end", "peak_id", "score", "dot", "FE", "pval", "qval", "rel_summit")
  GRanges(x)
})

layout(matrix(1:4, ncol = 2))
td = names(peak_gr)
names(td) = td
cut_type = "pval"
peak_gr_raw = peak_gr
peak_gr = lapply(td, function(nam){
  x = as.data.frame(peak_gr[[nam]])
  n = 100
  cuts = 0:n / n * max(x[[cut_type]])
  counts = sapply(cuts, function(y)sum(x[[cut_type]] >= y))
  plot(cuts, counts, type = "l")
  dx = rowMeans(cbind(cuts[-length(counts)], cuts[-1]))
  dy = counts[-length(counts)] - counts[-1]
  lines(dx, dy, col = "red")
  xm = dx[which(dy == max(dy))[1]]
  lines(rep(xm, 2), c(0, max(counts)), lty = 2, col = "red")
  title(nam)
  newx = subset(x, get(cut_type) >= xm)
  lines(c(0, max(cuts)), rep(nrow(newx), 2), lty = 2, col = "red")
  text(mean(c(0, max(cuts))), nrow(newx), labels = nrow(newx), adj = c(.5,-.5), col = "red")
  ddx = rowMeans(cbind(dx[-length(dx)], dx[-1]))
  ddy = dy[-length(dy)] - dy[-1]
  lines(ddx, ddy, col = "blue")
  dxm = ddx[which(ddy == max(ddy))[1]]
  lines(rep(dxm, 2), c(0, max(counts)), lty = 2, col = "blue")
  newdx = subset(x, get(cut_type) >= dxm)
  lines(c(0, max(cuts)), rep(nrow(newdx), 2), lty = 2, col = "blue")
  text(mean(c(0, max(cuts))), nrow(newdx), labels = nrow(newdx), adj = c(.5,1.5), col = "blue")
  return(GRanges(newdx))
})

peak_tab = overlapIntervalSets(peak_gr_raw[1:4])
peak_memb = elementMetadata(peak_tab)

cols = RColorBrewer::brewer.pal(4, "Dark2")
###venn diagrams
ggVenn(peak_memb[,1:2], circle.col = cols[1:2])
ggVenn(peak_memb[,3:4], circle.col = cols[3:4])
ggVenn(peak_memb[,c(1,3)], circle.col = cols[c(1, 3)])
ggVenn(peak_memb[,c(2,4)], circle.col = cols[c(2, 4)])

###barplots
ggBars(peak_memb) + scale_fill_manual(values = cols)

###pie
ggPie(peak_memb) + scale_fill_manual(values = cols)

###euler
ggEuler(peak_memb, line_width = 2)

bw_files = dir("inst/extdata/bigwigs/full", full.names = T)
names(bw_files) = bw_files %>% basename %>%
  sub("_pooled_F.+w", "", .)

names(peak_tab) = paste0("region_", seq_along(peak_tab))

# runx_10a = subset(peak_tab, H7_H3K4ME3)
qgr = centerFixedSizeGRanges(peak_tab, 20000)
qgr = read.table("/slipstream/home/joeboyd/H7_bivalent_enst_promoter_ext2500_5573genes.bed")

wsize = 50

bw_list = lapply(names(bw_files), function(nam){
  f = bw_files[nam]
  dt = fetchWindowedBigwig(bw_file = f, win_size = wsize, qgr = qgr)
  dt$sample = nam
  dt
})
all_bw_dt = rbindlist(bw_list)
bw_dt = all_bw_dt[sample == "MCF7_H3K4ME3"]
m_dt = bw_dt[, .(mFE = mean(FE)), by = .( x)]
spline_dt = bw_dt[, spline(x, FE), by = id]
# m_dt = m_dt[x %% 10 == 0]
ggplot(m_dt) + geom_line(aes(x = x, y = mFE))
ggplot(applySpline(m_dt, y_= "mFE")) + geom_line(aes(x = x, y = y))

spline_dt = applySpline(bw_dt, y_ = "FE", by_ = "id")
mcf7_k4_pos = names(subset(qgr, H7_H3K4ME3 & H7_H3K27ME3))
# mcf7_k4_pos = names(subset(qgr, MCF7_H3K4ME3 & MCF7_H3K27ME3))
tp2 = sample(mcf7_k4_pos, size = min(10, length(mcf7_k4_pos)))
ggplot(spline_dt[id %in% tp2]) + geom_line(aes(x = x, y = y)) + facet_grid(id ~ ., scales = "free")

o_dt = spline_dt[, .(ymax = max(y)), by = id]
lev = o_dt[order(ymax), id]
spline_dt$id = factor(spline_dt$id, levels = lev)

# ggplot(spline_dt[id %in% tp]) + geom_raster(aes(x = x, y = id, fill = y))
bw_dt$id = factor(bw_dt$id, levels = lev)
tp = mcf7_k4_pos

cmax_dt = centerAtMax(bw_dt, final_size = 4000)
o_dt = cmax_dt[, .(ymax = max(FE)), by = id]
lev = o_dt[order(ymax), id]
cmax_dt$id = factor(cmax_dt$id, levels = lev)

ggplot(cmax_dt[id %in% tp]) + geom_raster(aes(x = x, y = id, fill = FE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_distiller(type = "div", palette = "Spectral")


bw_dt[id %in% tp & x >= -2000 & x <= 2000]
all_bw_dt$id = factor(all_bw_dt$id, levels = lev)
all_bw_dt$sample = factor(all_bw_dt$sample, levels = c("H7_H3K4ME3", "MCF7_H3K4ME3", "H7_H3K27ME3", "MCF7_H3K27ME3"))
plot_dt = copy(all_bw_dt)


plot_dt[, lgFE := log2(FE)]
plot_dt[is.infinite(lgFE), lgFE := 0]

plot_dt[, hist(lgFE, main = sample), by = sample]

plot_dt[, zFE := (lgFE - mean(lgFE)) / sd(lgFE), ]
plot_dt[, hist(zFE, main = sample), by = sample]

dc_dt = dcast(plot_dt[abs(x) < 200], id ~ sample + x, value.var = "zFE")
dc_mat = as.matrix(dc_dt[,-1])
rownames(dc_mat) = dc_dt$id
rclusters = clustering_kmeans_nested_hclust(dc_mat, nclust = 6)
plot_dt$id = factor(plot_dt$id, levels = rclusters$id)
scale_floor = .5
scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))
cap = 3
plot_dt[zFE < -cap, zFE := -cap]
plot_dt[zFE > cap, zFE := cap]
ggplot(plot_dt[id %in% tp & abs(x) < 1000]) + geom_raster(aes(x = x, y = id, fill = zFE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)


# dc_dt = dcast(plot_dt, id ~ sample + x, value.var = "FE")
# dc_mat = as.matrix(dc_dt[,-1])
# rownames(dc_mat) = dc_dt$id
# km = kmeans(dc_mat, 6)
# id_lev = names(km$cluster)[order(km$cluster)]
# plot_dt$id = factor(plot_dt$id, levels = id_lev)
#plot_dt = plot_dt[sample %in% c("MCF7_H3K4ME3", "MCF7_H3K27ME3")]
plot_dt[, qFE := FE / quantile(FE, .99), by = .(sample) ]
plot_dt[qFE > 1, qFE := 1]
plot_dt = plot_dt[id %in% tp]
dc_dt = dcast(plot_dt, id ~ sample + x, value.var = "qFE")
dc_mat = as.matrix(dc_dt[,-1])
rownames(dc_mat) = dc_dt$id
km = kmeans(dc_mat, 6)
id_lev = names(km$cluster)[order(km$cluster)]
plot_dt$id = factor(plot_dt$id, levels = id_lev)

# plot_dt = plot_dt[sample %in% c("H7_H3K4ME3", "H7_H3K27ME3")]
cmax_plot_dt = centerAtMax(plot_dt, final_size = 10000, trim_to_valid = T)
cmax_plot_dt[, zFE := FE - mean, by = sample]
cmax_plot_dt[FE < 8, FE := 8]
scale_floor = .5
scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))
# ggplot(cmax_plot_dt[id %in% tp & x >= -2000 & x <= 2000]) + geom_raster(aes(x = x, y = id, fill = qFE)) +
p = ggplot(cmax_plot_dt) + geom_raster(aes(x = x, y = id, fill = qFE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)
p
# p + scale_fill_gradientn(colors = c("slateblue1", "slateblue1", "orange", "red"))

dc_dt = dcast(plot_dt[sample %in% c("MCF7_H3K4ME3", "MCF7_H3K27ME3") & abs(x) < 200], id ~ sample + x, value.var = "qFE")
dc_mat = as.matrix(dc_dt[,-1])
rownames(dc_mat) = dc_dt$id
rclusters = clustering_kmeans_nested_hclust(dc_mat, nclust = 6)
cmax_plot_dt$id = factor(cmax_plot_dt$id, levels = rclusters$id)

scale_floor = .1
scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))
cmax_plot_dt[FE > 15, FE := 15]
p = ggplot(cmax_plot_dt) + geom_raster(aes(x = x, y = id, fill = FE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)
p

scale_floor = .2
scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))
p = ggplot(cmax_plot_dt) + geom_raster(aes(x = x, y = id, fill = qFE)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank()) +
  scale_fill_distiller(type = "div", palette = "Spectral", values = scale_vals) + facet_grid(. ~ sample)
p



