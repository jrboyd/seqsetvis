library(data.table)
library(GenomicRanges)
library(pbapply)
library(ggplot2)
library(plotly)
load('gencode.v27.full.save')

#remove undesirable gene_types
# types = unique(dt$gene_type)
# pseudo_types = types[grepl("pseu", types)]
# dt = dt[!gene_type %in% pseudo_types]
# 
# types = unique(dt$gene_type)
# pseudo_types = types[grepl("^IG", types) | grepl("^TR", types) | grepl("^Mt", types)]
# dt = dt[!gene_type %in% pseudo_types]

dt$keep = F
nc_types = unique(dt$gene_type)
nc_types = c("antisense_RNA", "lincRNA", "sense_intronic",
             "bidirectional_promoter_lncRNA", "macro_lncRNA", "processed_transcript")
# write.table(nc_types, file = "nc_types.txt", quote = F, row.names = F, col.names = F)
dt[gene_type %in% c(nc_types, "protein_coding"), keep := T]
dt = dt[keep == T]

# dt[gene_type %in% c(nc_types, "protein_coding"), keep := T]
types = unique(dt$gene_type)

dt[strand == "+", end := start]
dt[strand == "-", start := end]
ext = 2000
dt[, start := start - ext]
dt[, end := end + ext]

gene_ref = copy(dt)
ref_report = copy(dt)
ref_gr = GRanges(dt)

mda_np_files = dir("/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/", 
                   pattern = "[12]_pooled", full.names = T)
mda_np_files = sapply(mda_np_files, function(x)dir(x, pattern = "passIDR.05.narrowPeak", full.names = T))
mda_grs = lapply(mda_np_files, function(x){
  dat = read.table(x)
  colnames(dat)= c("seqnames", "start", "end", "peak_id", "a", "b", "FE", "p", "q", "w")
  dat = subset(dat, FE > 6)
  GRanges(dat)
})
names(mda_grs) = sapply(strsplit(basename((mda_np_files)), "_"), function(x)x[2])

# names(mda_grs) = sapply(strsplit(basename((mda_np_files)), "_"), function(x)x[2])

hidden = lapply(names(mda_grs), function(nam){
  cn = paste("MDA231", nam, sep = "_")
  hits = unique(queryHits(findOverlaps(ref_gr, mda_grs[[nam]])))
  ref_report[[cn]] <<- F
  ref_report[[cn]][hits] <<- T
})

gene_t = "protein_coding"
layout(matrix(1:2))
limma::vennDiagram(ref_report[gene_type == gene_t, c("MDA231_Runx1", "MDA231_Runx2"), with = T], 
                   main = paste("MDA231 at", nrow(ref_report[gene_type == gene_t]), gene_t), names = names(mda_grs))
limma::vennDiagram(ref_report[gene_type != gene_t, c("MDA231_Runx1", "MDA231_Runx2"), with = T], 
                   main = paste("MDA231 at", nrow(ref_report[gene_type != gene_t]), "non-coding"), names = names(mda_grs))

to_prof_gr = GRanges(ref_report[MDA231_Runx1 | MDA231_Runx2])

runx_bws = c("RUNX1" = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MDA231_MK_runx/MDA231_Runx1_pooled_FE.bw",
             "RUNX2" = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MDA231_MK_runx/MDA231_Runx2_pooled_FE.bw")

load("JG_profs.2kb.save")

prof = copy(all_res)
prof = prof[, .(xrel = x - x[which(x == min(x))], x, y, group, index_o, index_se), by = .(gene_id)]
prof = prof[, x_tss_dist := xrel - floor(max(xrel) / 2)]
# prof = prof[, .(id = id - id[which(id == min(id))], y, group), by = .(gene_id, index_o, index_se)]

#problem with summit distance losing index_se shift
prof = prof[, .(x_summit_dist = xrel - xrel[which(y == max(y[xrel >= 2000 & xrel <= 6000]))], xrel, x_tss_dist, x, y, group, index_o, index_se), by = .(gene_id)]
prof = prof[x_summit_dist >= -2000 & x_summit_dist <= 2000]
table(prof$gene_id)

# setkey(gene_ref, gids)
setkey(ref_report, gene_id)
# prof[, c("root_id", "ver_id") := tstrsplit(gene_id, split = "\\.")]
prof$gene_type = ref_report[prof$gene_id]$gene_type
prof$seqnames = ref_report[prof$gene_id]$seqnames
prof[gene_type != "protein_coding", gene_type := "non_coding"]



# ggplot(prof[gene_id == "ENSG00000188976.10",]) + geom_line(aes(x = x_summit_dist, y = y, color = group))


#profiles around summits, by gene_type
summit_prof = prof[, .(y = mean(y)),by = .(x_summit_dist, group, index_se, gene_type)]
summit_prof = summit_prof[order(index_se)][order(x_summit_dist)][order(group)]
summit_prof[index_se == 2, x_summit_dist := x_summit_dist + 40]
ggplot(summit_prof) + geom_line(aes(x = x_summit_dist, y = y, color = group)) +
 facet_grid(. ~ gene_type) + labs(x = "distance to peak summit", y = "FE")
#profiles around summits, by chrm and gene_type
pdf("by_chrm.pdf", width = 4, height = 10)
chrsummit_prof = prof[, .(y = mean(y)),by = .(x_summit_dist, group, index_se, gene_type, seqnames)]
chrsummit_prof = chrsummit_prof[order(index_se)][order(x_summit_dist)]
ggplot(chrsummit_prof) + geom_line(aes(x = x_summit_dist, y = y, color = group)) +
  facet_grid(seqnames ~ gene_type) + labs(x = "distance to peak summit", y = "FE")
dev.off()

tss_prof = prof[, .(y = mean(y)),by = .(x_tss_dist, group, index_se, gene_type)]
ggplot(tss_prof) + geom_line(aes(x = x_tss_dist, y = y, color = group)) + 
  facet_grid(. ~ gene_type) + labs(x = "distance to tss", y = "FE")

ym = mean(prof$y) + 2*sd(prof$y)
ug = unique(prof$gene_id)
setkey(gene_ref, gene_id)
prof$gene_name = gene_ref[prof$gene_id]$gene_name


for(i in 1:2){#ceiling(length(ug) / 6)){
  js = -5:0+6*i
  p = ggplot(prof[gene_id %in% ug[js]]) + 
    geom_line(aes(x = x_tss_dist, y = y, color = group)) + 
    facet_grid(gene_name ~ group) + annotate("line", x = c(-2000,2000), y = rep(ym, 2)) + 
    labs(x = "relative position to tss", y = "FE") + 
    theme(strip.text.y = element_text(size = 8, colour = "black", angle = 0))
  print(p)
}
for(i in 1:2){#ceiling(length(ug) / 6)){
  js = -5:0+6*i
  p = ggplot(prof[gene_id %in% ug[js]]) + 
    geom_line(aes(x = x_summit_dist, y = y, color = group)) + 
    facet_grid(gene_name ~ group) + annotate("line", x = c(-2000,2000), y = rep(ym, 2)) + 
    labs(x = "relative position to peak summit", y = "FE") + 
    theme(strip.text.y = element_text(size = 8, colour = "black", angle = 0))
  print(p)
}


scat_prof = prof[x_summit_dist == 0 & index_se == 2]
scat_prof = dcast(scat_prof, formula = x_summit_dist + gene_type + gene_id + gene_name + seqnames  ~ group, value.var = "y")
# p = ggplot(scat_prof[sample(nrow(scat_prof))], aes(x = RUNX1, y = RUNX2, color = gene_type)) + geom_point(aes(label1 = gene_name, label2 = seqnames)) + labs(title = "FE at strongest peak within 2kb of tss")
# ggplotly(p, tooltip = c("label1", "label2"))
