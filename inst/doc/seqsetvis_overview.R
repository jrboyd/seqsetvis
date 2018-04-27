## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    dpi=60
)

## ----bioc install, eval=FALSE----------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("seqsetvis")

## ----load seqsetvis, message=FALSE-----------------------------------------
library(seqsetvis)

## ----load optional libs, message = FALSE-----------------------------------
library(GenomicRanges)
library(data.table)
library(cowplot)
theme_set(cowplot::theme_cowplot())

## ----overlap basic---------------------------------------------------------
olaps = ssvOverlapIntervalSets(CTCF_in_10a_narrowPeak_grs)
head(olaps)

## ----overlap GRangesList---------------------------------------------------
olaps_fromGRangesList = ssvOverlapIntervalSets(
    GenomicRanges::GRangesList(CTCF_in_10a_narrowPeak_grs))

## ----ssvMakeMembTable basic------------------------------------------------
head(ssvMakeMembTable(olaps))

## ----ssvMakeMembTable numeric----------------------------------------------
my_set_list = list(1:3, 2:3, 3:6)
ssvMakeMembTable(my_set_list)

## ----ssvMakeMembTable named numeric----------------------------------------
names(my_set_list) = c("first", "second", "third")
ssvMakeMembTable(my_set_list)

## ----ssvMakeMembTable character--------------------------------------------
my_set_list_char = lapply(my_set_list, function(x)letters[x])
ssvMakeMembTable(my_set_list_char)

## ----barplot, fig.width=4, fig.height=3------------------------------------
ssvFeatureBars(olaps)

## ----pie, fig.width=5, fig.height=3----------------------------------------
ssvFeaturePie(olaps)

## ----venn, fig.width=4, fig.height=3---------------------------------------
ssvFeatureVenn(olaps)

## ----euler, fig.width=4, fig.height=3--------------------------------------
ssvFeatureEuler(olaps)

## ----binary heatmap, fig.width=3, fig.height=4-----------------------------
ssvFeatureBinaryHeatmap(olaps)

## ----fetchWindowedBigwig, eval = FALSE-------------------------------------
#  bigwig_files = c(
#      system.file("extdata", "MCF10A_CTCF_FE_random100.bw",
#                  package = "seqsetvis"),
#      system.file("extdata", "MCF10AT1_CTCF_FE_random100.bw",
#                  package = "seqsetvis"),
#      system.file("extdata", "MCF10CA1_CTCF_FE_random100.bw",
#                  package = "seqsetvis")
#  )
#  names(bigwig_files) = sub("_FE_random100.bw", "", basename(bigwig_files))
#  # names(bigwig_files) = letters[1:3]
#  olap_gr = CTCF_in_10a_overlaps_gr
#  target_size = quantile(width(olap_gr), .75)
#  window_size = 50
#  target_size = round(target_size / window_size) * window_size
#  olap_gr = resize(olap_gr, target_size, fix = "center")
#  bw_gr = fetchWindowedBigwigList(bigwig_files, olap_gr, win_size = window_size)
#  
#  bw_gr

## --------------------------------------------------------------------------
olap_gr = CTCF_in_10a_overlaps_gr
bw_gr = CTCF_in_10a_profiles_gr

## ----factorize-------------------------------------------------------------
olap_groups = ssvFactorizeMembTable(mcols(olap_gr))

## ----lineplot basic, fig.width=6, fig.height=2.5---------------------------
# facet labels will display better if split into multiple lines
bw_gr$facet_label = sub("_", "\n", bw_gr$sample)
ssvSignalLineplot(bw_data = subset(bw_gr, id %in% 1:12), facet_ = "facet_label")

## ----lineplot region facet, fig.width=5, fig.height=3----------------------
ssvSignalLineplot(bw_data = subset(bw_gr, id %in% 1:4), facet_ = "id")

## ----lineplot aggregated, fig.width=5, fig.height=2------------------------
ssvSignalLineplotAgg(bw_data = bw_gr)

## ----lineplot aggregated smoothed, fig.width=5, fig.height=2---------------
ssvSignalLineplotAgg(bw_data = bw_gr, spline_n = 10)

## ----lineplot--------------------------------------------------------------
# append set info, modify aggregation group_ and add facet
olap_2groups = ssvFactorizeMembTable(ssvMakeMembTable(olap_gr)[, 1:2])
grouped_gr = GRanges(merge(bw_gr, olap_2groups))
grouped_gr = subset(grouped_gr, sample %in% c("MCF10A_CTCF", "MCF10AT1_CTCF"))
ssvSignalLineplotAgg(bw_data = grouped_gr, spline_n = 10,
                     group_ = c("sample", "group")) + 
    facet_wrap("group", ncol = 2) +
    labs(title = "Aggregated by peak call set", y = "FE", x = "bp from center")

## ----scatterplot basic, fig.width=3, fig.height=3--------------------------
ssvSignalScatterplot(bw_gr, x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF")

## ----scatterplot all sets, fig.width=8, fig.height=3-----------------------
ssvSignalScatterplot(bw_gr, x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF", 
                     color_table = olap_groups)

## ----scatterplot 2 sets, fig.width=6, fig.height=3-------------------------
# by subsetting the matrix returned by ssvMakeMembTable() we have a lot of 
# control over the coloring.
olap_2groups = ssvFactorizeMembTable(ssvMakeMembTable(olap_gr)[, 1:2])
ssvSignalScatterplot(bw_gr, x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF", 
                     color_table = olap_2groups)

## ----outside group, fig.width=5, fig.height=3------------------------------
olap_OutGroups = ssvFactorizeMembTable(
    ssvMakeMembTable(olap_gr)[, 3, drop = FALSE])
ssvSignalScatterplot(bw_gr, 
                     x_name = "MCF10A_CTCF", 
                     y_name = "MCF10AT1_CTCF", 
                     color_table = olap_OutGroups)

## ----scatterplot facet, fig.width=6, fig.height=4--------------------------
#tweaking group description will clean up plot labels a lot
olap_groups$group = gsub("_CTCF", "", olap_groups$group)
olap_groups$group = gsub(" & ", "\n", olap_groups$group)
ssvSignalScatterplot(bw_gr, x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF", 
                     color_table = olap_groups) + 
    facet_wrap("group") + guides(color = "none") +
    theme_linedraw() 

## ----banded quantiles------------------------------------------------------
ssvSignalBandedQuantiles(bw_gr, by_ = "sample", hsv_grayscale = TRUE, 
                         hsv_symmetric = TRUE, 
                         hsv_reverse = TRUE)

## ----heatmap basic, message=FALSE, fig.width=5-----------------------------
ssvSignalHeatmap(bw_gr, nclust = 3, facet_ = "facet_label") 

## ----heatmap perform pre-clustering----------------------------------------
bw_clust = ssvSignalClustering(bw_gr, nclust = 3)
bw_clust

## ----heatmap cluster selection---------------------------------------------
subset(bw_clust, cluster_id == 3)

## ----heatmap use pre-cluster, message=FALSE, fig.width=5-------------------
ssvSignalHeatmap(bw_clust, facet_ = "facet_label")

## ----setup np_files bw_files-----------------------------------------------
pkgdata_path = system.file("extdata", 
                           package = "seqsetvis")
cache_path = paste0(pkgdata_path, "/.cache")
# the next line is enough to initialize the cache
# BiocFileCache(cache = cache_path)
use_full_data = dir.exists(cache_path) & require(BiocFileCache)
if(use_full_data){
    library(BiocFileCache)
    ssv_bfc = BiocFileCache(cache = cache_path)
    bw_files = vapply(seq_along(CTCF_in_10a_bigWig_urls), function(i){
        rname = paste(names(CTCF_in_10a_bigWig_urls)[i], 
                      "bigwig", 
                      sep = ",")
        fpath = CTCF_in_10a_bigWig_urls[i]
        #bfcrpath calls bfcadd() if necessary and returns file path
        bfcrpath(ssv_bfc, rname = rname, fpath = fpath)
    }, "char")
    names(bw_files) = names(CTCF_in_10a_bigWig_urls)
    
    np_files = vapply(seq_along(CTCF_in_10a_narrowPeak_urls), function(i){
        rname = paste(names(CTCF_in_10a_narrowPeak_urls)[i], 
                      "narrowPeak", 
                      sep = ",")
        fpath = CTCF_in_10a_narrowPeak_urls[i]
        #bfcrpath calls bfcadd() if necessary and returns file path
        bfcrpath(ssv_bfc, rname = rname, fpath = fpath)
    }, "a")
    names(np_files) = names(CTCF_in_10a_narrowPeak_urls)
}else{
    bw_files = vapply(c("MCF10A_CTCF", "MCF10AT1_CTCF", "MCF10CA1_CTCF"), 
    function(x){
        system.file("extdata", paste0(x, "_FE_random100.bw"), 
                    package = "seqsetvis")
    }, "char")
    # set filepaths
    np_files = c(
        system.file("extdata", "MCF10A_CTCF_random100.narrowPeak", 
                    package = "seqsetvis"),
        system.file("extdata", "MCF10AT1_CTCF_random100.narrowPeak", 
                    package = "seqsetvis"),
        system.file("extdata", "MCF10CA1_CTCF_random100.narrowPeak", 
                    package = "seqsetvis")
    )
    
    names(np_files) = sub("_random100.narrowPeak", "", 
                          x = basename(np_files))
}

## ----load package narrowPeak-----------------------------------------------
# load peak calls
np_grs = easyLoad_narrowPeak(np_files)

## ----overlap peaks---------------------------------------------------------
olaps = ssvOverlapIntervalSets(np_grs)

## ----ctcf fig1,hold=TRUE, fig.align='center', fig.height=4, fig.width = 8----
p_bars = ssvFeatureBars(olaps, show_counts = FALSE) +  
    theme(legend.position = "left")
p_bars = p_bars + theme(axis.text.x = element_blank(), 
                axis.ticks.x = element_blank(),
                legend.justification = "center") +
    labs(fill = "cell line")
p_venn = ssvFeatureVenn(olaps, counts_txt_size = 4) + 
    guides(fill = "none", color = "none")
p_euler = ssvFeatureEuler(olaps) + 
    guides(fill = "none", color = "none")

cowplot::ggdraw() +
    cowplot::draw_plot(p_bars, x = 0, y = 0, width = .4, height = .7) +
    cowplot::draw_plot(p_venn, x = .4, y = .1, width = .3, height = .7) +
    cowplot::draw_plot(p_euler, x = 0.7, y = .1, width = 0.3, height = .7) +
    cowplot::draw_plot_label(c("CTCF binding in breast cancer cell lines", 
                               "A", "B", "C"),
                             x = c(.04, .17, 0.4, .7),
                             y = c(.92, .75, .75, .75), size = 10, hjust = 0)

## ----color change, eval=FALSE,hold = TRUE, fig.align='center', fig.height=4, fig.width = 8----
#  col_vals = c("MCF10A_CTCF" = 'red',
#               "MCF10AT1_CTCF" = "blue",
#               "MCF10CA1_CTCF" = "green")
#  sf = scale_fill_manual(values = col_vals)
#  sc = scale_color_manual(values = col_vals)
#  cowplot::ggdraw() +
#      cowplot::draw_plot(p_bars + sf,
#                         x = 0, y = 0,
#                         width = .4, height = .7) +
#      cowplot::draw_plot(p_venn + sf + sc,
#                         x = .4, y = .1,
#                         width = .3, height = .7) +
#      cowplot::draw_plot(p_euler + sf + sc,
#                         x = 0.7, y = .1,
#                         width = 0.3, height = .7) +
#      cowplot::draw_plot_label(c("CTCF binding in breast cancer cell lines",
#                                 "A", "B", "C"),
#                               x = c(.04, .17, 0.4, .7),
#                               y = c(.92, .75, .75, .75), size = 10, hjust = 0)

## ----load ChIPpeakAnno, message=FALSE--------------------------------------
library(ChIPpeakAnno)
data(TSS.human.GRCh38)
macs.anno <- annotatePeakInBatch(olaps, AnnotationData=TSS.human.GRCh38)

## ----distance filter-------------------------------------------------------
macs.anno = subset(macs.anno, distancetoFeature < 1000)

## ----subset 1--------------------------------------------------------------
subset(macs.anno, MCF10AT1_CTCF & MCF10A_CTCF & !MCF10CA1_CTCF)$feature

## ----subset 2--------------------------------------------------------------
subset(macs.anno, MCF10A_CTCF & !MCF10AT1_CTCF & !MCF10CA1_CTCF)$feature

## ----set fixed width, fig.height=3, fig.width=3----------------------------
window_size = 50
width_q75 = quantile(width(olaps), .75)
width_q75 = ceiling(width_q75 / window_size) * window_size
hist_res = hist(width(olaps))
lines(rep(width_q75, 2), c(0, max(hist_res$counts)), col = "red", lwd = 5)
text(width_q75, max(hist_res$counts), "fixedSize", adj = c(-.1, 1), col = "red")
## apply width
olaps_fixedSize = centerFixedSizeGRanges(olaps, width_q75)

## ----fetch package bw------------------------------------------------------
if(use_full_data){
    bw_gr = fetchWindowedBigwigList(file_paths = bw_files, 
                                    qgr = olaps_fixedSize,
                                    win_size = 50)
}

## ----ctcf scatterplots, fig.width=10, fig.height=4, message=FALSE----------
# shortening colnames will make group names less cumbersome in plot legend
colnames(mcols(olaps_fixedSize)) = sub("_CTCF", "", 
                                       colnames(mcols(olaps_fixedSize)))

all_groups = levels(ssvFactorizeMembTable(
    ssvMakeMembTable(olaps_fixedSize))$group)
all_colors = RColorBrewer::brewer.pal(length(all_groups), "Set1")
all_colors[5:7] = safeBrew(3, "Dark2")
names(all_colors) = all_groups
olap_groups_12 = ssvFactorizeMembTable(
    ssvMakeMembTable(olaps_fixedSize)[, 1:2])
p_12 = ssvSignalScatterplot(bw_gr, 
                            x_name = "MCF10A_CTCF", 
                            y_name = "MCF10AT1_CTCF", 
                            color_table = olap_groups_12) + 
    scale_color_manual(values = all_colors)

olap_groups_13 = ssvFactorizeMembTable(
    ssvMakeMembTable(olaps_fixedSize)[, c(1,3)])
p_13  = ssvSignalScatterplot(bw_gr, 
                             x_name = "MCF10A_CTCF", 
                             y_name = "MCF10CA1_CTCF", 
                     color_table = olap_groups_13) + 
    scale_color_manual(values = all_colors)


if(use_full_data){
    tp_12 = p_12 + scale_size_continuous(range = .1) + 
        scale_alpha_continuous(range = .1) + 
        geom_density2d(aes(color = group), h = 40, bins = 3)
    tp_13 = p_13 + scale_size_continuous(range = .1) + 
        scale_alpha_continuous(range = .1) + 
        geom_density2d(aes(color = group), h = 40, bins = 3)
    cowplot::plot_grid(tp_12 + labs(title = ""), 
                       tp_13 + labs(title = ""), 
                       label_y = .85, labels = "AUTO") 
}else{
    cowplot::plot_grid(p_12 + labs(title = ""), 
                       p_13 + labs(title = ""), 
                       label_y = .85, labels = "AUTO")
}


## ----ctcf heatmap, message=FALSE, fig.width=5------------------------------
bw_gr$facet_label = sub("_", "\n", bw_gr$sample)
clust_gr = ssvSignalClustering(bw_gr, nclust = 3, facet_ = "facet_label")
ssvSignalHeatmap(clust_gr, facet_ = "facet_label") + labs(fill = "FE", 
                                  y = "region", 
                                  x = "bp from center")

## ----ctcf recentered heatmap, message = FALSE, fig.width=10, fig.height=6----
center_gr = centerAtMax(clust_gr, view_size = 150, 
                        by_ = "id", check_by_dupes = FALSE)


p_center_hmap = ssvSignalHeatmap(center_gr, facet_ = "facet_label") + 
    labs(fill = "FE",
         y = "region",
         x = "bp from center")

## since center_gr still retains clustering information, clustering is not
## repeated by default, the following reclusters the data.
clust_center_gr = ssvSignalClustering(center_gr, nclust = 3)
p_center_hmap_reclust = ssvSignalHeatmap(clust_center_gr, 
                                         facet_ = "facet_label") + 
    labs(fill = "FE",
         y = "region",
         x = "bp from center")
cowplot::plot_grid(p_center_hmap + labs(title = "original clustering"), 
                   p_center_hmap_reclust + labs(title = "reclustered"))

## ----cluster annotation----------------------------------------------------
clust_df = as.data.frame(mcols(clust_gr))
clust_df = unique(clust_df[,c("id", "cluster_id")])

olap_clust_annot = olap_gr
mcols(olap_clust_annot) = data.frame(id = seq_along(olap_clust_annot))
olap_clust_annot = GRanges(merge(olap_clust_annot, clust_df))

olap_clust_annot = subset(olap_clust_annot, cluster_id %in% 1:2)

olap_clust_annot <- annotatePeakInBatch(olap_clust_annot, 
                                        AnnotationData=TSS.human.GRCh38)
olap_clust_annot$feature

## --------------------------------------------------------------------------
target_data = "MCF10A_CTCF"
chmm_win = 200 #window size is an important chromHMM parameter.
# 200 is the default window size and matches the state segementation

if(use_full_data){
    # set all file paths
    chmm_bw_file = bfcrpath(ssv_bfc, rnames = paste(target_data, "bigwig", 
                                                    sep = ","))
    chmm_np_file = bfcrpath(ssv_bfc, rnames = paste(target_data, "narrowPeak", 
                                                    sep = ","))
    chmm_seg_file = bfcrpath(ssv_bfc, rnames = "MCF7,segmentation",
                             fpath = chromHMM_demo_segmentation_url)
    
    query_chain = bfcquery(ssv_bfc, "hg19ToHg38,chain")
    if(nrow(query_chain) == 0){
        chain_hg19ToHg38_gz = bfcrpath(ssv_bfc, rnames = "hg19ToHg38,gz",
                                       fpath = chromHMM_demo_chain_url)
        ch_raw = readLines(gzfile(chain_hg19ToHg38_gz))
        ch_file = bfcnew(ssv_bfc, rname = "hg19ToHg38,chain")
        writeLines(ch_raw, con = ch_file)
        
    }
    chmm_chain_file = bfcrpath(ssv_bfc, rnames = "hg19ToHg38,chain")
    ch = rtracklayer::import.chain(chmm_chain_file)    
    # load segmentation data
    chmm_gr = rtracklayer::import.bed(chmm_seg_file)
    #cleanup state names.
    chmm_gr$name = gsub("\\+", " and ", chmm_gr$name)
    chmm_gr$name = gsub("_", " ", chmm_gr$name)
    #setup color to state mapping
    colDF = unique(mcols(chmm_gr)[c("name", "itemRgb")])
    state_colors = colDF$itemRgb
    names(state_colors) = colDF$name
    #liftover states from hg19 to hg38
    ch = rtracklayer::import.chain(chmm_chain_file)
    chmm_gr_hg38 = rtracklayer::liftOver(chmm_gr, ch)
    chmm_gr_hg38 = unlist(chmm_gr_hg38)
    chmm_grs_list = as.list(GenomicRanges::split(chmm_gr_hg38, 
                                                 chmm_gr_hg38$name))
    #transform narrowPeak ranges to summit positions
    chmm_np_grs = easyLoad_narrowPeak(chmm_np_file, file_names = target_data)
    chmm_summit_grs = lapply(chmm_np_grs, function(x){
        start(x) = start(x) + x$relSummit
        end(x) = start(x)
        x
    })
    qlist = append(chmm_summit_grs[1], chmm_grs_list)
    chmm_olaps = ssvOverlapIntervalSets(qlist, use_first = TRUE)
    #discard the columns for peak call and no_hit, not informative here.
    mcols(chmm_olaps)[[1]] = NULL
    chmm_olaps$no_hit = NULL
    #total width of genome assigned each state
    state_total_widths = sapply(chmm_grs_list, function(my_gr){
        sum(as.numeric(width(my_gr)))
    })
    #Expand state regions into 200 bp windows.
    state_wingrs = lapply(chmm_grs_list, function(my_gr){
        st = my_gr$name[1]
        wgr = unlist(slidingWindows(my_gr, chmm_win, chmm_win))
        wgr$state = st
        wgr
    })
    state_wingrs = unlist(GRangesList(state_wingrs))
    # fetch bw data for each state
    # it probably isn't useful to grab every single window for each state
    # so we can cap the number of each state carried forward
    max_per_state = 5000
    # flank size zooms out a bit from each chromHMM window
    flank_size = 400
    
    state_split = split(state_wingrs, state_wingrs$state)
    state_split = lapply(state_split, function(my_gr){
        samp_gr = sample(my_gr, min(length(my_gr), max_per_state))
        samp_gr = sort(samp_gr)
        names(samp_gr) = seq_along(samp_gr)
        samp_gr
    })
    state_gr = unlist(GRangesList(state_split))
    state_gr = resize(state_gr, width = chmm_win + 2 * flank_size, 
                      fix = "center")
    
    bw_states_gr = fetchWindowedBigwig(bw_file = chmm_bw_file, 
                                       qgr = state_gr, 
                                       win_size = 50)
    bw_states_gr$grp = sub("\\..+", "", bw_states_gr$id)
    bw_states_gr$grp_id = sub(".+\\.", "", bw_states_gr$id)
}else{
    max_per_state = 20
    flank_size = 400
    state_colors = chromHMM_demo_state_colors
    bw_states_gr = chromHMM_demo_bw_states_gr
    chmm_olaps = chromHMM_demo_overlaps_gr
    state_total_widths = chromHMM_demo_state_total_widths
}


## ----state raw, message = FALSE, fig.width=3, fig.height=3-----------------
olaps_df = as.data.frame(mcols(chmm_olaps))
colnames(olaps_df) = gsub("\\.", " ", colnames(olaps_df))

p_state_raw_count = ssvFeatureBars(olaps_df, show_counts = FALSE) +
    labs(fill = "state", x = "") +
    scale_fill_manual(values = state_colors) +
    theme_cowplot() + guides(fill = "none") +
    theme(axis.text.x = element_text(angle = 90, size = 8,
                                     hjust = 1, vjust = .5))
p_state_raw_count

## ----state enrichment, fig.width=3, fig.height=3---------------------------
state_width_fraction = state_total_widths / sum(state_total_widths)

state_peak_hits =  colSums(olaps_df)
state_peak_fraction = state_peak_hits / sum(state_peak_hits)

enrichment = state_peak_fraction / 
    state_width_fraction[names(state_peak_fraction)]
enrich_df = data.frame(state = names(enrichment), enrichment = enrichment)
p_state_enrichment =  ggplot(enrich_df) +
    geom_bar(aes(x = state, fill = state, y = enrichment), stat = "identity") +
    labs(x = "") + 
    theme_cowplot() + guides(fill = "none") +
    scale_fill_manual(values = state_colors) +
    theme(axis.text.x = element_text(angle = 90, size = 8,
                                     hjust = 1, vjust = .5))
p_state_enrichment

## ---- message=FALSE, fig.width=6-------------------------------------------
p_agg_tracks = ssvSignalLineplotAgg(bw_states_gr, 
                                    sample_ = "grp", 
                                    color_ = "grp")
gb = ggplot2::ggplot_build(p_agg_tracks)
yrng = range(gb$data[[1]]$y)
p_agg_tracks = p_agg_tracks +
    scale_color_manual(values = state_colors) +
    annotate("line", x = rep(-chmm_win/2, 2), y = yrng) +
    annotate("line", x = rep(chmm_win/2, 2), y = yrng) +
    labs(y = "FE", x = "bp", color = "state", 
         title = paste("Average FE by state,", target_data), 
         subtitle = paste("states sampled to", 
                          max_per_state, 
                          chmm_win, 
                          "bp windows each\n", 
                          flank_size, 
                          "bp flanking each side")) +
    theme(plot.title = element_text(hjust = 0))
p_agg_tracks

## ----state heatmap, fig.width=8--------------------------------------------
pdt = as.data.table(mcols(bw_states_gr))
pdt$grp_id = as.integer(pdt$grp_id)
# reassign grp_id to sort within each state set
dt_list = lapply(unique(pdt$grp), function(state){
    dt = pdt[grp == state]
    dtmax = dt[, .(ymax = y[which(x == x[order(abs(x))][1])]), by = grp_id]
    dtmax = dtmax[order(ymax, decreasing = TRUE)]
    dtmax[, grp_o := seq_len(.N)]
    dtmax$ymax = NULL
    dt = merge(dtmax, dt)
    dt[, grp_id := grp_o]
    dt$grp_o = NULL
    dt
})
# reassemble
pdt = rbindlist(dt_list)
# heatmap facetted by state and sorted in decreasing order
p_state_hmap = ggplot(pdt) + 
    geom_raster(aes(x = x, y = grp_id, fill = y)) + 
    scale_y_reverse() +
    facet_wrap("grp", nrow = 2) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          strip.text = element_text(size= 8)) +
    scale_fill_gradientn(colors = c("white", "orange", "red")) +
    labs(y = "", fill = "FE", 
         x = "bp from local summit", 
         title = paste(max_per_state, "random regions per state"))
p_state_hmap

## ---- fig.height=12, fig.width=8-------------------------------------------
bar_height = .25
line_height = .3
ggdraw() +
    draw_plot(p_state_raw_count + guides(fill = "none"), 
              x = 0, 
              y = 1 - bar_height, 
              width = .5, 
              height = bar_height) +
    draw_plot(p_state_enrichment + guides(fill = "none"), 
              x = .5, 
              y = 1 - bar_height, 
              width = .5, 
              height = bar_height) +
    draw_plot(p_agg_tracks, 
              x = 0, 
              y = 1 - bar_height - line_height, 
              width = 1, 
              height = line_height) +
    draw_plot(p_state_hmap, 
              x = 0, 
              y = 0, 
              width = 1, 
              height = 1 - bar_height - line_height) +
    draw_plot_label(LETTERS[1:4], 
                    c(0, 0.48, 0, 0), 
                    c(1, 1, 1 - bar_height, 1 - bar_height - line_height), 
                    size = 15)

