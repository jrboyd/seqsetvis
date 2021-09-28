testthat::context("with_sort_clust")

library(seqsetvis)
library(testthat)
library(GenomicRanges)

test_object = CTCF_in_10a_profiles_dt

clust_dt = ssvSignalClustering(test_object, nclust = 3)

resort_dt = within_clust_sort(clust_dt, within_order_strategy = "hclust")

cowplot::plot_grid(ssvSignalHeatmap(clust_dt), ssvSignalHeatmap(resort_dt))

