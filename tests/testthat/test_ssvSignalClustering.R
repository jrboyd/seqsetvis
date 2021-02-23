testthat::context("SignalClustering")

#new feature - ssvSignalClustering accepts memb_table
#new feature - ssvSignalHeatmap.ClusterBars displays cluster bars once on the left instead of in each facet
#new feature - ssvSignalClustering accepts centroids instead of nclust

library(seqsetvis)
library(testthat)
library(GenomicRanges)

test_that("ssvSignalHeatmap.ClusterBars", {
    p_heat = ssvSignalHeatmap.ClusterBars(CTCF_in_10a_profiles_dt)
    class(p_heat)
    expect_s3_class(p_heat, "ggplot")

    p_list = ssvSignalHeatmap.ClusterBars(CTCF_in_10a_profiles_dt, return_unassembled_plots = TRUE)
    class(p_list)
    class(p_list[[1]])
    expect_is(p_list, "list")
    expect_s3_class(p_list[[1]], "ggplot")
    expect_s3_class(p_list[[2]], "ggplot")
})

test_that("ssvSignalHeatmap fill_limits", {
    p_heat = ssvSignalHeatmap(CTCF_in_10a_profiles_dt, fill_limits = c(20, 60))
    expect_equal(max(p_heat$data$y), 60)
    expect_equal(min(p_heat$data$y), 20)

    p_heat2 = ssvSignalHeatmap.ClusterBars(CTCF_in_10a_profiles_dt, fill_limits = c(20, 60), return_unassembled_plots = TRUE)[[2]]
    expect_equal(max(p_heat2$data$y), 60)
    expect_equal(min(p_heat2$data$y), 20)
})

test_that("ssvSignalHeatmap memb_table", {
    memb_table = ssvMakeMembTable(CTCF_in_10a_narrowPeak_grs)
    setequal(rownames(memb_table), CTCF_in_10a_profiles_dt$id)
    clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt, memb_table = memb_table)
    p_heat = ssvSignalHeatmap(clust_dt, show_cluster_bars = FALSE)
    expect_s3_class(p_heat, "ggplot")

    p_heat.anno = add_cluster_annotation(p = p_heat, cluster_ids = clust_dt, xleft = -400, xright = -360, show_labels = FALSE, rect_colors = safeBrew(3))
    expect_s3_class(p_heat.anno, "ggplot")
})

test_that("ssvSignalHeatmap centroids from tidy", {
    centroids.tidy = CTCF_in_10a_profiles_dt[id %in% 1:3]
    data.table::setnames(centroids.tidy, "id", "cluster_id")

    cent_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt, k_centroids = centroids.tidy)
    p_heat = ssvSignalHeatmap(cent_dt, show_cluster_bars = TRUE)
    expect_s3_class(p_heat, "ggplot")
})

test_that("ssvSignalHeatmap centroids from wide", {
    centroids.tidy = CTCF_in_10a_profiles_dt[id %in% 1:3]
    data.table::setnames(centroids.tidy, "id", "cluster_id")

    wide_mat = make_clustering_matrix(centroids.tidy, row_ = "cluster_id")

    cent_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt, k_centroids = wide_mat)
    p_heat = ssvSignalHeatmap(cent_dt, show_cluster_bars = TRUE)
    expect_s3_class(p_heat, "ggplot")
})
