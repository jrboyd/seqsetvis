testthat::context("cluster sorting")

library(seqsetvis)
library(testthat)
library(GenomicRanges)

test_object = CTCF_in_10a_profiles_dt
set.seed(0)
clust_dt = ssvSignalClustering(test_object, nclust = 3)
show_plots = FALSE

test_that("within_clust_sort", {
    resort_dt = within_clust_sort(clust_dt, within_order_strategy = "hclust")
    resort_dt2 = within_clust_sort(clust_dt, within_order_strategy = "left")
    resort_dt3 = within_clust_sort(clust_dt, within_order_strategy = "right")

    if(show_plots){
        cowplot::plot_grid(ssvSignalHeatmap(clust_dt) + labs(title = "sort"),
                           ssvSignalHeatmap(resort_dt) + labs(title = "hclust"),
                           ssvSignalHeatmap(resort_dt2) + labs(title = "left"),
                           ssvSignalHeatmap(resort_dt3) + labs(title = "right"))
    }


    #test id levels still the same set
    expect_setequal(levels(clust_dt$id), levels(resort_dt$id))
    expect_setequal(levels(clust_dt$id), levels(resort_dt2$id))
    expect_setequal(levels(clust_dt$id), levels(resort_dt3$id))

    #test id levels are in different order
    expect_failure(expect_equal(levels(clust_dt$id), levels(resort_dt$id)))
    expect_failure(expect_equal(levels(clust_dt$id), levels(resort_dt2$id)))
    expect_failure(expect_equal(levels(clust_dt$id), levels(resort_dt3$id)))

    all_assign = lapply(list(clust_dt, resort_dt, resort_dt2, resort_dt3), function(x){
        unique(x[, .(id, cluster_id)])
    })
    all_membs = lapply(all_assign, function(x){
        split(as.character(x$id), x$cluster_id)
    })

    #test clusters haven't changed
    for(i in seq_along(all_membs[[1]])){
        expect_equal(all_membs[[1]][[i]], all_membs[[2]][[i]])
        expect_equal(all_membs[[1]][[i]], all_membs[[3]][[i]])
        expect_equal(all_membs[[1]][[i]], all_membs[[4]][[i]])
    }
})

test_that("reorder_clusters_manual", {
    reo_dt = reorder_clusters_manual(clust_dt, manual_order = c(2, 3))
    reo_dt.no_reapply = reorder_clusters_manual(clust_dt, manual_order = c(2, 3), reapply_cluster_names = FALSE)

    if(show_plots){
        cowplot::plot_grid(ssvSignalHeatmap(clust_dt) + labs(title = "original"),
                           ssvSignalHeatmap(reo_dt) + labs(title = "reordered"),
                           ssvSignalHeatmap(reo_dt.no_reapply) + labs(title = "no reapply"))
    }

    #test id levels still the same set
    expect_setequal(levels(clust_dt$id), levels(reo_dt$id))
    expect_setequal(levels(clust_dt$id), levels(reo_dt.no_reapply$id))

    #test id levels are in different order
    expect_failure(expect_equal(levels(clust_dt$id), levels(reo_dt$id)))
    expect_failure(expect_equal(levels(clust_dt$id), levels(reo_dt.no_reapply$id)))

    all_assign = lapply(list(clust_dt[order(id)], reo_dt[order(id)], reo_dt.no_reapply[order(id)]), function(x){
        unique(x[, .(id, cluster_id)])
    })
    all_membs = lapply(all_assign, function(x){
        split(as.character(x$id), x$cluster_id)
    })

    #test clusters haven't changed
    for(i in as.character(seq_along(all_membs[[1]]))){
        #this cluster should change
        expect_failure(expect_equal(all_membs[[1]][[i]], all_membs[[2]][[i]]))
        #this cluster should not
        expect_equal(all_membs[[1]][[i]], all_membs[[3]][[i]])
        # expect_setequal(all_membs[[1]][[i]], all_membs[[3]][[i]])
    }
})

test_that("reorder_clusters_hclust", {
    set.seed(0)
    clust_dt10.raw = ssvSignalClustering(test_object, nclust = 10)
    clust_dt10 = reorder_clusters_manual(clust_dt10.raw, c(3, 10), reapply_cluster_names = FALSE)
    clust_dt10.hclust = reorder_clusters_hclust(clust_dt10)
    clust_dt10.hclust.no_reapply = reorder_clusters_hclust(clust_dt10, reapply_cluster_names = FALSE)

    if(show_plots){
        cowplot::plot_grid(nrow = 1,
                           ssvSignalHeatmap(clust_dt10.raw) + labs(title = "raw"),
                           ssvSignalHeatmap(clust_dt10) + labs(title = "original"),
                           ssvSignalHeatmap(clust_dt10.hclust) + labs(title = "reordered"),
                           ssvSignalHeatmap(clust_dt10.hclust.no_reapply) + labs(title = "no reapply"))
    }

    #test id levels still the same set
    expect_setequal(levels(clust_dt10$id), levels(clust_dt10.hclust$id))
    expect_setequal(levels(clust_dt10$id), levels(clust_dt10.hclust.no_reapply$id))

    #test id levels are in different order
    expect_failure(expect_equal(levels(clust_dt10$id), levels(clust_dt10.hclust$id)))
    expect_failure(expect_equal(levels(clust_dt10$id), levels(clust_dt10.hclust.no_reapply$id)))

    all_assign = lapply(list(clust_dt10[order(id)], clust_dt10.hclust[order(id)], clust_dt10.hclust.no_reapply[order(id)]), function(x){
        unique(x[, .(id, cluster_id)])
    })
    all_membs = lapply(all_assign, function(x){
        split(as.character(x$id), x$cluster_id)
    })

    num_fail_1 = 0
    num_fail_2 = 0
    #test clusters haven't changed
    for(i in as.character(seq_along(all_membs[[1]]))){
        #this cluster should change
        if(length(all_membs[[1]][[i]]) != length(all_membs[[2]][[i]])){
            num_fail_1 = num_fail_1 + 1
        }else if(all(all_membs[[1]][[i]] != all_membs[[2]][[i]])){
            num_fail_1 = num_fail_1 + 1
        }
        if(length(all_membs[[1]][[i]]) != length(all_membs[[3]][[i]])){
            num_fail_2 = num_fail_2 + 1
        }else if(all(all_membs[[1]][[i]] != all_membs[[3]][[i]])){
            num_fail_2 = num_fail_2 + 1
        }

        # expect_equal()
        # #this cluster should not
        # expect_failure(expect_equal(all_membs[[1]][[i]], all_membs[[3]][[i]]))
        # # expect_setequal(all_membs[[1]][[i]], all_membs[[3]][[i]])
    }
    expect_gt(num_fail_1, 0)
    expect_equal(num_fail_2, 0)
})
