
valid_sort_strategies = c("hclust", "sort", "left", "right", "reverse")

#' Clustering as for a heatmap.  This is used internally by
#' \code{\link{ssvSignalHeatmap}} but can also be run before calling
#' ssvSignalHeatmap for greater control and access to clustering results
#' directly.
#'
#' Clustering is via k-means by default.  The number of clusters is determined
#' by nclust. Optionally, k-means can be initialized with a data.frame provided
#' to k_centroids. As an alternative to k-means, a membership table from
#' \code{\link{ssvMakeMembTable}} can be provided to determine logical clusters.
#'
#' Within each cluster, items will either be sorted by decreasing average signal
#' or hierachically clustered; this is controlled via within_order_strategy.
#'
#' @export
#'
#' @param bw_data a GRanges or data.table of bigwig signal. As returned from
#'   \code{\link{ssvFetchBam}} and \code{\link{ssvFetchBigwig}}
#' @param nclust Number of clusters.  Defaults to 6 if nclust, k_centroids, and
#'   memb_table are not set.
#' @param k_centroids data.frame of centroids for k-means clusters. Incompatible
#'   with nclust or memb_table.
#' @param memb_table Membership table as from \code{\link{ssvMakeMembTable}}.
#'   Logical groups from membership table will be clusters. Incompatible with
#'   nclust or k_centroids.
#' @param row_ variable name mapped to row, likely id or gene name for ngs data.
#'   Default is "id" and works with ssvFetch* output.
#' @param column_ varaible mapped to column, likely bp position for ngs data.
#'   Default is "x" and works with ssvFetch* output.
#' @param fill_ numeric variable to map to fill. Default is "y" and works with
#'   ssvFetch* output.
#' @param facet_ variable name to facet horizontally by. Default is "sample" and
#'   works with ssvFetch* output. Set to "" if data is not facetted.
#' @param cluster_ variable name to use for cluster info. Default is
#'   "cluster_id".
#' @param max_rows for speed rows are sampled to 500 by default, use Inf to plot
#'   full data
#' @param max_cols for speed columns are sampled to 100 by default, use Inf to
#'   plot full data
#' @param clustering_col_min numeric minimum for col range considered when
#'   clustering, default in -Inf
#' @param clustering_col_max numeric maximum for col range considered when
#'   clustering, default in Inf
#' @param within_order_strategy one of "hclust", "sort", "right", "left",
#'   "reverse".  If "hclust", hierarchical clustering will be used. If "sort", a
#'   simple decreasing sort of rosSums.  If "left", will atttempt to put high
#'   signal on left ("right" is opposite).  If "reverse" reverses existing order
#'   (should only be used after meaningful order imposed).
#' @param dcast_fill value to supply to dcast fill argument. default is NA.
#' @param iter.max Number of max iterations to allow for k-means. Default is 30.
#' @param fun.aggregate Function to aggregate when multiple values present for
#'   facet_, row_, and column_. The
#'   function should accept a single vector argument or be a character string
#'   naming such a function.
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @return data.table of signal profiles, ready for ssvSignalHeatmap
#' @examples
#' clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_gr)
#' ssvSignalHeatmap(clust_dt)
#'
#' clust_dt2 = ssvSignalClustering(CTCF_in_10a_profiles_gr, nclust = 2)
#' ssvSignalHeatmap(clust_dt2)
#'
#' #clustering can be targetted to a specific part of the region
#' clust_dt3 = ssvSignalClustering(CTCF_in_10a_profiles_gr, nclust = 2,
#'     clustering_col_min = -250, clustering_col_max = -150)
#' ssvSignalHeatmap(clust_dt3)
#' clust_dt4 = ssvSignalClustering(CTCF_in_10a_profiles_gr, nclust = 2,
#'     clustering_col_min = 150, clustering_col_max = 250)
#' ssvSignalHeatmap(clust_dt4)
ssvSignalClustering = function(bw_data,
                               nclust = NULL,
                               k_centroids = NULL,
                               memb_table = NULL,
                               row_ = "id",
                               column_ = "x",
                               fill_ = "y",
                               facet_ = "sample",
                               cluster_ = "cluster_id",
                               max_rows = 500,
                               max_cols = 100,
                               clustering_col_min = -Inf,
                               clustering_col_max = Inf,
                               within_order_strategy = valid_sort_strategies[2],
                               dcast_fill = NA,
                               iter.max = 30,
                               fun.aggregate = "mean"){
    message("clustering...")
    id__ = xbp = x = to_disp = y = hit = val = y = y_gap = group__ =  NULL#declare binding for data.table
    output_GRanges = FALSE
    if(is(bw_data, "GRanges")){
        bw_data = data.table::as.data.table(bw_data)
        output_GRanges = TRUE
    }
    stopifnot(is.data.table(bw_data))
    if(!is.null(k_centroids) & !is.null(memb_table)){
        stop("only one of k_centroids or memb_table is allowed.")
    }
    if(!is.null(k_centroids)){
        nclust = nrow(k_centroids)
    }
    if(!is.null(memb_table)){ #memb_table get handled to provide row_ to cluster_ mapping
        if(is.data.frame(memb_table)){
            if(!is.null(memb_table[[row_]]) & !is.null(memb_table[[cluster_]])){
                #in this case memb_table is already a valid row_ to cluster_ mapping
                #nothing needs done
            }else{
                if(is.null(memb_table[[row_]])){
                    memb_table[[row_]] = memb_table[["id"]]
                    memb_table[["id"]] = NULL
                }
                if(is.null(memb_table[[cluster_]])){
                    memb_table[[cluster_]] = memb_table[["group"]]
                    memb_table[["group"]] = NULL
                }
            }
        }else{
            memb_table = ssvFactorizeMembTable(memb_table)
            if(is.null(memb_table[[row_]])){
                memb_table[[row_]] = memb_table[["id"]]
                memb_table[["id"]] = NULL
            }
            if(is.null(memb_table[[cluster_]])){
                memb_table[[cluster_]] = memb_table[["group"]]
                memb_table[["group"]] = NULL
            }
        }
        if(is.null(memb_table[[row_]])) stop("Could not correctly handle memb_table. It should look like output of ssvMakeMembTable or ssvFactorizeMembTable.")
        if(is.null(memb_table[[cluster_]])) stop("Could not correctly handle memb_table. It should look like output of ssvMakeMembTable or ssvFactorizeMembTable.")
        nclust = length(unique(memb_table[[cluster_]]))
        cluster_assignment = memb_table[[cluster_]]
        names(cluster_assignment) = memb_table[[row_]]
    }else{
        cluster_assignment = NULL
    }
    if(is.null(nclust) & is.null(k_centroids) & is.null(memb_table)){
        nclust = 6
        # stop("one of nclust, k_centroids, or memb_table must be set.")
    }
    stopifnot(is.numeric(nclust))
    stopifnot(is.character(row_), is.character(column_), is.character(fill_),
              is.character(facet_), is.character(cluster_))
    stopifnot(row_ %in% colnames(bw_data), column_ %in% colnames(bw_data),
              fill_ %in% colnames(bw_data))
    stopifnot(facet_ %in% colnames(bw_data) || facet_ == "")
    stopifnot(is.numeric(max_rows), is.numeric(max_cols),
              is.numeric(clustering_col_min), is.numeric(clustering_col_max))

    plot_dt = data.table::copy(bw_data)
    dc_mat = make_clustering_matrix(
        plot_dt,
        row_ = row_,
        column_ = column_,
        fill_ = fill_,
        facet_ = facet_,
        max_rows = max_rows,
        max_cols = max_cols,
        clustering_col_min = clustering_col_min,
        clustering_col_max = clustering_col_max,
        dcast_fill = dcast_fill,
        fun.aggregate = fun.aggregate

    )

    if(!is.null(k_centroids)){
        if(!is.matrix(k_centroids)){
            if(!is.data.table(k_centroids)) stop("k_centroids must be a wide matrix of centroids or a tidy data.table")
            if(is.null(k_centroids[[cluster_]])) stop("k_centroids looks like a tidy data.table but is missing the cluster variable ", cluster_)
            k_centroids = make_clustering_matrix(
                k_centroids,
                row_ = cluster_,
                column_ = column_,
                fill_ = fill_, facet_ = facet_,
                max_rows = Inf,
                max_cols = max_cols,
                clustering_col_min = clustering_col_min,
                clustering_col_max = clustering_col_max,
                dcast_fill = dcast_fill
            )
        }
    }


    rclusters = clusteringKmeansNestedHclust(dc_mat,
                                             nclust = nclust,
                                             within_order_strategy,
                                             centroids = k_centroids,
                                             manual_mapping = cluster_assignment,
                                             iter.max = iter.max)
    rclusters = rclusters[rev(seq_len(nrow(rclusters))),]
    plot_dt = plot_dt[get(row_) %in% rclusters[["id__"]]]
    plot_dt[[row_]] = factor(plot_dt[[row_]], levels = rev(rclusters[["id__"]]))
    data.table::setkey(rclusters, id__)
    plot_dt[[cluster_]] = rclusters[list(plot_dt[[row_]]), group__]
    agg_cn = c(row_, column_, facet_)
    agg_cn = agg_cn[agg_cn != ""]
    agg_required = max(plot_dt[, .N, c(agg_cn)]$N) > 1
    if(agg_required){
        #apply fun.aggregate to plotted data
        constant_cn = c("seqnames", "start", "end", 'width', "strand", "cluster_id")
        constant_cn = constant_cn[constant_cn %in% colnames(plot_dt)]
        if(is.character(fun.aggregate)){
            fun.aggregate = get(fun.aggregate)
        }
        plot_dt.agg = plot_dt[, fun.aggregate(get(fill_)), c(unique(c(agg_cn, constant_cn)))]
        setnames(plot_dt.agg, "V1", fill_)
        message("Aggregation has occured (multiple samples per facet), some variables have been dropped.")
    }else{
        plot_dt.agg = plot_dt
    }
    if(output_GRanges){
        plot_dt.agg = GRanges(plot_dt.agg)
    }
    return(plot_dt.agg)
}

#' perform kmeans clustering on matrix rows and return reordered matrix along
#' with order matched cluster assignments. clusters are sorted using hclust on
#' centers
#'
#' @param mat numeric matrix to cluster.
#' @param nclust the number of clusters.
#' @param centroids optional matrix with same columns as mat and one centroid
#'   per row to base clusters off of.  Overrides any setting to nclust. Default
#'   of NULL results in randomly initialized k-means.
#' @param iter.max Number of max iterations to allow for k-means. Default is 30.
#'
#' @return data.table with group__ variable indicating cluster membership and id__
#'   variable that is a factor indicating order based on within cluster
#'   similarity
#' @export
#' @importFrom stats kmeans hclust dist
#' @examples
#' dt = data.table::copy(CTCF_in_10a_profiles_dt)
#' mat = data.table::dcast(dt, id ~ sample + x, value.var = "y" )
#' rn = mat$id
#' mat = as.matrix(mat[,-1])
#' rownames(mat) = rn
#' clust_dt = clusteringKmeans(mat, nclust = 3)
#' dt = merge(dt, clust_dt[, .(id = id__, group = group__)])
#' dt$id = factor(dt$id, levels = clust_dt$id)
#' dt[order(id)]
clusteringKmeans = function(mat, nclust, centroids = NULL, iter.max = 30) {
    if(!is.null(centroids)){
        if(!all(colnames(mat) == colnames(centroids))){
            stop("provided centroids matrix did not have identical column names to mat matrix.")
        }
        nclust = nrow(centroids)
        centers_arg = centroids
    }else{
        centers_arg = nclust
    }
    stopifnot(is.numeric(nclust))
    cluster_ordered = mat_name = NULL#declare binding for data.table
    if(nrow(mat) <= nclust){
        nclust = nrow(mat)
        mat_kmclust = list(
            centers = mat,
            cluster = seq(nrow(mat))
        )
        names(mat_kmclust$cluster) = seq(nrow(mat))

    }else{
        if(is.null(centroids)){
            if(nrow(unique(mat)) < nclust){
                nclust = nrow(unique(mat))
                warning("Reducing nclust to ", nclust,
                        " - maximum number of clusters allowed due to ",
                        "low uniqueness.")
                centers_arg = nclust
            }

        }
        mat_kmclust = stats::kmeans(mat, centers = centers_arg, iter.max = iter.max)
    }

    if(nclust < 3 || !is.null(centroids)){
        center_o = seq(nclust)
    }else{
        center_o = stats::hclust(stats::dist(mat_kmclust$centers))$order
    }
    center_reo = seq_along(center_o)
    names(center_reo) = center_o
    center_reo[as.character(mat_kmclust$cluster)]
    mat_dt = data.table::data.table(
        mat_name = names(mat_kmclust$cluster),
        cluster = mat_kmclust$cluster,
        cluster_ordered = center_reo[as.character(mat_kmclust$cluster)])
    mat_dt = mat_dt[order(cluster_ordered),
                    list(id__ = mat_name, group__ = cluster_ordered)]
    return(mat_dt)
}


#' perform kmeans clustering on matrix rows and return reordered matrix along
#' with order matched cluster assignments clusters are sorted using hclust on
#' centers the contents of each cluster are sorted using hclust
#' @param mat A wide format matrix
#' @param nclust the number of clusters
#' @param within_order_strategy one of "hclust", "sort", "right", "left",
#'   "reverse".  If "hclust", hierarchical clustering will be used. If "sort", a
#'   simple decreasing sort of rosSums.  If "left", will atttempt to put high
#'   signal on left ("right" is opposite).  If "reverse" reverses existing order
#'   (should only be used after meaningful order imposed).
#' @param centroids optional matrix with same columns as mat and one centroid
#'   per row to base clusters off of.  Overrides any setting to nclust. Default
#'   of NULL results in randomly initialized k-means.
#' @param manual_mapping optional named vector manually specififying cluster
#'   assignments. names should be item ids and values should be cluster names
#'   the items are assigned to. Default of NULL allows clustering to proceed.
#' @param iter.max Number of max iterations to allow for k-means. Default is 30.
#'
#' @export
#' @importFrom stats  hclust dist
#' @return data.table with 2 columns of cluster info. id__ column corresponds with
#'   input matrix rownames and is sorted within each cluster using hierarchical
#'   clusering group__ column indicates cluster assignment
#' @examples
#' dt = data.table::copy(CTCF_in_10a_profiles_dt)
#' mat = data.table::dcast(dt, id ~ sample + x, value.var = "y" )
#' rn = mat$id
#' mat = as.matrix(mat[,-1])
#' rownames(mat) = rn
#' clust_dt = clusteringKmeansNestedHclust(mat, nclust = 3)
#' clust_dt
clusteringKmeansNestedHclust = function(mat,
                                        nclust,
                                        within_order_strategy = valid_sort_strategies[2],
                                        centroids = NULL,
                                        manual_mapping = NULL,
                                        iter.max = 30) {
    stopifnot(is.numeric(nclust))
    stopifnot(within_order_strategy %in% valid_sort_strategies)
    group__ = id__ = within_o = NULL#declare binding for data.table

    if(nclust > (nrow(mat)-1) & is.null(centroids)){
        nclust = nrow(mat)-1
        message("nclust too high for number of items. Reducing to ", nclust)
    }

    if(is.null(manual_mapping)){
        mat_dt = clusteringKmeans(mat, nclust, centroids, iter.max)
    }else{
        manual_mapping
        mat_dt = data.table(id__ = names(manual_mapping), group__ = manual_mapping)
        mat_dt = mat_dt[order(group__)]
    }

    within_clust_sort.mat_dt(mat_dt = mat_dt,
                             mat = mat,
                             within_order_strategy = within_order_strategy)
}


within_clust_sort.mat_dt = function(mat_dt, mat, within_order_strategy){
    group__ = id__ = within_o = NULL
    stopifnot(within_order_strategy %in% valid_sort_strategies)
    mat_dt$within_o = as.numeric(-1)
    for (i in unique(mat_dt$group__)) {
        cmat = mat[mat_dt[group__ == i, id__], , drop = FALSE]
        if(within_order_strategy == "hclust"){
            if (nrow(cmat) > 2) {
                mat_dt[group__ == i, ]$within_o =
                    stats::hclust(stats::dist((cmat)))$order
            } else {
                mat_dt[group__ == i, ]$within_o = seq_len(nrow(cmat))
            }
        }else if(within_order_strategy == "sort"){
            mat_dt[group__ == i, ]$within_o = rank(-rowSums(cmat))
        }else if(within_order_strategy == "left"){
            mat_dt[group__ == i, ]$within_o = rank(apply(cmat, 1, function(x){weighted.mean(seq_along(x), x)}))
        }else if(within_order_strategy == "right"){
            mat_dt[group__ == i, ]$within_o = rank(-apply(cmat, 1, function(x){weighted.mean(seq_along(x), x)}))
        }else if(within_order_strategy == "reverse"){
            mat_dt[group__ == i, ]$within_o = rev(seq_len(nrow(cmat)))
        }
    }
    mat_dt = mat_dt[order(within_o), ][order(group__), ]
    data.table::set(mat_dt, j = "id__", value = factor(mat_dt$id__, levels = mat_dt$id__))
    mat_dt$within_o = NULL
    if(!is.factor(mat_dt$group__)){
        data.table::set(mat_dt, j = "group__", value = factor(mat_dt$group__, levels = unique(mat_dt$group__)))
    }
    mat_dt
}

#' within_clust_sort
#'
#' Without modifying cluster assignments, modify the order of rows within each
#' cluster based on within_order_strategy.
#'
#' This is particularly useful when you want to sort within each cluster by a
#' different variable from cluster assignment. Also if you've imported cluster
#' assigments but want to sort within each for the new data for a prettier
#' heatmap.
#'
#' TODO refactor shared code with clusteringKmeansNestedHclust
#'
#' @param clust_dt data.table output from \code{\link{ssvSignalClustering}}
#' @param row_ variable name mapped to row, likely id or gene name for ngs data.
#'   Default is "id" and works with ssvFetch* output.
#' @param column_ varaible mapped to column, likely bp position for ngs data.
#'   Default is "x" and works with ssvFetch* output.
#' @param fill_ numeric variable to map to fill. Default is "y" and works with
#'   ssvFetch* output.
#' @param facet_ variable name to facet horizontally by. Default is "sample" and
#'   works with ssvFetch* output. Set to "" if data is not facetted.
#' @param cluster_ variable name to use for cluster info. Default is
#'   "cluster_id".
#' @param clustering_col_min numeric minimum for col range considered when
#'   clustering, default in -Inf
#' @param clustering_col_max numeric maximum for col range considered when
#'   clustering, default in Inf
#' @param within_order_strategy one of "hclust", "sort", "right", "left",
#'   "reverse".  If "hclust", hierarchical clustering will be used. If "sort", a
#'   simple decreasing sort of rosSums.  If "left", will atttempt to put high
#'   signal on left ("right" is opposite).  If "reverse" reverses existing order
#'   (should only be used after meaningful order imposed).
#' @param dcast_fill value to supply to dcast fill argument. default is NA.
#'
#' @return data.table matching input clust_dt save for the reassignment of
#'   levels of row_ variable.
#' @export
#'
#' @examples
#' #clustering by relative value per region does a good job highlighting changes
#' #however, when then plotting raw values the order within clusters is not smooth
#' #this is a good situation to apply a separate sort within clusters.
#' prof_dt = CTCF_in_10a_profiles_dt
#' prof_dt = append_ynorm(prof_dt)
#' prof_dt[, y_relative := y_norm / max(y_norm), list(id)]
#'
#' clust_dt = ssvSignalClustering(prof_dt, fill_ = "y_relative")
#' clust_dt.sort = within_clust_sort(clust_dt)
#'
#' cowplot::plot_grid(
#'   ssvSignalHeatmap(clust_dt) + labs(title = "clustered by relative, sorted by relative"),
#'   ssvSignalHeatmap(clust_dt.sort) + labs(title = "clustered by relative, sorted by raw value")
#' )
#'
within_clust_sort = function(clust_dt,
                             row_ = "id",
                             column_ = "x",
                             fill_ = "y",
                             facet_ = "sample",
                             cluster_ = "cluster_id",
                             within_order_strategy = c("hclust", "sort", "left", "right")[2],
                             clustering_col_min = -Inf,
                             clustering_col_max = Inf,
                             dcast_fill = NA) {
    output_GRanges = FALSE
    if(is(clust_dt, "GRanges")){
        clust_dt = data.table::as.data.table(clust_dt)
        output_GRanges = TRUE
    }else{
        clust_dt = data.table::copy(clust_dt)
    }
    stopifnot(is.data.table(clust_dt))
    mat = make_clustering_matrix(
        clust_dt,
        row_ = row_,
        column_ = column_,
        fill_ = fill_,
        facet_ = facet_,
        max_rows = Inf,
        max_cols = Inf,
        clustering_col_min = clustering_col_min,
        clustering_col_max = clustering_col_max,
        dcast_fill = dcast_fill
    )


    group__ = id__ = within_o = NULL#declare binding for data.table

    mat_dt = unique(clust_dt[, list(id__ = get(row_), group__ = get(cluster_))])
    mat_dt = within_clust_sort.mat_dt(mat_dt, mat, within_order_strategy = within_order_strategy)
    #repair var names
    data.table::set(clust_dt, j = row_, value = factor(clust_dt[[row_]], levels = levels(mat_dt$id__)))
    return(clust_dt[])
}

#' copy_clust_info
#'
#' @param target A data.table or GRanges returned from ssvFetch*, the target to
#'   which cluster info will be added.
#' @param to_copy A data.table or GRanges returned from ssvSignalClustering,
#'   from which to copy cluster if.
#' @param row_ variable name mapped to row, likely id or gene name for ngs data.
#'   Default is "id" and works with ssvFetch* output.
#' @param cluster_ variable name to use for cluster info. Default is
#'   "cluster_id".
#'
#' @return data.table or GRanges (whichever target is) containing row order and
#'   cluster assignment derived from to_copy. Suitable for ssvSignalHeatmap and
#'   related functions.
#' @export
#'
#' @examples
#' #this takes cluster info from signal and applies to peak hits to
#' #create a heatmap of peak hits clustered by signal.
#' clust_dt1 = ssvSignalClustering(CTCF_in_10a_profiles_dt)
#' peak_hit_gr = ssvFetchGRanges(
#'   CTCF_in_10a_narrowPeak_grs,
#'   qgr = CTCF_in_10a_overlaps_gr
#' )
#' peak_hit_gr.clust = copy_clust_info(peak_hit_gr, clust_dt1)
#' peak_hit_gr.clust$hit = peak_hit_gr.clust$y > 0
#' ssvSignalHeatmap(peak_hit_gr.clust, fill_ = "hit") +
#'   scale_fill_manual(values = c("FALSE" = "gray90", "TRUE" = "black"))
copy_clust_info = function(target, to_copy, row_ = "id", cluster_ = "cluster_id"){
    output_GRanges = FALSE
    if(is(target, "GRanges")){
        target = data.table::as.data.table(target)
        output_GRanges = TRUE
    }else{
        target = data.table::copy(target)
    }
    to_copy = data.table::as.data.table(to_copy)

    stopifnot(!is.null(target[[row_]]))
    stopifnot(!is.null(to_copy[[row_]]))
    stopifnot(!is.null(to_copy[[cluster_]]))
    stopifnot(length(intersect(target[[row_]], to_copy[[row_]])) > 0)

    if(!all(target[[row_]] %in% to_copy[[row_]])){
        warning("Some '", row_, "' items missing from to_copy, items will be lost from target.")
    }

    target[[cluster_]] = NULL
    target = merge(target, unique(to_copy[, c(row_, cluster_), with = FALSE]), by = row_)
    data.table::set(target, j = row_, value = factor(target[[row_]], levels = levels(to_copy[[row_]])))
    if(output_GRanges){
        target = GenomicRanges::GRanges(target)
    }
    target
}


#' reorder_clusters_stepdown
#'
#' Attempts to reorder clusters so that rows with highest signal on the left
#' relative to the right appear at the top. Signal should have a roughly
#' diagonal pattern in a "stepdown" pattern.
#'
#' This can be down by column (step_by_column = TRUE) which averages across
#' facets.  By facet (step_by_column = FALSE, step_by_facet = TRUE) which
#' averages all columns per facet. Or both column and facet (step_by_column =
#' TRUE, step_by_facet = TRUE), which does no averaging so it looks at the full
#' matrix as plotted.
#'
#' @param clust_dt data.table output from \code{\link{ssvSignalClustering}}
#' @param row_ variable name mapped to row, likely id or gene name for ngs data.
#'   Default is "id" and works with ssvFetch* output.
#' @param column_ varaible mapped to column, likely bp position for ngs data.
#'   Default is "x" and works with ssvFetch* output.
#' @param fill_ numeric variable to map to fill. Default is "y" and works with
#'   ssvFetch* output.
#' @param facet_ variable name to facet horizontally by. Default is "sample" and
#'   works with ssvFetch* output. Set to "" if data is not facetted.
#' @param cluster_ variable name to use for cluster info. Default is
#'   "cluster_id".
#' @param reapply_cluster_names If TRUE, clusters will be renamed according to
#'   new order instead of their original names. Default is TRUE.
#' @param step_by_column If TRUE, column is considered for left-right cluster
#'   balance. Default is TRUE.
#' @param step_by_facet If TRUE, facet is considered for left-right cluster
#'   balance. Default is FALSE.
#'
#' @return data.table as output from \code{\link{ssvSignalClustering}}
#' @export
#'
#' @examples
#' clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt, nclust = 10)
#' new_dt = reorder_clusters_stepdown(clust_dt)
#' cowplot::plot_grid(
#'     ssvSignalHeatmap(clust_dt),
#'     ssvSignalHeatmap(new_dt)
#' )
reorder_clusters_stepdown = function(clust_dt,
                                     row_ = "id",
                                     column_ = "x",
                                     fill_ = "y",
                                     facet_ = "sample",
                                     cluster_ = "cluster_id",
                                     reapply_cluster_names = TRUE,
                                     step_by_column = TRUE,
                                     step_by_facet = FALSE
){
    new_dt = copy(clust_dt)
    if(facet_ != ""){
        if(!is.factor(new_dt[[facet_]])){
            new_dt[[facet_]] = factor(new_dt[[facet_]], levels = unique(new_dt[[facet_]]))
        }
    }
    agg_dt = new_dt[, list(y = mean(get(fill_))), c(cluster_, column_, ifelse(facet_ == "", character(), facet_))]
    agg_dt = agg_dt[order(get(column_))][order(get(facet_))]

    if(step_by_column & step_by_facet){
        wide_dt = dcast(agg_dt,
                        paste0(cluster_, "~", ifelse(facet_ == "", character(), paste0(facet_, "+")), column_),
                        value.var = fill_,
                        fun.aggregate = mean)
    }else if(step_by_column){
        wide_dt = dcast(agg_dt,
                        paste0(cluster_, "~", column_),
                        value.var = fill_,
                        fun.aggregate = mean)
    }else if(step_by_facet){
        wide_dt = dcast(agg_dt,
                        paste0(cluster_, "~", ifelse(facet_ == "", character(), facet_)),
                        value.var = fill_,
                        fun.aggregate = mean)
    }


    mat = as.matrix(wide_dt[,-1])
    rownames(mat) = wide_dt[[cluster_]]

    wm_o = sort(apply(mat, 1, function(x){
        weighted.mean(seq_along(x), x)
    }))
    reorder_clusters_manual(new_dt,
                            manual_order = names(wm_o),
                            cluster_ = cluster_, row_ = row_)
}



#' reorder_clusters_hclust
#'
#' Applies hierarchical clustering to centroids of clusters to reorder.
#'
#' @param clust_dt data.table output from \code{\link{ssvSignalClustering}}
#' @param hclust_result hclust result returned by a previous call of this function with
#'   identical paramters when return_hclust = TRUE.
#' @param row_ variable name mapped to row, likely id or gene name for ngs data.
#'   Default is "id" and works with ssvFetch* output.
#' @param column_ varaible mapped to column, likely bp position for ngs data.
#'   Default is "x" and works with ssvFetch* output.
#' @param fill_ numeric variable to map to fill. Default is "y" and works with
#'   ssvFetch* output.
#' @param facet_ variable name to facet horizontally by. Default is "sample" and
#'   works with ssvFetch* output. Set to "" if data is not facetted.
#' @param cluster_ variable name to use for cluster info. Default is
#'   "cluster_id".
#' @param reapply_cluster_names If TRUE, clusters will be renamed according to
#'   new order instead of their original names. Default is TRUE.
#' @param return_hclust If TRUE, return the result of hclust instead of the
#'   reordered clustering data.table. Default is FALSE.  Ignored if
#'   hclust_result is supplied.
#'
#' @return data.table as output from \code{\link{ssvSignalClustering}}
#' @export
#'
#' @examples
#' clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt, nclust = 10)
#' new_dt = reorder_clusters_hclust(clust_dt)
#' cowplot::plot_grid(
#'     ssvSignalHeatmap(clust_dt),
#'     ssvSignalHeatmap(new_dt)
#' )
reorder_clusters_hclust = function(clust_dt,
                                   hclust_result = NULL,
                                   row_ = "id",
                                   column_ = "x",
                                   fill_ = "y",
                                   facet_ = "sample",
                                   cluster_ = "cluster_id",
                                   reapply_cluster_names = TRUE,
                                   return_hclust = FALSE){
    new_dt = copy(clust_dt)
    if(is.null(hclust_result)){
        agg_dt = new_dt[, list(y = mean(get(fill_))), c(cluster_, column_, ifelse(facet_ == "", character(), facet_))]
        setnames(agg_dt, "y", fill_)
        wide_dt = dcast(agg_dt, paste0(cluster_, "~", column_, ifelse(facet_ == "", character(), paste0("+", facet_))), value.var = fill_)
        mat = as.matrix(wide_dt[,-1])
        rownames(mat) = wide_dt[[cluster_]]
        hclust_result = hclust(dist(mat))
        if(return_hclust){
            return(hclust_result)
        }
    }
    clust_lev = levels(new_dt[[cluster_]])[hclust_result$order]
    reorder_clusters_manual(clust_dt,
                            manual_order = clust_lev,
                            cluster_ = cluster_,
                            row_ = row_,
                            reapply_cluster_names = reapply_cluster_names)
}
#


#' reorder_clusters_manual
#'
#' Manually applies a new order (top to bottom) for cluster using the result of
#' ssvSignalClustering.
#'
#' @param clust_dt data.table output from \code{\link{ssvSignalClustering}}
#' @param manual_order New order for clusters  Does not need to include all
#'   clusters.  Any colors not included will be at the bottom in their original
#'   order.
#' @param row_ variable name mapped to row, likely id or gene name for ngs data.
#'   Default is "id" and works with ssvFetch* output.
#' @param cluster_ variable name to use for cluster info. Default is
#'   "cluster_id".
#' @param reapply_cluster_names If TRUE, clusters will be renamed according to
#'   new order instead of their original names. Default is TRUE.
#'
#' @return data.table as output from \code{\link{ssvSignalClustering}}
#' @export
#'
#' @examples
#' clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt, nclust = 3)
#' new_dt = reorder_clusters_manual(clust_dt = clust_dt, manual_order = 2)
#' cowplot::plot_grid(
#'     ssvSignalHeatmap(clust_dt),
#'     ssvSignalHeatmap(new_dt)
#' )
reorder_clusters_manual = function(clust_dt,
                                   manual_order,
                                   row_ = "id",
                                   cluster_ = "cluster_id",
                                   reapply_cluster_names = TRUE){
    stopifnot(manual_order %in% clust_dt[[cluster_]])
    new_dt = copy(clust_dt[order(get(row_))])
    if(!is.factor(new_dt[[cluster_]])){
        new_dt[[cluster_]] = factor(new_dt[[cluster_]], levels = rev(unique(new_dt[[cluster_]])))
    }
    manual_levels = c(manual_order, setdiff(levels(new_dt[[cluster_]]), manual_order))
    new_dt[[cluster_]] = factor(new_dt[[cluster_]], levels = manual_levels)
    new_dt = new_dt[order(get(cluster_))]
    new_dt[[row_]] = factor(new_dt[[row_]], levels = unique(new_dt[[row_]]))
    if(reapply_cluster_names){
        levels(new_dt[[cluster_]]) = seq_along(unique(new_dt[[cluster_]]))
    }
    new_dt
}



#' merge_clusters
#'
#' @param clust_dt data.table output from \code{\link{ssvSignalClustering}}
#' @param to_merge Clusters to merge. Must be items in clust_dt variable defined
#'   by cluster_ parameter.
#' @param row_ variable name mapped to row, likely id or gene name for ngs data.
#'   Default is "id" and works with ssvFetch* output.
#' @param cluster_ variable name to use for cluster info. Default is
#'   "cluster_id".
#' @param reapply_cluster_names If TRUE, clusters will be renamed according to
#'   new order instead of their original names. Default is TRUE.
#'
#' @return data.table as output from \code{\link{ssvSignalClustering}}
#' @export
#'
#' @examples
#' set.seed(0)
#' clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt, nclust = 6)
#' ssvSignalHeatmap(clust_dt)
#' agg_dt = clust_dt[, list(y = mean(y)), list(x, cluster_id, sample)]
#' ggplot(agg_dt, aes(x = x, y = y, color = sample)) +
#'   geom_path() +
#'   facet_grid(cluster_id~.)
#'
#' to_merge = c(2, 3, 5)
#' # debug(merge_clusters)
#' new_dt = merge_clusters(clust_dt, c(2, 3, 5), reapply_cluster_names = FALSE)
#' new_dt.relabel = merge_clusters(clust_dt, c(2, 3, 5), reapply_cluster_names = TRUE)
#' new_dt.relabel.sort = within_clust_sort(new_dt.relabel, within_order_strategy = "sort")
#'
#' table(clust_dt$cluster_id)
#' table(new_dt$cluster_id)
#'
#' cowplot::plot_grid(
#'   ssvSignalHeatmap(clust_dt) + labs(title = "original"),
#'   ssvSignalHeatmap(new_dt) + labs(title = "2,3,5 merged"),
#'   ssvSignalHeatmap(new_dt.relabel) + labs(title = "2,3,5 merged, renumbered"),
#'   ssvSignalHeatmap(new_dt.relabel.sort) + labs(title = "2,3,5 merged, renumbered and sorted")
#' )
#'
#'
#'
#'
merge_clusters = function(clust_dt,
                          to_merge,
                          row_ = "id",
                          cluster_ = "cluster_id",
                          reapply_cluster_names = TRUE){
    to_merge = as.character(to_merge)
    clust_lev = levels(clust_dt[[cluster_]])
    stopifnot(all(to_merge %in% clust_lev))
    first_clust = min(which(clust_lev %in% to_merge))
    new_order = union(clust_lev[seq(1, first_clust)], to_merge)

    new_dt = reorder_clusters_manual(clust_dt, manual_order = new_order, row_ = row_, cluster_ = cluster_, reapply_cluster_names = FALSE)
    new_name = paste(to_merge, collapse = "-")

    data.table::set(new_dt, i = which(new_dt[[cluster_]] %in% to_merge), j = cluster_, value = new_name)
    # new_dt[[row_]] = factor(new_dt[[row_]], levels = unique(new_dt[[row_]]))
    new_dt[[cluster_]] = factor(new_dt[[cluster_]], levels = unique(new_dt[[cluster_]]))
    if(reapply_cluster_names){
        levels(new_dt[[cluster_]]) = seq_along(unique(new_dt[[cluster_]]))
    }
    new_dt

}

#' reverse_clusters
#'
#' @param clust_dt data.table output from \code{\link{ssvSignalClustering}}
#' @param row_ variable name mapped to row, likely id or gene name for ngs data.
#'   Default is "id" and works with ssvFetch* output.
#' @param column_ varaible mapped to column, likely bp position for ngs data.
#'   Default is "x" and works with ssvFetch* output.
#' @param fill_ numeric variable to map to fill. Default is "y" and works with
#'   ssvFetch* output.
#' @param facet_ variable name to facet horizontally by. Default is "sample" and
#'   works with ssvFetch* output. Set to "" if data is not facetted.
#' @param cluster_ variable name to use for cluster info. Default is
#'   "cluster_id".
#' @param reverse_rows_within If TRUE, rows within clusters will be reversed as well. Default is TRUE.
#' @param reapply_cluster_names If TRUE, clusters will be renamed according to
#'   new order instead of their original names. Default is TRUE.
#'
#' @return data.table as output from \code{\link{ssvSignalClustering}}
#'
#' @export
#'
#' @examples
#' set.seed(0)
#' clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt, nclust = 3)
#' rev_dt = reverse_clusters(clust_dt)
#' rev_dt.no_relabel = reverse_clusters(clust_dt, reapply_cluster_names = FALSE)
#' rev_dt.not_rows = reverse_clusters(clust_dt, reverse_rows_within = FALSE)
#' cowplot::plot_grid(nrow = 1,
#'   ssvSignalHeatmap(clust_dt) + labs(title = "original"),
#'   ssvSignalHeatmap(rev_dt) + labs(title = "reversed"),
#'   ssvSignalHeatmap(rev_dt.no_relabel) + labs(title = "reversed, no relabel"),
#'   ssvSignalHeatmap(rev_dt.not_rows) + labs(title = "reversed, not rows")
#' )
reverse_clusters = function(clust_dt,
                            row_ = "id",
                            column_ = "x",
                            fill_ = "y",
                            facet_ = "sample",
                            cluster_ = "cluster_id",
                            reverse_rows_within = TRUE,
                            reapply_cluster_names = TRUE){
    new_dt = reorder_clusters_manual(
        clust_dt,
        manual_order = rev(levels(clust_dt[[cluster_]])),
        row_ = row_,
        cluster_ = cluster_,
        reapply_cluster_names = reapply_cluster_names
    )
    if(reverse_rows_within){
        new_dt = within_clust_sort(
            new_dt,
            row_ = row_,
            column_ = column_,
            fill_ = fill_,
            facet_ = facet_,
            cluster_ = cluster_, within_order_strategy = "reverse"
        )
    }
    new_dt
}

#' split_cluster
#'
#' Splits one specified cluster in number of new clusters determined by nclust
#'
#' @param clust_dt data.table output from \code{\link{ssvSignalClustering}}
#' @param to_split Cluster to split.
#' @param nclust Number of new clusters to create.
#' @param row_ variable name mapped to row, likely id or gene name for ngs data.
#'   Default is "id" and works with ssvFetch* output.
#' @param column_ varaible mapped to column, likely bp position for ngs data.
#'   Default is "x" and works with ssvFetch* output.
#' @param fill_ numeric variable to map to fill. Default is "y" and works with
#'   ssvFetch* output.
#' @param facet_ variable name to facet horizontally by. Default is "sample" and
#'   works with ssvFetch* output. Set to "" if data is not facetted.
#' @param cluster_ variable name to use for cluster info. Default is
#'   "cluster_id".
#' @param reapply_cluster_names If TRUE, clusters will be renamed according to
#'   new order instead of their original names. Default is TRUE.
#'
#' @return data.table as output from \code{\link{ssvSignalClustering}}
#' @export
#'
#' @examples
#' set.seed(0)
#' clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt, nclust = 3)
#' split_dt = split_cluster(clust_dt, to_split = 2, nclust = 3)
#' split_dt.no_rename = split_cluster(
#'   clust_dt,
#'   to_split = 2,
#'   nclust = 3,
#'   reapply_cluster_names = FALSE
#' )
#' cowplot::plot_grid(nrow = 1,
#'   ssvSignalHeatmap(clust_dt),
#'   ssvSignalHeatmap(split_dt),
#'   ssvSignalHeatmap(split_dt.no_rename)
#' )
#'
split_cluster = function(clust_dt,
                         to_split,
                         nclust = 2,
                         row_ = "id",
                         column_ = "x",
                         fill_ = "y",
                         facet_ = "sample",
                         cluster_ = "cluster_id",
                         reapply_cluster_names = TRUE){
    to_split = as.character(to_split)
    stopifnot(to_split %in% clust_dt[[cluster_]])
    split_dt = clust_dt[clust_dt[[cluster_]] == to_split]
    split_dt = ssvSignalClustering(split_dt,
                                   nclust = nclust,
                                   row_ = row_,
                                   column_ = column_,
                                   fill_ = fill_,
                                   facet_ = facet_,
                                   cluster_ = cluster_,
                                   max_rows = Inf,
                                   max_cols = Inf)
    new_lev = paste0(to_split, letters[seq_len(nclust)])
    levels(split_dt[[cluster_]]) = new_lev
    insert_at = which(levels(clust_dt[[cluster_]]) == to_split)

    new_dt =rbind(split_dt, clust_dt[clust_dt[[cluster_]] != to_split])
    new_dt$cluster_id = droplevels(new_dt$cluster_id)
    c_lev = levels(clust_dt[[cluster_]])
    new_lev = c(c_lev[seq_along(c_lev) < insert_at], new_lev, c_lev[seq_along(c_lev) > insert_at])

    stopifnot(setequal(new_lev, levels(new_dt$cluster_id)))
    new_dt$cluster_id = factor(new_dt$cluster_id, levels = new_lev)
    reorder_clusters_manual(new_dt,
                            manual_order = new_lev,
                            row_ = row_,
                            cluster_ = cluster_,
                            reapply_cluster_names = reapply_cluster_names)
}
