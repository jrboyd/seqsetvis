
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
#' @param row_ variable name mapped to row, likely id or gene name for ngs data. Default is "id" and works with ssvFetch* output.
#' @param column_ varaible mapped to column, likely bp position for ngs data. Default is "x" and works with ssvFetch* output.
#' @param fill_ numeric variable to map to fill. Default is "y" and works with ssvFetch* output.
#' @param facet_ variable name to facet horizontally by. Default is "sample" and works with ssvFetch* output. Set to "" if data is not facetted.
#' @param cluster_ variable name to use for cluster info. Default is "cluster_id".
#' @param max_rows for speed rows are sampled to 500 by default, use Inf to plot
#'   full data
#' @param max_cols for speed columns are sampled to 100 by default, use Inf to
#'   plot full data
#' @param clustering_col_min numeric minimum for col range considered when
#'   clustering, default in -Inf
#' @param clustering_col_max numeric maximum for col range considered when
#'   clustering, default in Inf
#' @param within_order_strategy one of "hclust" or "sort".  if hclust,
#'   hierarchical clustering will be used. if sort, a simple decreasing sort of
#'   rosSums.
#' @param dcast_fill value to supply to dcast fill argument. default is NA.
#' @param iter.max Number of max iterations to allow for k-means. Default is 30.
#'
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
ssvSignalClustering = function(bw_data, nclust = NULL,
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
                               within_order_strategy = c("hclust", "sort")[2],
                               dcast_fill = NA,
                               iter.max = 30){
    message("clustering...")
    id = xbp = x = to_disp = y = hit = val = y = y_gap = group =  NULL#declare binding for data.table
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
        fill_ = fill_, facet_ = facet_,
        max_rows = max_rows,
        max_cols = max_cols,
        clustering_col_min = clustering_col_min,
        clustering_col_max = clustering_col_max,
        dcast_fill = dcast_fill

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
    plot_dt = plot_dt[get(row_) %in% rclusters[["id"]]]
    plot_dt[[row_]] = factor(plot_dt[[row_]], levels = rclusters[["id"]])
    data.table::setkey(rclusters, id)
    plot_dt[[cluster_]] = rclusters[list(plot_dt[[row_]]), group]
    if(output_GRanges){
        plot_dt = GRanges(plot_dt)
    }
    return(plot_dt)
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
#' @return data.table with group variable indicating cluster membership and id
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
#' dt = merge(dt, clust_dt)
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
        if(nrow(unique(mat)) < nclust){
            nclust = nrow(unique(mat))
            warning("Reducing nclust to ", nclust,
                    " - maximum number of clusters allowed due to ",
                    "low uniqueness.")
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
                    list(id = mat_name, group = cluster_ordered)]
    return(mat_dt)
}


#' perform kmeans clustering on matrix rows and return reordered matrix along
#' with order matched cluster assignments clusters are sorted using hclust on
#' centers the contents of each cluster are sorted using hclust
#' @param mat A wide format matrix
#' @param nclust the number of clusters
#' @param within_order_strategy one of "hclust" or "sort".  if hclust,
#'   hierarchical clustering will be used. if sort, a simple decreasing sort of
#'   rosSums.
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
#' @return data.table with 2 columns of cluster info. id column corresponds with
#'   input matrix rownames and is sorted within each cluster using hierarchical
#'   clusering group column indicates cluster assignment
#' @examples
#' dt = data.table::copy(CTCF_in_10a_profiles_dt)
#' mat = data.table::dcast(dt, id ~ sample + x, value.var = "y" )
#' rn = mat$id
#' mat = as.matrix(mat[,-1])
#' rownames(mat) = rn
#' clust_dt = clusteringKmeansNestedHclust(mat, nclust = 3)
#' dt = merge(dt, clust_dt)
#' dt$id = factor(dt$id, levels = clust_dt$id)
#' dt[order(id)]
clusteringKmeansNestedHclust = function(mat, nclust, within_order_strategy = c("hclust", "sort")[2], centroids = NULL, manual_mapping = NULL, iter.max = 30) {
    stopifnot(is.numeric(nclust))
    stopifnot(within_order_strategy %in% c("hclust", "sort"))
    group = id = within_o = NULL#declare binding for data.table
    if(is.null(manual_mapping)){
        mat_dt = clusteringKmeans(mat, nclust, centroids, iter.max)
    }else{
        manual_mapping
        mat_dt = data.table(id = names(manual_mapping), group = manual_mapping)
        mat_dt = mat_dt[order(group)]
    }

    mat_dt$within_o = as.numeric(-1)
    for (i in unique(mat_dt$group)) {
        cmat = mat[mat_dt[group == i, id], , drop = FALSE]
        if(within_order_strategy == "hclust"){
            if (nrow(cmat) > 2) {
                mat_dt[group == i, ]$within_o =
                    stats::hclust(stats::dist((cmat)))$order
            } else {
                mat_dt[group == i, ]$within_o = seq_len(nrow(cmat))
            }
        }else if(within_order_strategy == "sort"){
            mat_dt[group == i, ]$within_o = rank(-rowSums(cmat))
        }


    }
    mat_dt = mat_dt[order(within_o), ][order(group), ]
    mat_dt$id = factor(mat_dt$id, levels = mat_dt$id)
    mat_dt$within_o = NULL
    if(!is.factor(mat_dt$group)){
        mat_dt$group = factor(mat_dt$group)
    }
    return(mat_dt)
}

#' within_clust_sort
#'
#' Without modifying cluster assignments, modify the order of rows within each cluster based on within_order_strategy.
#'
#' This is particularly useful when you want to sort within each cluster by a different variable from cluster assignment. Also if you've imported cluster assigments but want to sort within each for the new data for a prettier heatmap.
#'
#' TODO refactor shared code with clusteringKmeansNestedHclust
#'
#' @param clust_dt data.table output from \code{\link{ssvSignalClustering}}
#' @param row_ variable name mapped to row, likely id or gene name for ngs data. Default is "id" and works with ssvFetch* output.
#' @param column_ varaible mapped to column, likely bp position for ngs data. Default is "x" and works with ssvFetch* output.
#' @param fill_ numeric variable to map to fill. Default is "y" and works with ssvFetch* output.
#' @param facet_ variable name to facet horizontally by. Default is "sample" and works with ssvFetch* output. Set to "" if data is not facetted.
#' @param cluster_ variable name to use for cluster info. Default is "cluster_id".
#' @param clustering_col_min numeric minimum for col range considered when
#'   clustering, default in -Inf
#' @param clustering_col_max numeric maximum for col range considered when
#'   clustering, default in Inf
#' @param within_order_strategy one of "hclust" or "sort".  if hclust,
#'   hierarchical clustering will be used. if sort, a simple decreasing sort of
#'   rosSums.
#' @param dcast_fill value to supply to dcast fill argument. default is NA.
#'
#' @return data.table matching input clust_dt save for the reassignment of levels of row_ variable.
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
                             within_order_strategy = c("hclust", "sort")[2],
                             clustering_col_min = -Inf,
                             clustering_col_max = Inf,
                             dcast_fill = NA) {
    output_GRanges = FALSE
    if(is(clust_dt, "GRanges")){
        clust_dt = data.table::as.data.table(clust_dt)
        output_GRanges = TRUE
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

    stopifnot(within_order_strategy %in% c("hclust", "sort"))
    group = id = within_o = NULL#declare binding for data.table

    mat_dt = unique(clust_dt[, list(id = get(row_), group = get(cluster_))])
    mat_dt$within_o = as.numeric(-1)
    for (i in unique(mat_dt$group)) {
        cmat = mat[mat_dt[group == i, id], , drop = FALSE]
        if(within_order_strategy == "hclust"){
            if (nrow(cmat) > 2) {
                mat_dt[group == i, ]$within_o =
                    stats::hclust(stats::dist((cmat)))$order
            } else {
                mat_dt[group == i, ]$within_o = seq_len(nrow(cmat))
            }
        }else if(within_order_strategy == "sort"){
            mat_dt[group == i, ]$within_o = rank(-rowSums(cmat))
        }
    }
    mat_dt = mat_dt[order(within_o), ][order(group), ]
    mat_dt$id = factor(mat_dt$id, levels = mat_dt$id)
    mat_dt$within_o = NULL
    if(!is.factor(mat_dt$group)){
        mat_dt$group = factor(mat_dt$group)
    }

    clust_dt$id = factor(clust_dt$id, levels = rev(levels(mat_dt$id)))

    return(clust_dt)
}
