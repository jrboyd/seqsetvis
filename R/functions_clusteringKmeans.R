#' perform kmeans clustering on matrix rows and return reordered matrix along
#' with order matched cluster assignments. clusters are sorted using hclust on
#' centers
#'
#' @param mat numeric matrix to cluster.
#' @param nclust the number of clusters.
#' @param centroids optional matrix with same columns as mat and one centroid
#'   per row to base clusters off of.  Overrides any setting to nclust. Default
#'   of NULL results in randomly initialized k-means.
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

#' Title
#' TODO refactor shared code with clusteringKmeansNestedHclust
#' @param clust_dt
#' @param row_
#' @param column_
#' @param fill_
#' @param facet_
#' @param cluster_
#' @param within_order_strategy
#' @param clustering_col_min
#' @param clustering_col_max
#' @param dcast_fill
#'
#' @return
#' @export
#'
#' @examples
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

    mat_dt = unique(clust_dt[, .(id = get(row_), group = get(cluster_))])
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
