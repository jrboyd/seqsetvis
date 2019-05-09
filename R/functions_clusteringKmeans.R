#' perform kmeans clustering on matrix rows and return reordered matrix along
#' with order matched cluster assignments. clusters are sorted using hclust on
#' centers
#'
#' @param mat numeric matrix to cluster
#' @param nclust the number of clusters
#' @param seed DEPRECATED.  Call set.seed() prior to this funciton to allow
#'   reproducibility.
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
clusteringKmeans = function(mat, nclust, seed = NULL) {
    stopifnot(is.numeric(nclust))
    cluster_ordered = mat_name = NULL#declare binding for data.table
    if(!is.null(seed)){
        warning("'seed' parameter is now deprecated. ",
                "Please call set.seed() yourself if needed.")
    }
    if(nrow(mat) <= nclust){
        nclust = nrow(mat)
        mat_kmclust = list(
            centers = mat,
            cluster = seq(nrow(mat))
        )
        names(mat_kmclust$cluster) = seq(nrow(mat))

    }else{
        mat_kmclust = stats::kmeans(mat, centers = nclust, iter.max = 30)
    }

    if(nclust < 3){
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
#' with order matched cluster assignments
#' clusters are sorted using hclust on centers
#' the contents of each cluster are sorted using hclust
#' @param mat A wide format matrix
#' @param nclust the number of clusters
#' @param seed passed to set.seed() to allow reproducibility
#' @export
#' @importFrom stats  hclust dist
#' @return data.table with 2 columns of cluster info.
#' id column corresponds with input matrix rownames and is sorted within
#' each cluster using hierarchical clusering
#' group column indicates cluster assignment
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
clusteringKmeansNestedHclust = function(mat, nclust, seed = NULL) {
    stopifnot(is.numeric(nclust))
    group = id = within_o = NULL#declare binding for data.table
    mat_dt = clusteringKmeans(mat, nclust)
    mat_dt$within_o = as.integer(-1)
    for (i in seq_along(nclust)) {
        cmat = mat[mat_dt[group == i, id], , drop = FALSE]
        if (nrow(cmat) > 2) {
            mat_dt[group == i, ]$within_o =
                stats::hclust(stats::dist((cmat)))$order
        } else {
            mat_dt[group == i, ]$within_o = seq_len(nrow(cmat))
        }

    }
    mat_dt = mat_dt[order(within_o), ][order(group), ]
    mat_dt$within_o = NULL
    return(mat_dt)
}

