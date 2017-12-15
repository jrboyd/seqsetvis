library(rtracklayer)
library(pbapply)
library(data.table)

#' Transforms set of GRanges to all have the same size.
#'
#' \code{centerFixedSizeGranges} First calculates the central coordinate of each
#' GRange in \code{grs} and extends in both direction by half of \code{fixed_size}
#'
#' @param grs Set of GRanges with incosistent and/or incorrect size
#' @param fixed_size The final width of each GRange returned.
#' @return Set of GRanges after resizing all input GRanges, either shortened
#' or lengthened as required to match \code{fixed_size}
centerFixedSizeGRanges = function(grs, fixed_size = 2000){
  m = round(start(grs) + width(grs) / 2)
  ext = round(fixed_size / 2)
  start(grs) = m - ext
  end(grs) = m + fixed_size - ext - 1
  return(grs)
}

#' Fetch values from a bigwig appropriate for heatmaps etc.
#'
#' \code{fetchWindowedBigwig} Gets values for each region of the query GRanges (\code{qgr}).
#' Values correspond to the center of each window of size \code{win_size}.  A tidy formatted data.table
#' object is returned suitable for plotting using ggplots.
#'
#' @param bw_file The bigwig file to read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @return A tidy formatted data.table containing fetched values.
#'
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved from \code{bw_file}
#'
fetchWindowedBigwig = function(bw_file, qgr, win_size = 50){
  if(!all(width(qgr) %% win_size == 0)){
    stop(paste("all widths of qgr are not evenly divisible by win_size,", win_size))
  }
  if(!all(width(qgr) == width(qgr)[1])){
    warning(paste("Widths vary between GRanges.  Before aggregating or co-plotting, recommend one of:
                  1) modifying input qgr such that all GRanges have same width
                  2) filtering output so that the same values of x are present for every region"))
  }
  # suppressWarnings({
    bw_gr = rtracklayer::import.bw(bw_file, which = qgr)
  # })
  windows = slidingWindows(qgr, width = win_size, step = win_size)
  if(is.null(qgr$id)){
    if(!is.null(names(qgr))){
      qgr$id = names(qgr)
    }else{
      qgr$id = paste0("region_", seq_along(qgr))
    }
  }
  names(windows) = qgr$id
  # print(object.size(windows), units = "GB")
  windows = unlist(windows)
  windows$id = names(windows)
  mid_gr = function(gr){
    start(gr) + floor((width(gr) - 1)/2)
  }
  mids = mid_gr(windows)
  start(windows) = mids
  end(windows) = mids
  olaps = suppressWarnings(data.table::as.data.table(findOverlaps(query = windows, subject = bw_gr)))
  #patch up missing/out of bound data with 0
  missing_idx = setdiff(1:length(windows), olaps$queryHits)
  if(length(missing_idx) > 0){
    olaps = rbind(olaps, data.table::data.table(queryHits = missing_idx, subjectHits = length(bw_gr) + 1))[order(queryHits)]
    bw_gr = c(bw_gr, GRanges(seqnames(bw_gr)[length(bw_gr)], IRanges(1, 1), score = 0))
  }
  #set FE and output
  # windows = windows[olaps$queryHits]
  windows$FE = bw_gr[olaps$subjectHits]$score
  bw_dt = data.table::as.data.table(windows)
  bw_dt[, x := start - min(start) + win_size / 2, by = id]
  bw_dt[, x := x - round(mean(x)), by = id ]
  shift = round(win_size / 2)
  bw_dt[, start := start - shift + 1]
  bw_dt[, end := end + win_size - shift]
  return(bw_dt)
}

#' applies a spline smoothing to a tidy data.table containing x and y values.
#'
#' \code{applySpline} Is intended for two-dimensional tidy data.tables, as retured by \code{fetchWindowedBigwig}
#'
#' @param dt a tidy data.table containing two-dimensional data
#' @param x_ the variable name of the x-values
#' @param y_ the variable name of the y-values
#' @param by_ optionally, any variables that provide grouping to the data.  see details.
#' @param n the number of interpolation points to use per input point, see \code{?spline}
#'
#' @return a newly derived data.table that is \code{n} times longer than original.
#'
#' @details by_ is quite powerful.  If \code{by_ = c("gene_id", "sample_id")}, splines
#' will be calculated individually for each gene in each sample. alternatively if \code{by_ = c("gene_id")}
#' @seealso \code{\link{fetchWindowedBigwig}}
applySpline = function(dt, x_ = "x", y_ = "y", by_ = "", n = 8, floor_at_0 = T, ...){
  sdt = dt[, spline(get(x_), get(y_), n = .N*n, ...), by = by_]
  if(floor_at_0){
    sdt[y < 0, y := 0]
  }
  colnames(sdt)[colnames(sdt) == "x"] = x_
  colnames(sdt)[colnames(sdt) == "y"] = y_
  return(sdt)
}

#' centers profile of x and y.  default is to center by region but across all samples.
#'
#' \code{centerAtMax} locates the coordinate x of the maximum in y and shifts x such that it is zero at max y.
#'
#' @param x_ the variable name of the x-values
#' @param y_ the variable name of the y-values
#' @param by_ optionally, any variables that provide grouping to the data.  see details.
#' @param view_size the size in \code{x_} to consider for finding the max of \code{y_}
#' @details by_ is quite powerful.  If \code{by_ = c("gene_id", "sample_id")}, splines
#' will be calculated individually for each gene in each sample. alternatively if \code{by_ = c("gene_id")}
centerAtMax = function(dt, x_ = "x", y_ = "y", by_ = "id", view_size, trim_to_valid = T){
  dt = copy(dt)
  xmax = round(view_size / 2)
  closestToZero = function(x){
    x[order(abs(x))][1]
  }
  dt[, ymax :=  max(get(y_)[abs(get(x_)) <= xmax]), by = by_]
  dt[, xsummit :=  closestToZero(get(x_)[get(y_) == ymax]), by = by_]
  dt[, xnew := get(x_) - xsummit]
  dt[, ymax := NULL]
  dt[, xsummit := NULL]
  if(trim_to_valid){
    dt = dt[xnew >= -xmax & xnew <= xmax]
  }
  set(dt, j = x_, value = dt$xnew)
  dt$xnew = NULL
  return(dt)
}

#extracts data.table, one row per qgr from bw_gr

ggBandedQuantiles = function(bw_dt, hsv_min = 0, hsv_max = .7,
                             n_quantile = 18, quantile_min = .05,
                             quantile_max = .95, is_centered = T, win_size = 50){
  #hsv_min = 0; hsv_max = .7; n_quantile = 18; quantile_min = .05; quantile_max = .95; is_centered = T; win_size = 50
  require(data.table)
  q2do = 0:n_quantile/n_quantile
  q2do = round(quantile_min + q2do * (quantile_max - quantile_min), digits = 3)
  dt = bw_dt[, .(qs = .(.(quantile(FE, q2do)))), by = x]
  dt = cbind(dt, as.data.table(t(sapply(dt$qs, function(x)x[[1]]))))
  dt$qs = NULL
  dtm = melt(dt, id.vars = "x")
  setkey(dtm, variable)
  q2do_str = paste0(q2do * 100, "%")
  dt_low = dtm[q2do_str[-length(q2do_str)], .(x, low_q = variable, low = value)]
  dt_high = dtm[q2do_str[-1], .(x, high_q = variable, high = value)]
  dt_c = cbind(dt_low, dt_high[, -1])
  dt_c = dt_c[, .(x, low, high, q_range = paste0(sub("%", "", low_q), "-", high_q))]
  q_o = unique(dt_c$q_range)
  dt_c$q_range = factor(dt_c$q_range, levels = q_o)
  if(!is.null(win_size)){
    dt_c[, xe := x + as.integer(win_size)]
    dt_c = melt(dt_c, id.vars = c("low", "high", "q_range"), value.name = "x")
    dt_c = dt_c[order(variable, decreasing = T)][order(x)][order(q_range)]
  }
  if(is.numeric(is_centered)){
    dt_c[, x := x - is_centered]
  }else if(is.logical(is_centered)){
    if(is_centered){
      dt_c[, x := x - median(x)]
    }
  }
  cols = rainbow(length(q_o), start = hsv_min, end = hsv_max)
  names(cols) = q_o
  dt_c$rn = 1:nrow(dt_c)
  dt_c[, c("q_low", "q_high") := tstrsplit(sub("%", "", q_range), split = "-", keep = 1:2), by = rn]
  dt_c[, q_num := (as.numeric(q_low) + as.numeric(q_high)) / 2]
  dt_c[, c("rn", "q_low", "q_high") := NULL]
  ggplot(dt_c) +
    geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = q_range)) +
    labs(fill = "quantile band",
         y = "FE", x = "bp",
         title = "Enrichment Profiles",
         subtitle = "aggregated by quantile range") +
    scale_fill_manual(values = cols) +
    scale_color_manual(breaks = q2do, palette = cols)

}

#' perform kmeans clustering on matrix rows and return reordered matrix along with order matched cluster assignments
#' clusters are sorted using hclust on centers
clusteringKmeans = function(mat, nclust, seed = 0){
  set.seed(seed)
  mat_kmclust = kmeans(mat, centers = nclust, iter.max = 30)
  center_o = hclust(dist(mat_kmclust$centers))$order
  center_reo = 1:length(center_o)
  names(center_reo) = center_o
  center_reo[as.character(mat_kmclust$cluster)]
  mat_dt = data.table(mat_name = names(mat_kmclust$cluster),
                      cluster = mat_kmclust$cluster,
                      cluster_ordered = center_reo[as.character(mat_kmclust$cluster)])
  mat_dt = mat_dt[order(cluster_ordered), .(id = mat_name, group = cluster_ordered)]
  return(mat_dt)
}


#' perform kmeans clustering on matrix rows and return reordered matrix along with order matched cluster assignments
#' clusters are sorted using hclust on centers
#' the contents of each cluster are sorted using hclust
#' @param mat A wide format matrix
clusteringKmeansNestedHclust = function(mat, nclust){
  mat_dt = clusteringKmeans(mat, nclust)
  mat_dt$within_o = as.integer(-1)
  for(i in 1:length(nclust)){
    cmat = mat[mat_dt[group == i, id], , drop = F]
    if(nrow(cmat) > 2){
      mat_dt[group == i, ]$within_o = hclust(dist((cmat)))$order
    }else{
      mat_dt[group == i, ]$within_o = 1:nrow(cmat)
    }

  }
  mat_dt = mat_dt[order(within_o),][order(group),]
  mat_dt$within_o = NULL
  return(mat_dt)
}
