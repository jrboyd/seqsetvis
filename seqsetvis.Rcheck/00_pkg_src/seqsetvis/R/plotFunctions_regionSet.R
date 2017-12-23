#' plot profiles from bigwigs
#'
#' @param bw_dt data.table of bigwig signal
#' @param y_ the variable name in bw_dt for y axis in plot
#' @param x_ the variable name in bw_dt for x axis in plot
#' @param by_ the variable name in bw_dt to facet on
#' @param hsv_reverse logical, should color scale be reversed? default FALSE
#' @param hsv_saturation numeric [0, 1] saturation for color scale. default 1
#' @param hsv_value numeric [0, 1] value for color scale. default 1
#' @param hsv_grayscale logical, if TRUE gray() is used instead of rainbow(). default FALSE
#' @param hsv_hue_min numeric [0, hsv_hue_max) hue min of color scale
#' @param hsv_hue_max numeric (hsv_hue_min, 1] hue max of color scale
#' @param symm_colors if TRUE, colorscale is symmetrical, default FALSE.
#' @param n_quantile number of evenly size quantile bins
#' @param quantile_min the lowest quantile start
#' @param quantile_max the highest quantile end
#'
#' @return ggplot object using ribbon plots to show quantile distributions
regionSetPlotBandedQuantiles = function(bw_dt, y_ = "FE", x_ = "x", by_ = "fake",
                                        hsv_reverse = F,
                                        hsv_saturation = 1, hsv_value = 1,
                                        hsv_grayscale = F,
                                        hsv_hue_min = 0, hsv_hue_max = 0.7, symm_colors = F,
                                        n_quantile = 18, quantile_min = 0.05, quantile_max = 0.95
) {
  variable = value = V1 = low = high = low_q = high_q = q_range = rn = q_num = q_low = q_high = NULL #declare binding for data.table
  q2do = 0:n_quantile/n_quantile
  q2do = round(quantile_min + q2do * (quantile_max - quantile_min), digits = 3)

  if(by_ == "fake"){
    bw_dt[[by_]] = T
    todo = T
  }else{
    todo = unique(bw_dt[, get(by_)])
  }
  all_q = lapply(todo, function(td){
    # for(td in todo){

    dt = bw_dt[get(by_) == td, .(qs = .(.(quantile(get(y_), q2do)))), by = x_]
    dt = cbind(dt, as.data.table(t(sapply(dt$qs, function(x) x[[1]]))))
    dt$qs = NULL
    dtm = melt(dt, id.vars = x_)
    setkey(dtm, variable)
    q2do_str = paste0(q2do * 100, "%")
    dt_low = dtm[q2do_str[-length(q2do_str)], .(get(x_), low_q = variable, low = value)]
    dt_high = dtm[q2do_str[-1], .(get(x_), high_q = variable, high = value)]
    dt_c = cbind(dt_low, dt_high[, -1])
    dt_c = dt_c[, .(V1, low, high, q_range = paste0(sub("%", "", low_q), "-", high_q))]
    q_o = unique(dt_c$q_range)
    dt_c$q_range = factor(dt_c$q_range, levels = q_o)

    dt_c$rn = 1:nrow(dt_c)
    dt_c[, `:=`(c("q_low", "q_high"), tstrsplit(sub("%", "", q_range), split = "-", keep = 1:2)), by = rn]
    dt_c[, `:=`(q_num, (as.numeric(q_low) + as.numeric(q_high))/2)]
    dt_c[, `:=`(c("rn", "q_low", "q_high"), NULL)]
    dt_c$facet_group = td
    return(dt_c)
  })
  plot_dt = rbindlist(all_q)
  q_o = unique(plot_dt$q_range)
  if(symm_colors){
    ncol = ceiling(length(q_o) / 2)
    if(hsv_grayscale){
      gray_vals = hsv_hue_min + (hsv_hue_max - hsv_hue_min)*((seq_len(ncol) - 1)/(ncol - 1))
      cols = gray(gray_vals)
    }else{
      cols = rainbow(ncol, start = hsv_hue_min, end = hsv_hue_max, s = hsv_saturation, v = hsv_value)
    }
    if(hsv_reverse) cols = rev(cols)
    ci = c(seq_len(ncol), rev(seq_len(ncol) - ifelse(length(q_o) %% 2 == 0, 0, 1)))
    cols = cols[ci]
  }else{
    ncol = length(q_o)
    if(hsv_grayscale){
      gray_vals = hsv_hue_min + (hsv_hue_max - hsv_hue_min)*((seq_len(ncol) - 1)/(ncol - 1))
      cols = gray(gray_vals)
    }else{
      cols = rainbow(ncol, start = hsv_hue_min, end = hsv_hue_max, s = hsv_saturation, v = hsv_value)
    }

    if(hsv_reverse) cols = rev(cols)
  }

  names(cols) = q_o

  #copy factor level order if any
  if(is.factor(bw_dt[, get(by_)])){
    plot_dt$facet_group = factor(plot_dt$facet_group, levels = levels(bw_dt[, get(by_)]))
  }
  p = ggplot(plot_dt) + geom_ribbon(aes_(x = ~V1, ymin = ~low, ymax = ~high, fill = ~q_range)) +
    labs(fill = "quantile band",
         y = y_, x = "bp",
         title = "Enrichment Profiles",
         subtitle = "aggregated by quantile range") +
    scale_fill_manual(values = cols) +
    scale_color_manual(breaks = q2do, palette = cols)
  if(length(todo) > 1){
    p = p +
      facet_grid(facet_group ~ ., scales = "free_y") +
      theme(strip.text.y = element_text(size = 7))
  }
  p
}

#' maps signal from 2 sample profiles to the x and y axis. axes are standard or "volcano" min XY vs fold-change Y/X
#'
#' @param bw_dt data.table of sample profiles
#' @param x_name sample name to map to x-axis, must be stored in variable specified in \code{xy_variable}
#' @param y_name sample name to map to y-axis, must be stored in variable specified in \code{xy_variable}
#' @param plotting_group work in progress, data.table that specifies groups for color scale
#' @param value_variable variable name that stores numeric values for plotting, default is "FE"
#' @param xy_variable variable name that stores sample, must contain entires for \code{x_name} and \code{y_name}
#' @param value_function a function to apply to \code{value_variable} in all combintations of \code{by_} per \code{x_name} and \code{y_name}
#' @param by_ variables that store individual measurement ids
#' @param plot_type standard or volcano, default is "standard"
#' @param show_help if TRUE overlay labels to aid plot interpretation, default is FALSE
#' @param fixed_coords if TRUE coordinate system is 1:1 ratio, default is TRUE
#'
#' @return ggplot of points comparing signal from 2 samples
regionSetPlotScatter = function(bw_dt, x_name, y_name,
                                plotting_group = NULL,
                                value_variable = "FE", xy_variable = "sample",
                                value_function = max,
                                by_ = "id",
                                plot_type = c("standard", "volcano")[1],
                                show_help = F, fixed_coords = T){
  xval = yval = xvolcano = id = yvolcano = NULL #declare binding for data.table
  plot_dt = merge(bw_dt[get(xy_variable) == x_name, .(xval = value_function(get(value_variable))), by = by_],
                  bw_dt[get(xy_variable) == y_name, .(yval = value_function(get(value_variable))), by = by_])
  if(is.null(plotting_group)){
    plot_dt$plotting_group = factor("none")
  }else{
    plot_dt = merge(plot_dt, plotting_group)
  }
  if(plot_type == "standard"){
    if(fixed_coords){
      xlim = c(0, plot_dt[, max(xval, yval)])
      ylim = xlim
    }else{
      xlim = c(0, plot_dt[, max(xval)])
      ylim = c(0, plot_dt[, max(yval)])
    }

    p = ggplot(plot_dt) +
      geom_point(aes(x = xval, y = yval, col = plotting_group)) +
      labs(x = x_name, y = y_name, title = "Max FE in regions") +
      ylim(ylim) + xlim(xlim)

    if(show_help){
      pos1 = .2 * max(lim)
      pos2 = .8 * max(lim)
      p = p + annotate("segment", x = 0, xend = 0, y = pos1, yend = pos2, arrow = arrow())
      p = p + annotate("label", x = 0, y = mean(lim), label = gsub(" ", "\n", paste(y_name, "binding")), hjust = 0)
      p = p + annotate("segment", x = pos1, xend = pos2, y = 0, yend = 0, arrow = arrow())
      p = p + annotate("label", x = mean(lim), y = 0, label = paste(x_name, "binding"), vjust = 0)
      p = p + annotate("label", x = 0, y = 0, label = "no\nbinding", hjust = 0, vjust = 0)
      p = p + annotate("segment", x = pos1, xend = pos2, y = pos1, yend = pos2, arrow = arrow())
      p = p + annotate("label", x = mean(lim), y = mean(lim), label = gsub(" ", "\n", paste("both binding")))
    }
  }else if(plot_type == "volcano"){
    plot_dt[, xvolcano := log2(max(yval, 1) / max(xval, 1)), by = id]
    plot_dt[, yvolcano := log2(max(min(yval, xval), 1)), by = id]
    xmax = plot_dt[, max(abs(c(xvolcano)))]
    lim = c(-xmax, xmax)
    p = ggplot(plot_dt) +
      geom_point(aes(x = xvolcano, y = yvolcano, col = plotting_group)) +
      labs(x = paste("log2 fold-change of", y_name, "over", x_name),
           y = paste("log2 min of", y_name, "and", x_name), title = "Max FE in regions") +
      xlim(lim)
    if(show_help){
      pos1 = .2 * max(lim)
      pos2 = .8 * max(lim)
      p = p + annotate("segment", y = 1, yend = 1, x = pos1, xend = pos2, arrow = arrow())
      p = p + annotate("label", y = 1, x = max(lim)/2, label = gsub(" ", "\n", paste(y_name, "binding")), vjust = 0)
      p = p + annotate("segment", y = 1, yend = 1, x = -pos1, xend = -pos2, arrow = arrow())
      p = p + annotate("label", y = 1, x = -max(lim)/2, label = gsub(" ", "\n", paste(x_name, "binding")), vjust = 0)
      p = p + annotate("label", x = 0, y = 1, label = "no\nbinding", hjust = .5, vjust = 0)
      ylim = range(plot_dt$yvolcano)
      ypos1 = 1 + (max(ylim) - min(ylim)) * .2
      ypos2 = 1 + (max(ylim) - min(ylim)) * .8
      p = p + annotate("segment", x = 0, xend = 0, y = ypos1, yend = ypos2, arrow = arrow())
      p = p + annotate("label", x = mean(lim), y = mean(ylim), label = gsub(" ", "\n", paste("both binding")))
      p
    }
  }
  if(is.null(plotting_group)){
    p = p + guides(color = "none")
  }
  if(fixed_coords) p = p + coord_fixed()
  p
}

#' heatmap
#'
#' @param bw_dt data.table of signals
#' @param max_disp numeric to sample number of rows down to, default is 2000
#'
#' @return ggplot heatmap of signal profiles, facetted by sample
regionSetPlotHeatmap = function(bw_dt, max_disp = 2000){
  id = xbp = x = to_disp = FE = hit = val = y = y_gap = NULL#declare binding for data.table
  plot_dt = copy(bw_dt)
  #sample down to max_disp
  uniq = unique(plot_dt$id)
  hit_disp = sample(uniq, min(max_disp, length(uniq)))
  plot_dt = plot_dt[id %in% hit_disp]

  #this should be warning or error, not silent removal
  # k = !duplicated(plot_dt[, paste(hit, x, sample)])
  # plot_dt = plot_dt[k,]

  # plot_dt[FE > 20, FE := 20]
  ###convert bp to x and add space between samples
  plot_dt[, xbp := x]
  plot_dt$sample = factor(plot_dt$sample, levels = to_disp)
  plot_dt[, x := xbp / 50 + 1]
  col_size = max(plot_dt$x)
  plot_dt[, x := x + col_size * (as.numeric(sample) - 1)]
  ##add space between samples
  col_gap_size = round(max(plot_dt$x) * .02)
  plot_dt[, x := x + col_gap_size * floor((x-1) / col_size)]
  ##text location for labels in center of each column
  x_txt = col_size / 2 + (col_size + col_gap_size) * (1:length(to_disp)-1)

  ###sort by row mean first
  hit_val = plot_dt[ , .(val = mean(FE)), by = hit]
  total_hit_o = as.character(hit_val[order(val), hit])
  setkey(plot_dt, hit)
  plot_dt = plot_dt[.(total_hit_o)]

  ###perform clustering on rows
  nclust = 6
  nclust = min(nclust, length(uniq)/2)
  if(nclust < 2) return(NULL)
  # print(paste("nclust", nclust))
  dc_dt = dcast(plot_dt, formula = hit ~ x, value.var = "FE")
  p_mat = as.matrix(dc_dt[,-1])
  p_mat[is.na(p_mat)] = 0
  rownames(p_mat) = dc_dt$hit
  km = kmeans(p_mat, centers = nclust)
  hc_centers = hclust(dist(km$centers, method = "euc"))
  hc_reo = 1:nclust
  names(hc_reo) = hc_centers$order

  ###debug stuff
  # layout(matrix(1))
  # plot(hc_centers)
  # heatmap.2(km$centers[names(hc_reo),], Colv = F, Rowv = F, trace = "n")
  # layout(matrix(1:10, ncol = 2))
  # par(mai = rep(0,4))
  # for(i in names(hc_reo)){
  #   test_hit = km$cluster == i
  #   test_hit = names(test_hit)[test_hit]
  #   plot(colMeans(p_mat[test_hit,]))
  #   title(i)
  # }

  # require(ggdendro)
  # dendro_data(hc_centers)
  o = order(hc_reo[as.character(km$cluster)])


  # plot(clusts)
  hit_o = names(km$cluster)[o]
  plot_dt$hit = factor(plot_dt$hit, levels = hit_o)
  plot_dt$y = as.numeric(factor(plot_dt$hit))

  ymax = max(plot_dt$y)
  row_gap_size = round(ymax * .02)

  clusts = hc_reo[as.character(km$cluster[o])]
  clust_counts = table(clusts)
  cs_counts = cumsum(clust_counts)
  starts = c(1, cs_counts[-nclust] + 1)
  ends = c(cs_counts)
  plot_dt$y_gap = 0
  setkey(plot_dt, y)
  for(i in 2:nclust){
    s = starts[i]
    e = ends[i]
    gap = row_gap_size * (i - 1)
    plot_dt[.(s:e), y_gap := gap ]
  }
  plot_dt[, y := y + y_gap]
  plot_dt$y_gap = NULL
  # setkey(plot_dt, hit)
  ggplot(plot_dt) +
    geom_raster(aes(x = x, y = y, fill = FE)) +
    # facet_grid(. ~ sample) +
    scale_y_reverse() + labs(x = "", y = "") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank()) +
    annotate("text", x = x_txt, y = - .02 * ymax, label = to_disp, hjust = .5, vjust = 0)


}



#' perform kmeans clustering on matrix rows and return reordered matrix along with order matched cluster assignments
#' clusters are sorted using hclust on centers
#'
#' @param mat numeric matrix to clutser
#' @param nclust the number of clusters
#' @param seed passed to set.seed() to allow reproducibility
clusteringKmeans = function(mat, nclust, seed = 0) {
  cluster_ordered = mat_name = NULL#declare binding for data.table

  set.seed(seed)
  mat_kmclust = kmeans(mat, centers = nclust, iter.max = 30)
  center_o = hclust(dist(mat_kmclust$centers))$order
  center_reo = 1:length(center_o)
  names(center_reo) = center_o
  center_reo[as.character(mat_kmclust$cluster)]
  mat_dt = data.table(mat_name = names(mat_kmclust$cluster), cluster = mat_kmclust$cluster, cluster_ordered = center_reo[as.character(mat_kmclust$cluster)])
  mat_dt = mat_dt[order(cluster_ordered), list(id = mat_name, group = cluster_ordered)]
  return(mat_dt)
}


#' perform kmeans clustering on matrix rows and return reordered matrix along with order matched cluster assignments
#' clusters are sorted using hclust on centers
#' the contents of each cluster are sorted using hclust
#' @param mat A wide format matrix
#' @param nclust the number of clusters
clusteringKmeansNestedHclust = function(mat, nclust) {
  group = id = within_o = NULL#declare binding for data.table
  mat_dt = clusteringKmeans(mat, nclust)
  mat_dt$within_o = as.integer(-1)
  for (i in 1:length(nclust)) {
    cmat = mat[mat_dt[group == i, id], , drop = F]
    if (nrow(cmat) > 2) {
      mat_dt[group == i, ]$within_o = hclust(dist((cmat)))$order
    } else {
      mat_dt[group == i, ]$within_o = 1:nrow(cmat)
    }

  }
  mat_dt = mat_dt[order(within_o), ][order(group), ]
  mat_dt$within_o = NULL
  return(mat_dt)
}

# clusteringKmeansNestedHclust.dt = function(dt, nclust, rownames_ = "id", colnames_ = "") {
#
# }
