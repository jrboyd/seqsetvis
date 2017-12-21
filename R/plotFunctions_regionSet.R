regionSetPlotBandedQuantiles = function(bw_dt, y_ = "FE", x_ = "x", by_ = "fake",
                                        hsv_reverse = F,
                                        hsv_saturation = 1, hsv_value = 1,
                                        hsv_grayscale = F,
                                        hsv_hue_min = 0, hsv_hue_max = 0.7, symm_colors = F,
                                        n_quantile = 18, quantile_min = 0.05, quantile_max = 0.95,
                                        is_centered = T, win_size = 50) {
  # hsv_min = 0; hsv_max = .7; n_quantile = 18; quantile_min = .05; quantile_max = .95; is_centered = T; win_size = 50
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
#
# regionSetPlotHeatmap
#
regionSetPlotScatter = function(bw_dt, x_name, y_name,
                                value_var = "FE", xy_var = "sample",
                                by_ = "id",

                                plot_val = c("max", "mean", "median"),
                                plot_type = c("standard", "volcano")[1],
                                show_help = F, fixed_coords = T){
  bw_dt = copy(all_bw_dt)
  x_name = "H7_H3K27ME3"
  y_name = "H7_H3K4ME3"
  bw_dt[get(xy_var) == x_name, max()]
  #isolate because these trigger updates of xy_plot_dt() which already triggers reactivity
  if(ptype == "standard"){
    bw_dt[, max(xval, yval)]
    lim = c(0, bw_dt[, max(xval, yval)])
    p = ggplot(bw_dt) +
      geom_point(aes(x = xval, y = yval, col = plotting_group)) +
      labs(x = x_, y = y_, title = "Max FE in regions") +
      ylim(lim) + xlim(lim)
    if(fixed_coords) p = p + coord_fixed()
    if(show_help){
      pos1 = .2 * max(lim)
      pos2 = .8 * max(lim)
      p = p + annotate("segment", x = 0, xend = 0, y = pos1, yend = pos2, arrow = arrow())
      p = p + annotate("label", x = 0, y = mean(lim), label = gsub(" ", "\n", paste(y_, "binding")), hjust = 0)
      p = p + annotate("segment", x = pos1, xend = pos2, y = 0, yend = 0, arrow = arrow())
      p = p + annotate("label", x = mean(lim), y = 0, label = paste(x_, "binding"), vjust = 0)
      p = p + annotate("label", x = 0, y = 0, label = "no\nbinding", hjust = 0, vjust = 0)
      p = p + annotate("segment", x = pos1, xend = pos2, y = pos1, yend = pos2, arrow = arrow())
      p = p + annotate("label", x = mean(lim), y = mean(lim), label = gsub(" ", "\n", paste("both binding")))
      p
    }
  }else if(ptype == "volcano"){
    bw_dt[, xvolcano := log2(max(yval, 1) / max(xval, 1)), by = id]
    bw_dt[, yvolcano := max(min(yval, xval), 1), by = id]
    xmax = bw_dt[, max(abs(c(xvolcano)))]
    lim = c(-xmax, xmax)
    p = ggplot(bw_dt) +
      geom_point(aes(x = xvolcano, y = yvolcano, col = plotting_group)) +
      labs(x = paste("log2 fold-change of", y_, "over", x_),
           y = paste("log2 min of", y_, "and", x_), title = "Max FE in regions") +
      xlim(lim)
    if(fixed_coords) p = p + coord_fixed()
    if(show_help){
      pos1 = .2 * max(lim)
      pos2 = .8 * max(lim)
      p = p + annotate("segment", y = 1, yend = 1, x = pos1, xend = pos2, arrow = arrow())
      p = p + annotate("label", y = 1, x = max(lim)/2, label = gsub(" ", "\n", paste(y_, "binding")), vjust = 0)
      p = p + annotate("segment", y = 1, yend = 1, x = -pos1, xend = -pos2, arrow = arrow())
      p = p + annotate("label", y = 1, x = -max(lim)/2, label = gsub(" ", "\n", paste(x_, "binding")), vjust = 0)
      p = p + annotate("label", x = 0, y = 1, label = "no\nbinding", hjust = .5, vjust = 0)
      ylim = range(bw_dt$yvolcano)
      ypos1 = 1 + (max(ylim) - min(ylim)) * .2
      ypos2 = 1 + (max(ylim) - min(ylim)) * .8
      p = p + annotate("segment", x = 0, xend = 0, y = ypos1, yend = ypos2, arrow = arrow())
      p = p + annotate("label", x = mean(lim), y = mean(ylim), label = gsub(" ", "\n", paste("both binding")))
      p
    }
  }
  p
}

#' perform kmeans clustering on matrix rows and return reordered matrix along with order matched cluster assignments
#' clusters are sorted using hclust on centers
clusteringKmeans = function(mat, nclust, seed = 0) {
  set.seed(seed)
  mat_kmclust = kmeans(mat, centers = nclust, iter.max = 30)
  center_o = hclust(dist(mat_kmclust$centers))$order
  center_reo = 1:length(center_o)
  names(center_reo) = center_o
  center_reo[as.character(mat_kmclust$cluster)]
  mat_dt = data.table(mat_name = names(mat_kmclust$cluster), cluster = mat_kmclust$cluster, cluster_ordered = center_reo[as.character(mat_kmclust$cluster)])
  mat_dt = mat_dt[order(cluster_ordered), .(id = mat_name, group = cluster_ordered)]
  return(mat_dt)
}


#' perform kmeans clustering on matrix rows and return reordered matrix along with order matched cluster assignments
#' clusters are sorted using hclust on centers
#' the contents of each cluster are sorted using hclust
#' @param mat A wide format matrix
#' @param nclust the number of clusters
clusteringKmeansNestedHclust = function(mat, nclust) {
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

clusteringKmeansNestedHclust.dt = function(dt, nclust, rownames_ = "id", colnames_ = "") {

}
