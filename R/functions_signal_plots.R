#' plot profiles from bigwigs
#' @export
#' @param bw_data a GRanges or data.table of bigwig signal.
#' As returned from \code{\link{fetchWindowedBigwig}} and \code{\link{fetchWindowedBigwigList}}
#' @param y_ the variable name in bw_data for y axis in plot
#' @param x_ the variable name in bw_data for x axis in plot
#' @param by_ the variable name in bw_data to facet on
#' @param hsv_reverse logical, should color scale be reversed? default FALSE
#' @param hsv_saturation numeric [0, 1] saturation for color scale. default 1
#' @param hsv_value numeric [0, 1] value for color scale. default 1
#' @param hsv_grayscale logical, if TRUE gray() is used instead of rainbow(). default FALSE
#' @param hsv_hue_min numeric [0, hsv_hue_max) hue min of color scale
#' @param hsv_hue_max numeric (hsv_hue_min, 1] hue max of color scale
#' @param hsv_symmetric if TRUE, colorscale is symmetrical, default FALSE.
#' @param n_quantile number of evenly size quantile bins
#' @param quantile_min the lowest quantile start
#' @param quantile_max the highest quantile end
#'
#' @return ggplot object using ribbon plots to show quantile distributions
#' @import ggplot2
#' @importFrom grDevices gray rainbow
#' @importFrom stats quantile
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @examples
#' #rainbow colors
#' ssvSignalBandedQuantiles(CTCF_in_10a_profiles_gr)
#' ssvSignalBandedQuantiles(CTCF_in_10a_profiles_gr, hsv_symmetric = TRUE,
#'     hsv_reverse = TRUE)
#' #grayscale
#' ssvSignalBandedQuantiles(CTCF_in_10a_profiles_gr, hsv_grayscale = TRUE)
#' ssvSignalBandedQuantiles(CTCF_in_10a_profiles_gr, hsv_grayscale = TRUE,
#'     hsv_symmetric = TRUE, hsv_reverse = TRUE)
#' #using "by_" per sample
#' ssvSignalBandedQuantiles(CTCF_in_10a_profiles_gr, hsv_grayscale = TRUE,
#'     hsv_symmetric = TRUE, hsv_reverse = TRUE, by_ = "sample")
#' #adding spline smoothing
#' splined = applySpline(CTCF_in_10a_profiles_gr, n = 10,
#'     by_ = c("id", "sample"))
#' ssvSignalBandedQuantiles(splined, n_quantile = 50,
#'     quantile_min = .25, quantile_max = .75,
#'     hsv_symmetric = TRUE, hsv_reverse = TRUE, by_ = "sample")
ssvSignalBandedQuantiles = function(bw_data, y_ = "y", x_ = "x", by_ = "fake",
                                    hsv_reverse = FALSE,
                                    hsv_saturation = 1, hsv_value = 1,
                                    hsv_grayscale = FALSE,
                                    hsv_hue_min = 0, hsv_hue_max = 0.7, hsv_symmetric = FALSE,
                                    n_quantile = 18, quantile_min = 0.05, quantile_max = 0.95
) {
    if(class(bw_data)[1] == "GRanges"){
        bw_data = data.table::as.data.table(bw_data)
    }
    stopifnot(data.table::is.data.table(bw_data))
    stopifnot(is.character(x_), is.character(y_), is.character(by_))
    stopifnot(x_ %in% colnames(bw_data), y_ %in% colnames(bw_data))
    stopifnot(by_ %in% colnames(bw_data) || by_ == "fake")
    stopifnot(is.logical(hsv_reverse), is.logical(hsv_grayscale),
              is.logical(hsv_symmetric))
    num_args = c(hsv_saturation, hsv_value, hsv_hue_min,
                 hsv_hue_max, quantile_min, quantile_max)
    stopifnot(all(is.numeric(num_args)))
    stopifnot(all(num_args >= 0), all(num_args <= 1))
    stopifnot(n_quantile > 2)

    variable = value = V1 = low = high = low_q = high_q = NULL
    q_range = rn = q_num = q_low = q_high = NULL #declare binding for data.table
    q2do = 0:n_quantile/n_quantile
    q2do = round(quantile_min + q2do * (quantile_max - quantile_min), digits = 3)

    if(by_ == "fake"){
        bw_data[[by_]] = TRUE
        todo = TRUE
    }else{
        todo = unique(bw_data[, get(by_)])
    }

    calc_quantiles = function(td){#, x_, y_, by_, q2do){
        dt = bw_data[get(by_) == td, list(qs = list(list(stats::quantile(get(y_), q2do)))), by = x_]
        qmat = matrix(unlist(dt$qs), nrow = nrow(dt), ncol = length(q2do), byrow = TRUE)
        colnames(qmat) = paste0(q2do*100, "%")
        # sapply(dt$qs, function(x)x[[1]])
        dt = cbind(dt, qmat)
        dt$qs = NULL
        dtm = data.table::melt(dt, id.vars = x_)
        data.table::setkey(dtm, variable)
        q2do_str = paste0(q2do * 100, "%")
        dt_low = dtm[q2do_str[-length(q2do_str)], list(get(x_), low_q = variable, low = value)]
        dt_high = dtm[q2do_str[-1], list(get(x_), high_q = variable, high = value)]
        dt_c = cbind(dt_low, dt_high[, -1])
        dt_c = dt_c[, list(V1, low, high, q_range = paste0(sub("%", "", low_q), "-", high_q))]
        q_o = unique(dt_c$q_range)
        dt_c$q_range = factor(dt_c$q_range, levels = q_o)

        dt_c$rn = seq_len(nrow(dt_c))
        dt_c[, `:=`(c("q_low", "q_high"), data.table::tstrsplit(sub("%", "", q_range), split = "-", keep = seq_len(2))), by = rn]
        dt_c[, `:=`(q_num, (as.numeric(q_low) + as.numeric(q_high))/2)]
        dt_c[, `:=`(c("rn", "q_low", "q_high"), NULL)]
        dt_c$facet_group = td
        return(dt_c)
    }

    all_q = lapply(todo, calc_quantiles)
    plot_dt = data.table::rbindlist(all_q)
    q_o = unique(plot_dt$q_range)
    if(hsv_symmetric){
        ncol = ceiling(length(q_o) / 2)
        if(hsv_grayscale){
            gray_vals = hsv_hue_min + (hsv_hue_max - hsv_hue_min)*((seq_len(ncol) - 1)/(ncol - 1))
            cols = grDevices::gray(gray_vals)
        }else{
            cols = grDevices::rainbow(ncol, start = hsv_hue_min, end = hsv_hue_max, s = hsv_saturation, v = hsv_value)
        }
        if(hsv_reverse) cols = rev(cols)
        ci = c(seq_len(ncol), rev(seq_len(ncol) - ifelse(length(q_o) %% 2 == 0, 0, 1)))
        cols = cols[ci]
    }else{
        ncol = length(q_o)
        if(hsv_grayscale){
            gray_vals = hsv_hue_min + (hsv_hue_max - hsv_hue_min)*((seq_len(ncol) - 1)/(ncol - 1))
            cols = grDevices::gray(gray_vals)
        }else{
            cols = grDevices::rainbow(ncol, start = hsv_hue_min, end = hsv_hue_max, s = hsv_saturation, v = hsv_value)
        }

        if(hsv_reverse) cols = rev(cols)
    }

    names(cols) = q_o

    #copy factor level order if any
    if(is.factor(bw_data[, get(by_)])){
        plot_dt$facet_group = factor(plot_dt$facet_group, levels = levels(bw_data[, get(by_)]))
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
            facet_grid(facet_group ~ .) +
            theme(strip.text.y = element_text(size = 7))
    }
    p
}

#' maps signal from 2 sample profiles to the x and y axis. axes are standard or "volcano" min XY vs fold-change Y/X
#' @export
#' @param bw_data a GRanges or data.table of bigwig signal.
#' As returned from \code{\link{fetchWindowedBigwig}} and \code{\link{fetchWindowedBigwigList}}
#' @param x_name sample name to map to x-axis, must be stored in variable specified in \code{xy_variable}
#' @param y_name sample name to map to y-axis, must be stored in variable specified in \code{xy_variable}
#' @param color_table data.frame with 2 columns, one of which must be named "group" and gets mapped to color.
#' The other column must be the same as by_ parameter and is used for merging.
#' @param value_variable variable name that stores numeric values for plotting, default is "y"
#' @param xy_variable variable name that stores sample, must contain entires for \code{x_name} and \code{y_name}
#' @param value_function a function to apply to \code{value_variable} in all combintations of \code{by_} per \code{x_name} and \code{y_name}
#' @param by_ variables that store individual measurement ids
#' @param plot_type standard or volcano, default is "standard"
#' @param show_help if TRUE overlay labels to aid plot interpretation, default is FALSE
#' @param fixed_coords if TRUE coordinate system is 1:1 ratio, default is TRUE
#'
#' @import ggplot2
#' @return ggplot of points comparing signal from 2 samples
#' @examples
#' ssvSignalScatterplot(CTCF_in_10a_profiles_gr,
#'     x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF")
#' ssvSignalScatterplot(CTCF_in_10a_profiles_gr,
#'     x_name = "MCF10A_CTCF", y_name = "MCF10CA1_CTCF")
#'
#' ssvSignalScatterplot(CTCF_in_10a_profiles_gr,
#'     x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF",
#'     value_function = median) + labs(title = "median FE in regions")
#'
#' ssvSignalScatterplot(CTCF_in_10a_profiles_gr,
#'     x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF",
#'     plot_type = "volcano")
#'
#' ssvSignalScatterplot(CTCF_in_10a_profiles_gr,
#'     x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF",
#'     plot_type = "volcano", show_help = TRUE)
ssvSignalScatterplot = function(bw_data, x_name, y_name,
                                color_table = NULL,
                                value_variable = "y", xy_variable = "sample",
                                value_function = max,
                                by_ = "id",
                                plot_type = c("standard", "volcano")[1],
                                show_help = FALSE, fixed_coords = TRUE){
    xval = yval = xvolcano = id = yvolcano = group = NULL #declare binding for data.table
    if(class(bw_data)[1] == "GRanges"){
        bw_data = data.table::as.data.table(bw_data)
    }
    stopifnot(data.table::is.data.table(bw_data))
    stopifnot(xy_variable %in% colnames(bw_data))
    if(!any(x_name == bw_data[[xy_variable]])){
        stop(x_name, "not found in", xy_variable, "variable of data.table")
    }
    if(!any(y_name == bw_data[[xy_variable]])){
        stop(y_name, "not found in", xy_variable, "variable of data.table")
    }
    stopifnot(is.character(x_name), is.character(y_name),
              is.character(xy_variable), is.character(by_),
              is.character(plot_type))
    stopifnot(plot_type %in% c("standard", "volcano"))
    stopifnot(is.logical(show_help), is.logical(fixed_coords))
    stopifnot(is.function(value_function))

    plot_dt = merge(bw_data[get(xy_variable) == x_name, list(xval = value_function(get(value_variable))), by = by_],
                    bw_data[get(xy_variable) == y_name, list(yval = value_function(get(value_variable))), by = by_])
    if(!is.null(color_table)){
        plot_dt = merge(plot_dt, color_table, by = by_)
    }
    # if(is.null(plotting_group)){
    #     plot_dt$plotting_group = factor("none")
    # }else{
    #     plot_dt = merge(plot_dt, plotting_group)
    # }
    if(plot_type == "standard"){
        if(fixed_coords){
            xlim = c(0, plot_dt[, max(xval, yval)])
            ylim = xlim
        }else{
            xlim = c(0, plot_dt[, max(xval)])
            ylim = c(0, plot_dt[, max(yval)])
        }

        p = ggplot(plot_dt, aes(x = xval, y = yval))
        if(is.null(color_table)){
            p = p + geom_point(mapping = aes(alpha = 1, size = 1.5))
        }else{
            p = p + geom_point(mapping = aes(alpha = 1, size = 1.5, color = group))
        }
        p = p +
            scale_alpha_identity() +
            scale_size_identity() +
            guides(alpha = "none", size = "none") +
            labs(x = x_name, y = y_name, title = paste(value_variable, "aggregated per", paste(by_, collapse = ", "))) +
            ylim(ylim) + xlim(xlim)

        if(show_help){
            xpos1 = .2 * max(xlim)
            xpos2 = .8 * max(xlim)
            ypos1 = .2 * max(ylim)
            ypos2 = .8 * max(ylim)
            p = p + annotate("segment", x = 0, xend = 0, y = ypos1, yend = ypos2, arrow = arrow())
            p = p + annotate("label", x = 0, y = mean(xlim), label = gsub(" ", "\n", paste(y_name, "binding")), hjust = 0)
            p = p + annotate("segment", x = xpos1, xend = xpos2, y = 0, yend = 0, arrow = arrow())
            p = p + annotate("label", x = mean(xlim), y = 0, label = paste(x_name, "binding"), vjust = 0)
            p = p + annotate("label", x = 0, y = 0, label = "no\nbinding", hjust = 0, vjust = 0)
            p = p + annotate("segment", x = xpos1, xend = xpos2, y = ypos1, yend = ypos2, arrow = arrow())
            p = p + annotate("label", x = mean(xlim), y = mean(ylim), label = gsub(" ", "\n", paste("both binding")))
        }
    }else if(plot_type == "volcano"){
        plot_dt[, xvolcano := log2(max(yval, 1) / max(xval, 1)), by = id]
        plot_dt[, yvolcano := log2(max(min(yval, xval), 1)), by = id]
        xmax = plot_dt[, max(abs(c(xvolcano)))]
        lim = c(-xmax, xmax)
        p = ggplot(plot_dt, aes(x = xvolcano, y = yvolcano))
        if(is.null(color_table)){
            p = p + geom_point(mapping = aes(alpha = 1, size = 1.5))
        }else{
            p = p + geom_point(mapping = aes(alpha = 1, size = 1.5, color = group))
        }
        p = p +
            scale_size_identity() +
            scale_alpha_identity() +
            guides(alpha = "none", size = "none") +
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
    # if(is.null(plotting_group)){
    #     p = p + guides(color = "none")
    # }
    if(fixed_coords) p = p + coord_fixed()
    p
}

#' clustering as for a heatmap
#' @export
#' @param bw_data a GRanges or data.table of bigwig signal.
#' As returned from \code{\link{fetchWindowedBigwig}} and \code{\link{fetchWindowedBigwigList}}
#' @param nclust number of clusters
#' @param row_ variable name mapped to row, likely peak id or gene name for ngs data
#' @param column_ varaible mapped to column, likely bp position for ngs data
#' @param cluster_ variable name to use for cluster info
#' @param fill_ numeric variable to map to fill
#' @param facet_ variable name to facet horizontally by
#' @param max_rows for speed rows are sampled to 500 by default, use Inf to plot full data
#' @param max_cols for speed columns are sampled to 100 by default, use Inf to plot full data
#' @param clustering_col_min numeric minimum for col range considered when clustering, default in -Inf
#' @param clustering_col_max numeric maximum for col range considered when clustering, default in Inf
#' @rawNamespace import(data.table, except = c(shift, first, second))
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
ssvSignalClustering = function(bw_data, nclust = 6,
                               row_ = "id",
                               column_ = "x", fill_ = "y", facet_ = "sample",
                               cluster_ = "cluster_id",
                               max_rows = 500, max_cols = 100,
                               clustering_col_min = -Inf, clustering_col_max = Inf){
    id = xbp = x = to_disp = y = hit = val = y = y_gap = group =  NULL#declare binding for data.table
    output_GRanges = FALSE
    if(class(bw_data)[1] == "GRanges"){
        bw_data = data.table::as.data.table(bw_data)
        output_GRanges = TRUE
    }
    stopifnot(is.data.table(bw_data))
    stopifnot(is.numeric(nclust) || nclust < 2)
    stopifnot(is.character(row_), is.character(column_), is.character(fill_),
              is.character(facet_), is.character(cluster_))
    stopifnot(row_ %in% colnames(bw_data), column_ %in% colnames(bw_data),
              fill_ %in% colnames(bw_data))
    stopifnot(facet_ %in% colnames(bw_data) || facet_ == "")
    stopifnot(is.numeric(max_rows), is.numeric(max_cols),
              is.numeric(clustering_col_min), is.numeric(clustering_col_max))

    plot_dt = data.table::copy(bw_data)
    raw_nc = length(unique(plot_dt[[column_]]))
    if(raw_nc > max_cols){
        new_scale = (seq_len(max_cols)-1) / (max_cols - 1)
        old_scale = (seq_len(raw_nc)-1) / (raw_nc - 1)
        kept = vapply(new_scale, function(v){
            which.min(abs(v - old_scale))
        }, 1)
        kept = sort(unique(plot_dt[[column_]]))[kept]
        plot_dt = plot_dt[x %in% kept]
        warning(raw_nc - max_cols,
                " columns were discarded according to max_cols: ",
                max_cols)
    }
    row_ids = unique(plot_dt[[row_]])
    raw_nr = length(row_ids)
    if(raw_nr > max_rows){
        set.seed(0)
        row_ids = sample(row_ids, max_rows)
        plot_dt = plot_dt[get(row_) %in% row_ids]
        warning(raw_nr - max_rows,
                " rows were discarded according to max_cols: ",
                max_rows)
    }

    dc_formula = paste(row_, "~", paste(c(facet_, column_), collapse = " + "))
    if(is.numeric(plot_dt[[column_]])){
        dc_dt = data.table::dcast(plot_dt[get(column_) > clustering_col_min &
                                              get(column_) < clustering_col_max],
                                  formula = dc_formula,
                                  value.var = fill_)
    }else{
        dc_dt = data.table::dcast(plot_dt,
                                  formula = dc_formula,
                                  value.var = fill_)
    }

    dc_mat = as.matrix(dc_dt[,-1])
    rownames(dc_mat) = dc_dt[[row_]]
    rclusters = clusteringKmeansNestedHclust(dc_mat, nclust = nclust)
    rclusters = rclusters[rev(seq_len(nrow(rclusters))),]

    plot_dt[[row_]] = factor(plot_dt[[row_]], levels = rclusters[[row_]])
    data.table::setkey(rclusters, id)
    plot_dt[[cluster_]] = rclusters[list(plot_dt$id), group]
    if(output_GRanges){
        plot_dt = GRanges(plot_dt)
    }
    return(plot_dt)
}

#' heatmap style representation of membership table.
#' instead of clustering, each column is sorted starting from the left.
#' @export
#' @param bw_data a GRanges or data.table of bigwig signal.
#' As returned from \code{\link{fetchWindowedBigwig}} and \code{\link{fetchWindowedBigwigList}}
#' @param nclust number of clusters
#' @param perform_clustering should clustering be done? default is auto.
#' auto considers if row_ has been ordered by being a factor and if cluster_ is a numeric.
#' @param row_ variable name mapped to row, likely peak id or gene name for ngs data
#' @param column_ varaible mapped to column, likely bp position for ngs data
#' @param fill_ numeric variable to map to fill
#' @param facet_ variable name to facet horizontally by
#' @param cluster_ variable name to use for cluster info
#' @param max_rows for speed rows are sampled to 500 by default, use Inf to plot full data
#' @param max_cols for speed columns are sampled to 100 by default, use Inf to plot full data
#' @param clustering_col_min numeric minimum for col range considered when clustering, default in -Inf
#' @param clustering_col_max numeric maximum for col range considered when clustering, default in Inf
#'
#' @import ggplot2
#' @return ggplot heatmap of signal profiles, facetted by sample
#' @examples
#' #the simplest use
#' ssvSignalHeatmap(CTCF_in_10a_profiles_gr)
#'
#' #clustering can be done manually beforehand
#' clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_gr, nclust = 3)
#' ssvSignalHeatmap(clust_dt)
ssvSignalHeatmap = function(bw_data,
                            nclust = 6,
                            perform_clustering = c("auto", "yes", "no")[1],
                            row_ = "id",
                            column_ = "x",
                            fill_ = "y",
                            facet_ = "sample",
                            cluster_ = "cluster_id",
                            max_rows = 500,
                            max_cols = 100,
                            clustering_col_min = -Inf,
                            clustering_col_max = Inf){
    id = xbp = x = to_disp = y = hit = val = y = y_gap = cluster_id = NULL#declare binding for data.table
    if(class(bw_data)[1] == "GRanges"){
        bw_data = data.table::as.data.table(bw_data)
    }
    stopifnot(is.data.table(bw_data))
    stopifnot(is.numeric(nclust) || nclust < 2)
    stopifnot(is.character(row_), is.character(column_), is.character(fill_),
              is.character(facet_), is.character(cluster_))
    stopifnot(row_ %in% colnames(bw_data), column_ %in% colnames(bw_data),
              fill_ %in% colnames(bw_data))
    stopifnot(facet_ %in% colnames(bw_data) || facet_ == "")
    stopifnot(is.numeric(max_rows), is.numeric(max_cols),
              is.numeric(clustering_col_min), is.numeric(clustering_col_max))
    #determine if user wants clustering
    do_cluster = perform_clustering == "yes"
    if(perform_clustering == "auto"){
        if(is.factor(bw_data[[row_]]) & is.numeric(bw_data[[cluster_]])){
            do_cluster = FALSE
        }else{
            do_cluster = TRUE
        }
    }
    if(do_cluster){
        message("clustering...")
        plot_dt = ssvSignalClustering(bw_data = bw_data,
                                      nclust = nclust,
                                      row_ = row_,
                                      column_ = column_,
                                      fill_ = fill_,
                                      facet_ = facet_,
                                      cluster_ = cluster_,
                                      max_rows = max_rows,
                                      max_cols = max_cols,
                                      clustering_col_min = clustering_col_min,
                                      clustering_col_max = clustering_col_max)
    }else{
        plot_dt = bw_data
    }
    message("making plot...")
    scale_floor = .1
    scale_vals = c(0, scale_floor + 0:10/10*(1 - scale_floor))
    p = ggplot(plot_dt) +
        geom_raster(aes_string(x = column_, y = row_, fill = fill_)) +
        facet_grid(paste(". ~", facet_)) +
        theme(axis.line = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.background = element_blank()) +
        scale_fill_distiller(type = "div", palette = "Spectral",
                             values = scale_vals)

    xs = sort(unique(plot_dt$x), decreasing = FALSE)
    xleft = min(xs) - (xs[2] - xs[1])/2
    xright = max(xs) + (xs[2] - xs[1])/2
    if(xleft < 0 & xright > 0){
        xbr = c(xleft, 0, xright)
    }else{
        xbr = c(xleft, xright)
    }

    p = p + scale_x_continuous(breaks = xbr) +
        theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = .5))

    rclust = plot_dt[, list(cluster_id = unique(get(cluster_))), by = get(row_)]
    ends = cumsum(rev(table(rclust$cluster_id)))
    starts = c(1, ends[-length(ends)] + 1)
    starts = starts - .5
    ends = ends + .5

    xfactor = diff(range(plot_dt$x))
    xfloor = min(plot_dt$x)

    df_rects = data.frame(xmin = xfloor - .12*xfactor,
                          xmax = xfloor - .03*xfactor,
                          ymin = starts,
                          ymax = ends)
    df_rects$fill = c("black", "gray")[seq_len(nrow(df_rects))%%2+1]
    df_rects$color = c("gray", "black")[seq_len(nrow(df_rects))%%2+1]
    df_rects = df_rects[rev(seq_len(nrow(df_rects))),]
    for(i in seq_len(nrow(df_rects))){
        p = p + annotate("rect",
                         xmin = df_rects$xmin[i],
                         xmax = df_rects$xmax[i],
                         ymin= df_rects$ymin[i],
                         ymax = df_rects$ymax[i],
                         fill = df_rects$fill[i])
        p = p + annotate("text",
                         x = mean(c(df_rects$xmin[i], df_rects$xmax[i])),
                         y = mean(c(df_rects$ymin[i], df_rects$ymax[i])),
                         label = i,
                         color = df_rects$color[i])
    }
    p
}

#' construct line type plots where each region in each sample is represented
#' @export
#' @param bw_data a GRanges or data.table of bigwig signal.
#' As returned from \code{\link{fetchWindowedBigwig}} and \code{\link{fetchWindowedBigwigList}}
#' @param x_ variable name mapped to x aesthetic, x by default.
#' @param y_ variable name mapped to y aesthetic, y by default.
#' @param color_ variable name mapped to color aesthetic, sample by default.
#' @param sample_ variable name, along with region_ used to group and
#' facet by default, change group_ or facet_ to override.
#' @param region_ variable name, along with sample_ used to group and facet
#' by default, change group_ or facet_ to override.
#' @param group_ group aesthetic keeps lines of geom_path from mis-connecting.
#' auto_grp by default which combines sample_ and region_.
#' probably shouldn't change.
#' @param line_alpha alpha value for lines. default is 1.
#' @param facet_ facetting divides up plots.
#' auto_facet by default which combines sample_ and region_.
#' if overriding facet_method with facet_grid, make sure to include ~ between
#' two variables, ie. "a~b", ".~b", "a~."
#' @param facet_method ggplot2 facetting method or wrapper for same,
#' facet_wrap by default.
#' @param spline_n if not NULL, applySpline will be called with n = spline_n.
#' default is NULL.
#' @import ggplot2
#' @return ggplot of signal potentially facetted by region and sample
#' @examples
#' bw_gr = CTCF_in_10a_profiles_gr
#' ssvSignalLineplot(subset(bw_gr, bw_gr$id %in% seq_len(3)), facet_ = "sample")
#' ssvSignalLineplot(subset(bw_gr, bw_gr$id %in% seq_len(3)),
#'     facet_ = "sample~.",
#'     facet_method = facet_grid)
#' ssvSignalLineplot(subset(bw_gr, bw_gr$id %in% seq_len(3)),
#'     facet_ = paste("sample", "~", "id"), facet_method = facet_grid)
#' ssvSignalLineplot(subset(bw_gr, bw_gr$id %in% seq_len(3)))
#' ssvSignalLineplot(subset(bw_gr, bw_gr$id %in% seq_len(3)), facet_ = "id")
#' ssvSignalLineplot(subset(bw_gr, bw_gr$id %in% seq_len(3)),
#'     facet_ = "id", spline_n = 10)
ssvSignalLineplot = function(bw_data, x_ = "x", y_ = "y", color_ = "sample",
                             sample_ = "sample", region_ = "id",
                             group_ = "auto_grp", line_alpha = 1,
                             facet_ = "auto_facet",
                             facet_method = facet_wrap, spline_n = NULL){
    auto_grp = auto_facet = NULL
    if(class(bw_data)[1] == "GRanges"){
        bw_data = data.table::as.data.table(bw_data)
    }
    stopifnot(is.data.table(bw_data))
    stopifnot(is.character(x_), is.character(y_), is.character(color_),
              is.character(sample_), is.character(region_),
              is.character(group_), is.character(facet_))
    stopifnot(x_ %in% colnames(bw_data), y_ %in% colnames(bw_data),
              color_ %in% colnames(bw_data), sample_ %in% colnames(bw_data),
              region_ %in% colnames(bw_data))
    stopifnot(group_ %in% colnames(bw_data) || group_ == "auto_grp")
    facet_names = strsplit(facet_, " ?[~+] ?")[[1]]
    facet_names = facet_names[facet_names != "."]
    stopifnot(all(facet_names %in% colnames(bw_data)) || facet_ == "auto_facet")

    stopifnot(is.function(facet_method))
    stopifnot(is.numeric(spline_n) || is.null(spline_n))
    bw_data[,auto_grp := paste(get(sample_), get(region_))]
    bw_data[,auto_facet := paste(get(sample_), get(region_))]
    if(!is.null(spline_n)){
        p_dt = applySpline(bw_data, n = spline_n, x_ = x_, y_ = y_,
                           by_ = c(group_, sample_, region_, facet_))
    }else{
        p_dt = bw_data
    }
    ggplot(p_dt) + geom_path(aes_string(x = x_,
                                        y = y_,
                                        col = color_,
                                        group = group_),
                             alpha = line_alpha) +
        facet_method(facet_)
}

#' aggregate line signals in a single line plot
#' @export
#' @param bw_data a GRanges or data.table of bigwig signal.
#' As returned from \code{\link{fetchWindowedBigwig}} and \code{\link{fetchWindowedBigwigList}}
#' @param x_ variable name mapped to x aesthetic, x by default.
#' @param y_ variable name mapped to y aesthetic, y by default.
#' @param sample_ variable name, along with region_ used to group by default,
#' @param color_ variable name mapped to color aesthetic, sample_ by default.
#' change group_ to override.
#' @param group_ group aesthetic keeps lines of geom_path from mis-connecting.
#' auto_grp by default uses sample_. probably shouldn't change.
#' @param agg_fun the aggregation function to apply by sample_ and x_,
#' default is mean
#' @param spline_n if not NULL, applySpline will be called with n = spline_n.
#' default is NULL.
#' @import ggplot2
#' @return ggplot of signal aggregated with agg_fun() by sample.
#' @examples
#' bw_gr = CTCF_in_10a_profiles_gr
#' ssvSignalLineplotAgg(bw_gr) +
#'     labs(title = "agg regions by sample.")
#' ssvSignalLineplotAgg(CTCF_in_10a_profiles_gr, spline_n = 10) +
#'     labs(title = "agg regions by sample, with spline smoothing.")
#' ssvSignalLineplotAgg(subset(bw_gr, bw_gr$id %in% seq_len(10)),
#'     sample_ = "id", color_ = "id") +
#'     labs(title = "agg samples by region id (weird)")
#' ssvSignalLineplotAgg(subset(bw_gr, bw_gr$id %in% seq_len(10)), sample_ = "id",
#'     color_ = "id", spline_n = 10) +
#'     labs(title = "agg samples by region id (weird), with spline smoothing")
ssvSignalLineplotAgg = function(bw_data, x_ = "x", y_ = "y",
                                sample_ = "sample",
                                color_ = sample_,
                                group_ = "auto_grp", agg_fun = mean,
                                spline_n = NULL){
    auto_grp = NULL
    if(class(bw_data)[1] == "GRanges"){
        bw_data = data.table::as.data.table(bw_data)
    }
    stopifnot(is.data.table(bw_data))
    stopifnot(is.character(x_), is.character(y_),
              is.character(sample_), is.character(color_),
              is.character(group_))
    stopifnot(x_ %in% colnames(bw_data), y_ %in% colnames(bw_data),
              color_ %in% colnames(bw_data), sample_ %in% colnames(bw_data))
    stopifnot(group_ %in% colnames(bw_data) || group_ == "auto_grp")

    stopifnot(is.function(agg_fun))
    stopifnot(is.numeric(spline_n) || is.null(spline_n))

    bw_data[, auto_grp := paste(get(sample_))]
    agg_dt = bw_data[, list(y = agg_fun(get(y_))),
                     by = c(unique(c(group_, color_, sample_, x_)))]

    if(!is.null(spline_n)){
        p_dt = applySpline(agg_dt, n = spline_n, x_ = x_, y_ = y_,
                           by_ = c(group_, sample_))
    }else{
        p_dt = agg_dt
    }
    p_dt = p_dt[order(get(x_))]
    ggplot(p_dt) +
        geom_path(aes_string(x = x_, y = y_, col = color_, group = group_))
}


