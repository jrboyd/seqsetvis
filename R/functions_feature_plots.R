# plotting functions for region/feature set comparison go here
# all plotting functions accept a variety of objects as first argument and
# handle them via the S4  generic ssvMakeMembTable
# ssvMakeMembTable has implemented methods for:
#   - membership table
#       columns are sets
#       rows are members
#       values are logicals indicating membership
#   - list containing character vectors of set members
#   - list of GRanges

#' ssvFeatureVenn
#'
#' ggplot implementation of vennDiagram from limma package.  Currently limited
#' at 3 sets.  [ssvFeatureUpset] and [ssvFeatureBinaryHeatmap] are good options for
#' more than 3 sets. [ssvFeatureEuler] can work too but can take a very long time
#' to run for more than 5 or so.
#'
#' @export
#'
#' @param object will be passed to [ssvMakeMembTable] for
#' conversion to membership matrix
#' @param group_names useful if names weren't provided or were lost in
#' creating membership matrix
#' @param counts_txt_size font size for count numbers
#' @param counts_as_labels if TRUE, geom_label is used instead of geom_text.
#' can be easier to read.
#' @param show_outside_count if TRUE, items outside of all sets are counted
#' outside. can be confusing.
#' @param line_width uses size aesthetic to control line width of circles.
#' @param circle_colors colors to use for circle line colors. Uses Dark2 set
#' from RColorBrewer by default.
#' @param fill_alpha alpha value to use for fill, defaults to .3.
#' @param line_alpha numeric value from 0 to 1. Alpha value for circle line
#' @param counts_color character. single color to use for displaying counts
#' @param counts_as_percent if TRUE, convert counts to percentages in plots.
#' @param percentage_digits The number of digits to round percentages to, default is 1.
#' @param percentage_suffix The character to append to percentages, default is "%".
#' @param n_points integer.  number of points to approximate circle with.
#' default is 200.
#' @param return_data logical.  If TRUE, return value is no longer ggplot and
#' is instead the data used to generate that plot. Default is FALSE.
#' @return ggplot venn diagram
#' @import ggplot2
#' @import S4Vectors
#' @importFrom limma vennCounts
#' @examples
#' ssvFeatureVenn(list(1:3, 2:6))
#' ssvFeatureVenn(CTCF_in_10a_overlaps_gr)
#' ssvFeatureVenn(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
#'
#' ssvFeatureVenn(list(1:3, 2:6),
#'   counts_as_percent = TRUE,
#'   percentage_digits = 2)
#'
#' ssvFeatureVenn(list(1:3, 2:6),
#'   counts_as_percent = TRUE,
#'   percentage_digits = 0,
#'   percentage_suffix = " %",
#'   counts_txt_size = 12)
ssvFeatureVenn = function(object,
                          group_names = NULL,
                          counts_txt_size = 5,
                          counts_as_labels = FALSE,
                          show_outside_count = FALSE,
                          line_width = 3,
                          circle_colors = NULL,
                          fill_alpha = .3,
                          line_alpha = 1,
                          counts_color = NULL,
                          counts_as_percent = FALSE,
                          percentage_digits = 1,
                          percentage_suffix = "%",
                          n_points = 200,
                          return_data = FALSE) {
    size = group = x = y = label = NULL#declare binding for data.table

    object = ssvMakeMembTable(object)
    stopifnot(is.numeric(counts_txt_size), is.numeric(line_width),
              is.numeric(fill_alpha))
    stopifnot(is.logical(counts_as_labels), is.logical(show_outside_count))
    all(apply(object, 2, class) == "logical")
    object <- limma::vennCounts(object)
    set_counts <- object[, "Counts"]
    if(counts_as_percent){
        set_counts = set_counts / sum(set_counts)
        set_counts = paste0(round(set_counts * 100, digits = percentage_digits), percentage_suffix)
    }

    nsets <- ncol(object) - 1
    if (nsets > 3)
        stop("Can't plot Venn diagram for more than 3 sets")
    if (is.null(group_names))
        group_names <- factor(colnames(object)[seq_len(nsets)],
                              levels = colnames(object)[seq_len(nsets)])
    if (is.null(counts_color))
        counts_color <- grDevices::rgb(0,0,0)

    xcentres <- switch(nsets,
                       0,
                       c(-1, 1),
                       c(-1, 1, 0),
                       c(1, 2, 1, 2)
    )
    ycentres <- switch(nsets,
                       0,
                       c(0, 0),
                       c(1, 1, -2)/sqrt(3)
    )
    r <- 1.5
    p = ggellipse(xcentres = xcentres,
                  ycentres = ycentres,
                  r = r,
                  circle_colors = circle_colors,
                  group_names = group_names,
                  line_alpha =  line_alpha,
                  fill_alpha = fill_alpha,
                  line_width = line_width,
                  n_points = n_points)

    df_text <- switch(nsets,
                      {
                          df = data.frame(x = 0,
                                          y = 0,
                                          label = set_counts[-1],
                                          size = counts_txt_size,
                                          col = counts_color)
                          if (show_outside_count) {
                              df = rbind(df, data.frame(x = 2.3, y = -2.1,
                                                        label = set_counts[1],
                                                        size = counts_txt_size,
                                                        col = counts_color))
                          }
                          df
                      }, {
                          df = data.frame(x = c(1.5, -1.5, 0),
                                          y = c(0.1, 0.1, 0.1),
                                          label = set_counts[-1],
                                          size = counts_txt_size,
                                          col = counts_color)
                          if (show_outside_count) {
                              df = rbind(df, data.frame(x = 2.3,
                                                        y = -2.1,
                                                        label = set_counts[1],
                                                        size = counts_txt_size,
                                                        col = counts_color))
                          }
                          df
                      }, {
                          df = data.frame(x = c(0, 1.5, 0.75, -1.5,
                                                -0.75, 0, 0),
                                          y = c(-1.7, 1, -0.35, 1, -0.35,
                                                0.9, 0),
                                          label = set_counts[-1],
                                          size = counts_txt_size,
                                          col = counts_color)
                          if (show_outside_count) {
                              df = rbind(df, data.frame(x = 2.5, y = -3,
                                                        label = set_counts[1],
                                                        size = counts_txt_size,
                                                        col = counts_color))
                          }
                          df
                      })

    counts_method = ifelse(counts_as_labels, geom_label, geom_text)
    p = p + counts_method(data = df_text,
                          mapping = aes(x = x,
                                        y = y,
                                        label = label,
                                        size = size),
                          show.legend = FALSE)
    if(return_data){
        return(
            data.table(xcentres = xcentres,
                       ycentres = ycentres,
                       r = r,
                       r2 = r,
                       phi = rep(0, length(xcentres)),
                       group_names = group_names)
        )
    }
    return(p)
}


#' ssvFeatureUpset
#'
#' Uses the UpSetR package to create an [UpSetR::upset] plot of region overlaps.
#'
#' @param object will be passed to \code{\link{ssvMakeMembTable}} for conversion
#'   to membership matrix
#' @param return_UpSetR If TRUE, return the UpSetR object, The default is FALSE
#'   and results in a ggplotified version compatible with cowplot etc.
#' @param nsets Number of sets to look at
#' @param nintersects Number of intersections to plot. If set to NA, all
#'   intersections will be plotted.
#' @param order.by How the intersections in the matrix should be ordered by.
#'   Options include frequency (entered as "freq"), degree, or both in any
#'   order.
#' @param ... Additional parameters passed to \code{\link[UpSetR]{upset}} in the
#'   UpSetR  package.
#'
#' @return ggplot version of UpSetR plot
#' @export
#' @import UpSetR
#' @import ggplotify
#'
#' @examples
#' ssvFeatureUpset(list(1:3, 2:6))
#' ssvFeatureUpset(CTCF_in_10a_overlaps_gr)
#' ssvFeatureUpset(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
ssvFeatureUpset = function(object,
                           return_UpSetR = FALSE,
                           nsets = NULL,
                           nintersects = 15,
                           order.by = "freq",
                           ...){
    object = ssvMakeMembTable(object)
    up_df = as.matrix(object)
    up_df = ifelse(up_df, 1, 0)
    up_df = as.data.frame(up_df)
    if(is.null(nsets)) nsets = ncol(up_df)
    if(ncol(up_df) > 1){
        p_up = UpSetR::upset(up_df,
                             nsets = nsets,
                             nintersects = nintersects,
                             order.by = order.by,
                             ...)
    }else{
        p_up = ggplot() + theme_void() +
            labs(title = "cannot run UpSetR on only 1 set.")
    }
    if(return_UpSetR){
        return(p_up)
    }
    p = ggplotify::as.ggplot(p_up)
    p
}


#' Try to load a bed-like file and convert it to a GRanges object
#'
#' @export
#' @param object A membership table
#' @param line_width numeric, passed to size aesthetic to control line width
#' @param shape shape argument passed to eulerr::euler
#' @param n_points number of points to use for drawing ellipses, passed to
#' eulerr:::ellipse
#' @param fill_alpha numeric value from 0 to 1. Alpha value for circle fill
#' @param line_alpha numeric value from 0 to 1. Alpha value for circle line
#' @param circle_colors colors to choose from for circles.  passed to ggplot2
#' color scales.
#' @param return_data logical.  If TRUE, return value is no longer ggplot and
#' is instead the data used to generate that plot. Default is FALSE.
#' @return ggplot of venneuler results
#' @import ggplot2
#' @import eulerr
#' @import S4Vectors
#' @examples
#' ssvFeatureEuler(list(1:3, 2:6))
#' ssvFeatureEuler(CTCF_in_10a_overlaps_gr)
#' ssvFeatureEuler(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
ssvFeatureEuler = function(object,
                           line_width = 2,
                           shape = c("circle", "ellipse")[1],
                           n_points = 200,
                           fill_alpha = .3,
                           line_alpha = 1,
                           circle_colors = NULL,
                           return_data = FALSE) {
    x = y = group = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    stopifnot(is.numeric(line_width),
              is.numeric(n_points),
              is.numeric(fill_alpha))
    stopifnot(is.character(shape))
    stopifnot(shape %in% c("circle", "ellipse"))
    cn = colnames(object)
    eu = eulerr::euler(object, shape = shape)
    dd = eu$ellipses

    h <- dd$h
    k <- dd$k
    a <- dd$a
    b <- dd$b
    phi <- dd$phi

    if(return_data){
        return(
            data.table(xcentres = h,
                       ycentres = k,
                       r = a,
                       r2 = b,
                       phi = phi,
                       group_names = cn)
        )
    }

    p = ggellipse(xcentres = h,
                  ycentres = k,
                  r = a,
                  r2 = b,
                  phi = phi,
                  circle_colors = circle_colors,
                  group_names = cn,
                  line_alpha =  line_alpha,
                  fill_alpha = fill_alpha,
                  line_width = line_width,
                  n_points = n_points)
    p
}


#' bar plots of set sizes
#' @export
#' @param object passed to ssvMakeMembTable for conversion to membership table
#' @param show_counts logical.  should counts be displayed at the center of each
#' bar. default is TRUE
#' @param bar_colors character. rcolor or hex colors. default of NULL
#' uses RColorBrewer Dark2. Will repeat to match number of samples.
#' @param counts_text_colors character. rcolor or hex colors. default of NULL
#' uses RColorBrewer Dark2. Will repeat to match number of samples.
#' @param return_data logical.  If TRUE, return value is no longer ggplot and
#' is instead the data used to generate that plot. Default is FALSE.
#' @param count_label_size Font size bar count labels. Default is 8.
#' @return ggplot of bar plot of set sizes
#' @import ggplot2
#' @import S4Vectors
#' @examples
#' ssvFeatureBars(list(1:3, 2:6))
#' ssvFeatureBars(CTCF_in_10a_overlaps_gr, count_label_size = 10)
#' ssvFeatureBars(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
ssvFeatureBars = function(object,
                          show_counts = TRUE,
                          bar_colors = NULL,
                          counts_text_colors = NULL,
                          return_data = FALSE,
                          count_label_size = 8) {
    group = count = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    stopifnot(is.logical(show_counts))
    stopifnot(is.null(bar_colors) || all(is.character(bar_colors)))
    n_bars = ncol(object)
    #bar colors
    if (is.null(bar_colors)) {
        bar_colors = safeBrew(n_bars, "Dark2")
    } else {
        not_hex = substr(bar_colors, 0, 1) != "#"
        bar_colors[not_hex] = col2hex(bar_colors[not_hex])
    }
    if (length(bar_colors) < n_bars)
        bar_colors <- rep(bar_colors, length.out = n_bars)
    names(bar_colors) = colnames(object)
    #text colors
    if (is.null(counts_text_colors)) {
        counts_text_colors = rep("black", n_bars)
    } else {
        not_hex = substr(counts_text_colors, 0, 1) != "#"
        counts_text_colors[not_hex] = col2hex(counts_text_colors[not_hex])
    }
    if (length(counts_text_colors) < n_bars)
        counts_text_colors <- rep(counts_text_colors, length.out = n_bars)
    names(counts_text_colors) = colnames(object)

    hit_counts = colSums(object)
    hit_counts_df = data.frame(count = hit_counts,
                               group = factor(names(hit_counts),
                                              levels = names(hit_counts)))
    if(return_data){
        return(as.data.table(hit_counts_df))
    }
    p <- ggplot(hit_counts_df, aes(x = group, y = count, fill = group)) +
        labs(x = "") +
        geom_bar(stat = "identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank()) +
        scale_fill_manual(values = bar_colors)
    if(show_counts)
        p = p + annotate("text",
                         y = hit_counts/2,
                         x = seq_along(hit_counts),
                         label = hit_counts,
                         color = counts_text_colors,
                         size = count_label_size/ggplot2::.pt)
    return(p)
}

#' ssvFeaturePie
#'
#' Generate a ggplot pie plot of set sizes.
#'
#' @export
#' @param object object that ssvMakeMembTable can convert to logical matrix
#'   membership
#' @param slice_colors colors to use for pie slices
#' @param return_data logical.  If TRUE, return value is no longer ggplot and is
#'   instead the data used to generate that plot. Default is FALSE.
#' @import ggplot2
#' @import S4Vectors
#' @return ggplot pie graph of set sizes
#' @examples
#' ssvFeaturePie(list(1:3, 2:6))
#' ssvFeaturePie(CTCF_in_10a_overlaps_gr)
#' ssvFeaturePie(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
ssvFeaturePie = function(object, slice_colors = NULL, return_data = FALSE) {
    count = group = NULL#declare binding for data.table

    stopifnot(is.null(slice_colors) || all(is.character(slice_colors)))

    object = ssvMakeMembTable(object)
    hit_counts = colSums(object)
    hit_counts_df = data.frame(count = hit_counts,
                               group = factor(names(hit_counts),
                                              levels = names(hit_counts)))

    n_slices = ncol(object)
    if (is.null(slice_colors)) {
        slice_colors = safeBrew(n_slices, "Dark2")
    } else {
        not_hex = substr(slice_colors, 0, 1) != "#"
        slice_colors[not_hex] = col2hex(slice_colors[not_hex])
    }
    if (length(slice_colors) < n_slices)
        slice_colors <- rep(slice_colors, length.out = n_slices)
    names(slice_colors) = colnames(object)

    if(return_data){
        return(as.data.table(hit_counts_df))
    }

    p <- ggplot(hit_counts_df, aes(x = "", y = count, fill = group)) +
        labs(x = "") +
        geom_bar(width = 1, stat = "identity") +
        guides(x = "none") +
        theme(axis.text.x = element_blank(),
              axis.line = element_blank(),
              panel.background = element_blank(),
              axis.ticks = element_blank()) +
        coord_polar("y", start = 0) +
        scale_y_reverse() +
        scale_fill_manual(values = slice_colors)
    hc = rev(hit_counts)/sum(hit_counts)
    hc = c(0, cumsum(hc)[-length(hc)]) + hc/2
    p = p + annotate(geom = "text",
                     y = hc * sum(hit_counts),
                     x = 1.1,
                     label = rev(hit_counts_df$count))
    return(p)
}



#' ssvFeatureBinaryHeatmap
#'
#' Outputs a ggplot binary heatmap, where color indicates TRUE and the other
#' indicates FALSE in a membership table. The heatmap is sorted, TRUE at the
#' top, by column left to right. Changes to column order can reveal different
#' patterns.
#'
#' As a svg output, the final plot can be unwieldy. The default of
#' `raster_approximation` = TRUE is easier to work with, especially for larger
#' membership tables.
#'
#' @export
#' @param object passed to ssvMakeMembTable
#' @param raster_approximation If TRUE, instead of standard ggplot, write
#' temporary raster png image and redraw that as plot background. default is
#' FALSE
#' @param true_color character. rcolor or hex color used for TRUE values.
#' default is "black".
#' @param false_color character. rcolor or hex color used for TRUE values.
#' default is "#EFEFEF", a gray.
#' @param raster_width_min raster width will be minimum multiple of number of
#' columns over this number.  ignored if raster_approximation is FALSE.
#' @param raster_height_min raster height will be minimum multiple of number of
#' rows over this number  ignored if raster_approximation is FALSE
#' @param return_data logical.  If TRUE, return value is no longer ggplot and
#' is instead the data used to generate that plot. Default is TRUE
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @return ggplot using geom_tile of membership table sorted from left to right.
#' @import png
#' @import ggplot2
#' @importFrom grid rasterGrob
#' @examples
#' ssvFeatureBinaryHeatmap(list(1:3, 2:6))
#' # horizontal version
#' ssvFeatureBinaryHeatmap(list(1:3, 2:6)) + coord_flip() +
#'   theme(axis.text.x = element_blank(), axis.text.y = element_text())
#' ssvFeatureBinaryHeatmap(CTCF_in_10a_overlaps_gr)
#' ssvFeatureBinaryHeatmap(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
#' ssvFeatureBinaryHeatmap(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,3:2])
ssvFeatureBinaryHeatmap = function(object,
                                   raster_approximation = TRUE,
                                   true_color = "black",
                                   false_color = "#EFEFEF",
                                   raster_width_min = 1000,
                                   raster_height_min = 1000,
                                   return_data = FALSE) {
    groups = bool = value = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    stopifnot(is.logical(raster_approximation))

    mat = ifelse(as.matrix(object), 1, 0)
    for (i in rev(seq_len(ncol(mat)))){
        mat = mat[order(mat[, i], decreasing = FALSE), , drop = FALSE]
    }
    hdt = data.table::as.data.table(cbind(mat, row = seq_len(nrow(mat))))
    hdt = data.table::melt(hdt, id.vars = "row", variable.name = "groups")

    if (raster_approximation) {
        stopifnot(is.numeric(raster_width_min), is.numeric(raster_height_min))
        dmat = as.matrix(
            data.table::dcast(
                hdt,
                value.var = "value",
                formula = rev(row) ~ groups))[, 1 + seq_len(ncol(object)),
                                              drop = FALSE]
        low_v = 1
        dmat = ifelse(dmat > low_v, low_v, dmat)
        # dmat = -1*(dmat - low_v)
        if (raster_width_min >= ncol(dmat)) {
            # expand cols as necessary to meet target
            raster_width_multiplier = ceiling(raster_width_min/ncol(dmat))
            comp_col = rep(seq_len(ncol(dmat)), each = raster_width_multiplier)
        } else {
            # evenly spaced sample down to target
            comp_col = round(seq_len(raster_width_min)/raster_width_min *
                                 ncol(dmat))
        }
        if (raster_height_min >= nrow(dmat)) {
            # expand rows as necessary to meet target
            raster_height_multiplier = ceiling(raster_height_min/nrow(dmat))
            comp_row = rep(seq_len(nrow(dmat)), each = raster_height_multiplier)
        } else {
            # evenly spaced sample down to target
            comp_row = round(seq_len(raster_height_min)/raster_height_min *
                                 nrow(dmat))
        }
        # dmat = (dmat * -.95 + .95)
        comp_mat = dmat[comp_row, comp_col]
        comp_ar = array(comp_mat, c(nrow(comp_mat), ncol(comp_mat), 3))
        comp_ar[, , 1] = ifelse(comp_mat == 1,
                                col2rgb(true_color)[1,1]/255,
                                col2rgb(false_color)[1,1]/255)
        comp_ar[, , 2] = ifelse(comp_mat == 1,
                                col2rgb(true_color)[2,1]/255,
                                col2rgb(false_color)[2,1]/255)
        comp_ar[, , 3] = ifelse(comp_mat == 1,
                                col2rgb(true_color)[3,1]/255,
                                col2rgb(false_color)[3,1]/255)
        png::writePNG(comp_ar, target = "tmp.png")
        png_grob = grid::rasterGrob(png::readPNG("tmp.png"),
                                    width = unit(1, "npc"),
                                    height = unit(1, "npc"),
                                    interpolate = FALSE)
        file.remove("tmp.png")
        p = ggplot(hdt) +
            geom_tile(aes(x = groups, y = row, fill = NA, col = NULL)) +
            scale_fill_manual(values = c(`TRUE` = "black",
                                         `FALSE` = "#EFEFEF")) +
            scale_alpha(0) +
            labs(fill = "", y = "") +
            guides(fill = "none") +
            theme(
                axis.line = element_blank(),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                plot.background = element_blank(),
                panel.background = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title = element_blank()
            ) +
            annotation_custom(png_grob, xmin = 0.5, xmax = ncol(dmat) + 0.5,
                              ymin = 0.5, ymax = nrow(dmat) + 0.5)
    } else {
        hdt[, `:=`(bool, value == 1)]
        p = ggplot(hdt) +
            geom_tile(aes(x = groups, y = row, fill = bool, col = NULL)) +
            scale_fill_manual(values = c(`TRUE` = true_color,
                                         `FALSE` = false_color)) +
            labs(fill = "", y = "") +
            guides(fill = "none") +
            theme(
                axis.line = element_blank(),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                plot.background = element_blank(),
                panel.background = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title = element_blank()
            )
    }
    if(return_data){
        return(hdt)
    }
    return(p)
}
