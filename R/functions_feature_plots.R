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

#' ggplot implementation of vennDiagram from limma package.  currently limited at 3 sets
#' @export
#' @param object will be passed to \code{\link{ssvMakeMembTable}} for conversion to membership matrix
#' @param group_names useful if names weren't provided or were lost in creating membership matrix
#' @param counts_txt_size font size for count numbers
#' @param counts_as_labels if TRUE, geom_label is used instead of geom_text.  can be easier to read.
#' @param show_outside_count if TRUE, items outside of all sets are counted outside. can be confusing.
#' @param line_width uses size aesthetic to control line width of circles.
#' @param circle_colors colors to use for circle line colors. Uses Dark2 set from RColorBrewer by default.
#' @param fill_alpha alpha value to use for fill, defaults to .3.
#' @param line_alpha numeric [0,1], alpha value for circle line
#' @param counts_color character. single color to use for displaying counts
#' @param n_points integer.  number of points to approximate circle with.
#' default is 200.
#' @return ggplot venn diagram
#' @import ggplot2
#' @import S4Vectors
#' @importFrom limma vennCounts
#' @examples
#' ssvFeatureVenn(list(1:3, 2:6))
#' ssvFeatureVenn(CTCF_in_10a_overlaps_gr)
#' ssvFeatureVenn(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
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
                          n_points = 200) {
    size = group = x = y = label = NULL#declare binding for data.table

    object = ssvMakeMembTable(object)
    stopifnot(is.numeric(counts_txt_size), is.numeric(line_width),
              is.numeric(fill_alpha))
    stopifnot(is.logical(counts_as_labels), is.logical(show_outside_count))
    all(apply(object, 2, class) == "logical")
    object <- limma::vennCounts(object)
    set_counts <- object[, "Counts"]
    nsets <- ncol(object) - 1
    if (nsets > 3)
        stop("Can't plot Venn diagram for more than 3 sets")
    if (is.null(group_names))
        group_names <- factor(colnames(object)[seq_len(nsets)], levels = colnames(object)[seq_len(nsets)])
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
                                          label = set_counts[-1], size = counts_txt_size, col = counts_color)
                          if (show_outside_count) {
                              df = rbind(df, data.frame(x = 2.3, y = -2.1, label = set_counts[1], size = counts_txt_size, col = counts_color))
                          }
                          df
                      }, {
                          df = data.frame(x = c(1.5, -1.5, 0),
                                          y = c(0.1, 0.1, 0.1),
                                          label = set_counts[-1], size = counts_txt_size, col = counts_color)
                          if (show_outside_count) {
                              df = rbind(df, data.frame(x = 2.3, y = -2.1, label = set_counts[1], size = counts_txt_size, col = counts_color))
                          }
                          df
                      }, {
                          df = data.frame(x = c(0, 1.5, 0.75, -1.5, -0.75, 0, 0),
                                          y = c(-1.7, 1, -0.35, 1, -0.35, 0.9, 0),
                                          label = set_counts[-1],
                                          size = counts_txt_size, col = counts_color)
                          if (show_outside_count) {
                              df = rbind(df, data.frame(x = 2.5, y = -3, label = set_counts[1], size = counts_txt_size, col = counts_color))
                          }
                          df
                      })

    counts_method = ifelse(counts_as_labels, geom_label, geom_text)
    p = p + counts_method(data = df_text, mapping = aes(x = x, y = y, label = label, size = size))
    return(p)
}

# from https://gist.github.com/trinker/31edc08d0a4ec4c73935a23040c2f6cb p_load(dplyr, venneuler, ggforce, textshape)
#' Try to load a bed-like file and convert it to a GRanges object
#' @export
#' @param object A membership table
#' @param line_width numeric, passed to size aesthetic to control line width
#' @param shape shape argument passed to eulerr::euler
#' @param n_points number of points to use for drawing ellipses, passed to  eulerr:::ellipse
#' @param fill_alpha numeric [0,1], alpha value for circle fill
#' @param line_alpha numeric [0,1], alpha value for circle line
#' @param circle_colors colors to choose from for circles.  passed to ggplot2 color scales.
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
                           circle_colors = NULL) {
    x = y = group = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    stopifnot(is.numeric(line_width), is.numeric(n_points), is.numeric(fill_alpha))
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
#' uses RColorBrewer Dark2.
#' @return ggplot of bar plot of set sizes
#' @import ggplot2
#' @import S4Vectors
#' @examples
#' ssvFeatureBars(list(1:3, 2:6))
#' ssvFeatureBars(CTCF_in_10a_overlaps_gr)
#' ssvFeatureBars(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
ssvFeatureBars = function(object, show_counts = TRUE, bar_colors = NULL) {
    group = count = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    stopifnot(is.logical(show_counts))
    stopifnot(is.null(bar_colors) || all(is.character(bar_colors)))
    n_bars = ncol(object)
    if (is.null(bar_colors)) {
        bar_colors = safeBrew(n_bars, "Dark2")
    } else {
        not_hex = substr(bar_colors, 0, 1) != "#"
        bar_colors[not_hex] = col2hex(bar_colors[not_hex])
    }
    if (length(bar_colors) < n_bars)
        bar_colors <- rep(bar_colors, length.out = n_bars)
    names(bar_colors) = colnames(object)
    hit_counts = colSums(object)
    hit_counts_df = data.frame(count = hit_counts, group = factor(names(hit_counts), levels = names(hit_counts)))
    p <- ggplot(hit_counts_df, aes(x = group, y = count, fill = group)) +
        labs(x = "") +
        geom_bar(stat = "identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank()) +
        scale_fill_manual(values = bar_colors)
    if(show_counts)
        p = p + annotate("text", y = hit_counts/2, x = seq_along(hit_counts), label = hit_counts)
    return(p)
}


# bar plots of set sizes, does fancy facetting
#
# @param object passed to ssvMakeMembTable for conversion to membership table
# @param peak_gr peak sets before overlapping
# @param flip_group_stage logical, different facetting.  default is FALSE
#
# @return ggplot of bar plot of set sizes
# ssvFeatureBars2 = function(object, peak_gr, flip_group_stage = FALSE) {
#   x = count = group = stage = NULL#declare binding for data.table
#   object = ssvMakeMembTable(object)
#   hit_counts = colSums(object)
#   hit_counts_df = data.frame(count = hit_counts, group = factor(names(hit_counts), levels = names(hit_counts)), stage = "merged")
#   raw_counts = sapply(peak_gr, length)
#   raw_counts_df = data.frame(count = raw_counts, group = factor(names(raw_counts), levels = names(raw_counts)), stage = "raw")
#   counts_df = rbind(hit_counts_df, raw_counts_df)
#   counts_df$stage = factor(counts_df$stage, levels = c("raw", "merged"))
#   counts_df$x = paste(counts_df$group, counts_df$stage)
#   counts_df$x = factor(counts_df$x, levels = counts_df$x[order(counts_df$stage)][order(counts_df$group)])
#   bar_width = 0.8
#   if (!flip_group_stage) {
#     p <- ggplot(counts_df, aes(x = x, y = count, fill = group)) + labs(x = "") + geom_bar(width = bar_width, stat = "identity",
#                                                                                           position = "dodge") + scale_x_discrete(labels = counts_df$stage[order(counts_df$x)]) + theme_bw() + theme(axis.text.x = element_text(angle = 90,
#                                                                                                                                                                                                                                hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
#     p = p + annotate("text", y = (counts_df$count[order(counts_df$x)])/2, x = seq_along(counts_df$count), label = counts_df$count[order(counts_df$x)])
#   } else {
#     p <- ggplot(counts_df, aes(x = group, y = count, fill = stage)) + labs(x = "") + geom_bar(width = bar_width, stat = "identity",
#                                                                                               position = "dodge") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank(),
#                                                                                                                                        panel.grid.minor.x = element_blank())
#     xs = rep(seq_along(levels(counts_df$group)), each = 2) + rep(c(-1, 1), length(levels(counts_df$group))) * bar_width/4
#     p = p + annotate("text", y = (counts_df$count[order(counts_df$x)])/2, x = xs, label = counts_df$count[order(counts_df$x)])
#   }
#   {
#     p <- ggplot(counts_df, aes(x = stage, y = count, fill = group)) + labs(x = "") + geom_bar(width = bar_width, stat = "identity") +
#       theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank(),
#                          panel.grid.minor.x = element_blank())
#   }
#   p
#   return(p)
# }

# pie plot with facetting
# @param object passed to ssvMakeMembTable
# @param peak_gr peak sets before overlapping
#
# @return ggplot of pie charts
# ssvFeaturePie2 = function(object, peak_gr) {
#   stage = count = group = NULL#declare binding for data.table
#   object = ssvMakeMembTable(object)
#   hit_counts = colSums(object)
#   hit_counts_df = data.frame(count = hit_counts, group = factor(names(hit_counts), levels = names(hit_counts)), stage = "merged")
#   raw_counts = sapply(peak_gr, length)
#   raw_counts_df = data.frame(count = raw_counts, group = factor(names(raw_counts), levels = names(raw_counts)), stage = "raw")
#   counts_df = rbind(hit_counts_df, raw_counts_df)
#   counts_df$stage = factor(counts_df$stage, levels = c("raw", "merged"))
#   counts_df$x = paste(counts_df$group, counts_df$stage)
#   counts_df$x = factor(counts_df$x, levels = counts_df$x[order(counts_df$stage)][order(counts_df$group)])
#
#   counts_df = counts_df[rev(order(counts_df$group)), ]
#   counts_df = counts_df[order(counts_df$stage), ]
#   p <- ggplot(counts_df, aes(x = stage, y = count, fill = (group))) + labs(x = "") + geom_bar(width = 1, stat = "identity") +
#     guides(x = "none") + theme(axis.text.x = element_blank(), panel.background = element_blank(), axis.ticks = element_blank()) +
#     coord_polar("y", start = 0) + scale_y_reverse()
#   ys = unlist(lapply(levels(counts_df$stage), function(x) {
#     stage_df = counts_df[counts_df$stage == x, ]
#     hc = (stage_df$count)/sum(stage_df$count)
#     hc = c(0, cumsum(hc)[-length(hc)]) + hc/2
#     hc * sum(stage_df$count)
#   }))
#   xs = rep(c(1.1, 2.1), each = length(levels(counts_df$group)))
#   txt = counts_df$count
#   p = p + annotate(geom = "text", y = ys, x = xs, label = txt)
#   return(p)
# }

#' pie plot of set sizes
#' @export
#' @param object object that ssvMakeMembTable can convert to logical matrix membership
#' @import ggplot2
#' @import S4Vectors
#' @return ggplot pie graph of set sizes
#' @examples
#' ssvFeaturePie(list(1:3, 2:6))
#' ssvFeaturePie(CTCF_in_10a_overlaps_gr)
#' ssvFeaturePie(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
ssvFeaturePie = function(object) {
    count = group = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    hit_counts = colSums(object)
    hit_counts_df = data.frame(count = hit_counts, group = factor(names(hit_counts), levels = names(hit_counts)))

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
        scale_fill_brewer(palette = "Dark2")
    hc = rev(hit_counts)/sum(hit_counts)
    hc = c(0, cumsum(hc)[-length(hc)]) + hc/2
    p = p + annotate(geom = "text", y = hc * sum(hit_counts), x = 1.1, label = rev(hit_counts_df$count))
    return(p)
}



#' binary heatmap indicating membership.
#' heatmap is sorted by column left to right.
#' change column order to reveal patterns
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
#' @rawNamespace import(data.table, except = c(shift, first, second))
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
ssvFeatureBinaryHeatmap = function(object, raster_approximation = FALSE,
                                   true_color = "black",
                                   false_color = "#EFEFEF",
                                   raster_width_min = 1000, raster_height_min = 1000) {
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
        dmat = as.matrix(data.table::dcast(hdt, value.var = "value", formula = rev(row) ~ groups))[, 1 + seq_len(ncol(object)), drop = FALSE]
        low_v = 1
        dmat = ifelse(dmat > low_v, low_v, dmat)
        # dmat = -1*(dmat - low_v)
        if (raster_width_min >= ncol(dmat)) {
            # expand cols as necessary to meet target
            raster_width_multiplier = ceiling(raster_width_min/ncol(dmat))
            comp_col = rep(seq_len(ncol(dmat)), each = raster_width_multiplier)
        } else {
            # evenly spaced sample down to target
            comp_col = round(seq_len(raster_width_min)/raster_width_min * ncol(dmat))
        }
        if (raster_height_min >= nrow(dmat)) {
            # expand rows as necessary to meet target
            raster_height_multiplier = ceiling(raster_height_min/nrow(dmat))
            comp_row = rep(seq_len(nrow(dmat)), each = raster_height_multiplier)
        } else {
            # evenly spaced sample down to target
            comp_row = round(seq_len(raster_height_min)/raster_height_min * nrow(dmat))
        }
        # dmat = (dmat * -.95 + .95)
        comp_mat = dmat[comp_row, comp_col]
        comp_ar = array(comp_mat, c(nrow(comp_mat), ncol(comp_mat), 3))
        comp_ar[, , 1] = ifelse(comp_mat == 1, col2rgb(true_color)[1,1]/255, col2rgb(false_color)[1,1]/255)
        comp_ar[, , 2] = ifelse(comp_mat == 1, col2rgb(true_color)[2,1]/255, col2rgb(false_color)[2,1]/255)
        comp_ar[, , 3] = ifelse(comp_mat == 1, col2rgb(true_color)[3,1]/255, col2rgb(false_color)[3,1]/255)
        png::writePNG(comp_ar, target = "tmp.png")
        png_grob = grid::rasterGrob(png::readPNG("tmp.png"),
                                    width = unit(1, "npc"),
                                    height = unit(1, "npc"),
                                    interpolate = FALSE)
        file.remove("tmp.png")
        p = ggplot(hdt) +
            geom_tile(aes(x = groups, y = row, fill = NA, col = NULL)) +
            scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "#EFEFEF")) +
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
    return(p)
}
