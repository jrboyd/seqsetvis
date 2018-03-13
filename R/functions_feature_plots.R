# plotting functions for set comparison go here
# all plotting functions accept a variety of objects as first argument and
# handle them via the S4  generic ssvMakeMembTable
# ssvMakeMembTable has implemented methods for:
#   - membership table
#       columns are sets
#       rows are members
#       values are logicals indicating membership
#   - list containing character vectors of set members
#   - TODO GrangesOverlapSets object
#   - TODO list of Granges to convert to GRangesOverlapSets

#' ggplot implementation of vennDiagram from limma package.  currently limited at 3 sets
#' @export
#' @param object will be passed to \code{\link{ssvMakeMembTable}} for conversion to membership matrix
#' @param group_names useful if names weren't provided or were lost in creating membership matrix
#' @param counts_txt_size font size for count numbers
#' @param counts_as_labels if TRUE, geom_label is used instead of geom_text.  can be easier to read.
#' @param show_outside_count if TRUE, items outside of all sets are counted outside. can be confusing.
#' @param lwd uses size aesthetic to control line width of circles.
#' @param circle_color colors to use for circle line colors. Uses Dark2 set from RColorBrewer by default.
#' @param fill_circles if TRUE, fill circles
#' @param fill_alpha alpha value to use for fill, defaults to .5.
#' @param counts_color single color to use for displaying counts
#' @return ggplot venn diagram
#' @import ggplot2
#' @importFrom limma vennCounts
#' @importFrom ggforce geom_circle
#' @examples
#' ssvFeatureVenn(list(1:3, 2:6))
#' ssvFeatureVenn(CTCF_in_10a_overlaps_gr)
#' ssvFeatureVenn(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
ssvFeatureVenn = function(object, group_names = NULL, counts_txt_size = 5,
                          counts_as_labels = FALSE, show_outside_count = FALSE, lwd = 3,
                          circle_color = NULL, fill_circles = TRUE,
                          fill_alpha = ifelse(fill_circles, 0.5, 0), counts_color = NULL) {
    size = group = x = y = label = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    all(apply(object, 2, class) == "logical")
    object <- limma::vennCounts(object)
    set_counts <- object[, "Counts"]
    nsets <- ncol(object) - 1
    if (nsets > 3)
        stop("Can't plot Venn diagram for more than 3 sets")
    if (is.null(group_names))
        group_names <- factor(colnames(object)[1:nsets], levels = colnames(object)[1:nsets])
    if (is.null(circle_color)) {
        circle_color = safeBrew(nsets, "Dark2")
    } else {
        not_hex = substr(circle_color, 0, 1) != "#"
        circle_color[not_hex] = col2hex(circle_color[not_hex])
    }
    if (length(circle_color) < nsets)
        circle_color <- rep(circle_color, length.out = nsets)
    if (is.null(counts_color))
        counts_color <- grDevices::rgb(0,0,0)
    col_scale = circle_color
    names(col_scale) = group_names
    ahex = substr(grDevices::rgb(red = 1, blue = 1, green = 1, alpha = fill_alpha), start = 8, stop = 9)
    fill_scale = paste0(col_scale, ahex)
    names(fill_scale) = names(col_scale)

    p = ggplot() +
        labs(fill = "", color = "") +
        scale_color_manual(values = col_scale) +
        scale_size_identity() +
        scale_fill_manual(values = fill_scale) +
        theme_minimal() +
        theme(plot.background = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              legend.position = "top") +
        guides(fill = guide_legend(override.aes = list(shape = 21)))
    p = p + coord_fixed()
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
    # xtext <- switch(nsets,
    #                 -1.2,
    #                 c(-1.2, 1.2),
    #                 c(-1.2, 1.2, 0))
    # ytext <- switch(nsets,
    #                 1.8,
    #                 c(1.8, 1.8),
    #                 c(2.4, 2.4, -3))

    df_circles = data.frame(xcentres, ycentres, r, size = lwd, col = circle_color, group = group_names)


    p = p + ggforce::geom_circle(data = df_circles, mapping = aes(x0 = xcentres, y0 = ycentres, r = r, size = size, col = group,
                                                                  fill = group))
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
    p = p + counts_method(data = df_text, mapping = aes(x = x, y = y, label = label, size = counts_txt_size))
    return(p)
}

#' bar plots of set sizes
#' @export
#' @param object passed to ssvMakeMembTable for conversion to membership table
#'
#' @return ggplot of bar plot of set sizes
#' @import ggplot2
#' @examples
#' ssvFeatureBars(list(1:3, 2:6))
#' ssvFeatureBars(CTCF_in_10a_overlaps_gr)
#' ssvFeatureBars(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
ssvFeatureBars = function(object) {
    group = count = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    hit_counts = colSums(object)
    hit_counts_df = data.frame(count = hit_counts, group = factor(names(hit_counts), levels = names(hit_counts)))
    p <- ggplot(hit_counts_df, aes(x = group, y = count, fill = group)) +
        labs(x = "") +
        geom_bar(width = 1, stat = "identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank()) +
        scale_fill_brewer(palette = "Dark2")
    p = p + annotate("text", y = hit_counts/2, x = 1:length(hit_counts), label = hit_counts)
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

# from https://gist.github.com/trinker/31edc08d0a4ec4c73935a23040c2f6cb p_load(dplyr, venneuler, ggforce, textshape)
#' Try to load a bed-like file and convert it to a GRanges object
#' @export
#' @param object A membership table
#' @param line_width numeric, passed to size aesthetic to control line width
#' @param shape shape argument passed to eulerr::euler
#' @param n number of points to use for drawing ellipses, passed to  eulerr:::ellipse
#' @param fill_circles logical, should circles be filled or empty?
#' @param alpha numeric [0,1], alpha value for circle fill
#' @param col_scale colors to choose from for circles.  passed to ggplot2 color scales.
#' @return ggplot of venneuler results
#' @import ggplot2
#' @import eulerr
#' @examples
#' ssvFeatureEuler(list(1:3, 2:6))
#' ssvFeatureEuler(CTCF_in_10a_overlaps_gr)
#' ssvFeatureEuler(S4Vectors::mcols(CTCF_in_10a_overlaps_gr)[,2:3])
ssvFeatureEuler = function(object, line_width = 2, shape = c("circle", "ellipse")[1],
                           n = 200, fill_circles = T, alpha = .15,
                           col_scale = NULL) {
    x = y = group = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    cn = colnames(object)
    eu = eulerr::euler(object, shape = shape)
    dd = eu$ellipses

    h <- dd$h
    k <- dd$k
    a <- dd$a
    b <- dd$b
    phi <- dd$phi
    # internal function from eulerr
    eulerr_ellipse = function (h, k, a, b = a, phi = 0, n = 200L) {
        theta <- seq.int(0, 2 * pi, length.out = n)
        m <- length(h)
        out <- vector("list", m)
        for (i in seq_along(h)) {
            out[[i]]$x <- h[i] + a[i] * cos(theta) * cos(phi[i]) -
                b[i] * sin(theta) * sin(phi[i])
            out[[i]]$y <- k[i] + b[i] * sin(theta) * cos(phi[i]) +
                a[i] * cos(theta) * sin(phi[i])
        }
        out
    }

    e <- eulerr_ellipse(h, k, a, b, phi, n)
    names(e) = colnames(object)
    ell_dt = lapply(e, function(ell)data.table::data.table(x = ell$x, y = ell$y))
    ell_dt = data.table::rbindlist(ell_dt, use.names = T, idcol = "group")
    ggplot(ell_dt) + geom_polygon(aes(x = x, y = y, fill = group), alpha = .3) + theme_void()


    if(is.null(col_scale)){
        col_scale = safeBrew(ncol(object), "Dark2")
    }

    if (fill_circles) {
        p = ggplot(ell_dt, aes(x = x, y = y, fill = group, col = group, size = 2, alpha = alpha))
    } else {
        p = ggplot(ell_dt, aes(x = x, y = y, fill = NA, col = group, size = 2))
    }
    p = p + geom_polygon() + labs(fill = "", color = "") + scale_size_identity() + scale_shape_identity() + scale_alpha_identity() +
        scale_fill_manual(values = col_scale) + scale_color_manual(values = col_scale) + theme_minimal() + theme(plot.background = element_blank(),
                                                                                                                 axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "top") +
        guides(fill = guide_legend(override.aes = list(shape = 21))) + coord_fixed()  #+
    p
}

#' binary heatmap indicating membership.
#' heatmap is sorted by column left to right.  change column order to reveal patterns
#' @export
#' @param object passed to ssvMakeMembTable
#' @param raster_approximation instead of standard plot, write temporary raster image and redraw that as plot background.
#' @param raster_width_min raster width will be minimun multiple of number of columns over this number
#' @param raster_height_min raster height will be minimun multiple of number of rows over this number
#' @rawNamespace import(data.table, except = shift)
#' @import png
#' @import ggplot2
#' @importFrom grid rasterGrob
ssvFeatureBinaryHeatmap = function(object, raster_approximation = TRUE, raster_width_min = 1000, raster_height_min = 1000) {
    groups = bool = value = NULL#declare binding for data.table
    object = ssvMakeMembTable(object)
    mat = ifelse(as.matrix(object), 0, 1)
    for (i in rev(1:ncol(mat))) {
        mat = mat[order(mat[, i], decreasing = TRUE), , drop = FALSE]
    }
    hdt = data.table::as.data.table(cbind(mat, row = 1:nrow(mat)))
    hdt = data.table::melt(hdt, id.vars = "row", variable.name = "groups")

    if (raster_approximation) {
        dmat = as.matrix(data.table::dcast(hdt, value.var = "value", formula = rev(row) ~ groups))[, 1 + seq_len(ncol(object)), drop = FALSE]
        low_v = 0.9372549
        dmat = ifelse(dmat > low_v, low_v, dmat)
        if (raster_width_min >= ncol(dmat)) {
            # expand cols as necessary to meet target
            raster_width_multiplier = ceiling(raster_width_min/ncol(dmat))
            comp_col = rep(1:ncol(dmat), each = raster_width_multiplier)
        } else {
            # evenly spaced sample down to target
            comp_col = round(1:raster_width_min/raster_width_min * ncol(dmat))
        }
        if (raster_height_min >= nrow(dmat)) {
            # expand rows as necessary to meet target
            raster_height_multiplier = ceiling(raster_height_min/nrow(dmat))
            comp_row = rep(1:nrow(dmat), each = raster_height_multiplier)
        } else {
            # evenly spaced sample down to target
            comp_row = round(1:raster_height_min/raster_height_min * nrow(dmat))
        }
        # dmat = (dmat * -.95 + .95)
        comp_mat = dmat[comp_row, comp_col]
        png::writePNG(comp_mat, target = "tmp.png")
        png_grob = grid::rasterGrob(png::readPNG("tmp.png"),
                                    width = unit(1, "npc"),
                                    height = unit(1, "npc"),
                                    interpolate = FALSE)
        file.remove("tmp.png")
        p = ggplot(hdt) +
            geom_tile(aes(x = groups, y = row, fill = NA, col = NULL)) +
            scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "#EFEFEF")) +
            scale_alpha(0) +
            labs(fill = "binding", y = "", title = "Marks per region clustering") +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                plot.background = element_blank(),
                panel.background = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text = element_text(size = rel(1.5)),
                legend.key.size = unit(0.12, units = "npc"),
                legend.text = element_text(size = rel(2)),
                legend.title = element_text(size = rel(3)),
                axis.title = element_blank()
            ) +
            annotation_custom(png_grob, xmin = 0.5, xmax = ncol(dmat) + 0.5,
                              ymin = 0.5, ymax = nrow(dmat) + 0.5)
    } else {
        hdt[, `:=`(bool, value == 1)]
        p = ggplot(hdt) +
            geom_tile(aes(x = groups, y = row, fill = bool, col = NULL)) +
            scale_fill_manual(values = c(`FALSE` = "black", `TRUE` = "#EFEFEF")) +
            labs(fill = "binding", y = "", title = "Marks per region clustering") +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                plot.background = element_blank(),
                panel.background = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text = element_text(size = rel(1.5)),
                legend.key.size = unit(0.12, units = "npc"),
                legend.text = element_text(size = rel(2)),
                legend.title = element_text(size = rel(3)),
                axis.title = element_blank()
            )
    }
    return(p)
}