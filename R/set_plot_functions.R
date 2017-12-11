col2hex = function(color_name) {
    rgb(t(col2rgb(color_name))/255)
}

set_list2memb = function(set_list) {
    if (is.null(names(set_list))) {
        names(set_list) = paste0("set_", LETTERS[seq_along(set_list)])
    }
    rn = unique(unlist(set_list))
    cn = names(set_list)
    memb = matrix(F, nrow = length(rn), ncol = length(cn))
    rownames(memb) = rn
    colnames(memb) = cn
    for (column in cn) {
        memb[set_list[[column]], column] = T
    }
    return(memb)
}

input2membership_matrix = function(object) {
    if (class(object) == "GRanges") {
        object = elementMetadata(object)
    }
    if (class(object) == "DataFrame") {
        object = as.data.frame(object)
    }
    
    if (class(object) == "list") {
        if (all(sapply(object, class) == "factor")) {
            object = lapply(object, as.character)
        }
        if (all(sapply(object, class) == "character")) {
            object = set_list2memb(object)
        } else {
            stop(paste("can't handle list of non-character classes: ", paste(sapply(object, class), 
                collapse = ", ")))
        }
        
    }
    if (is.matrix(object)) {
        object = as.data.frame(object)
    }
    if (is.null(colnames(object))) {
        colnames(object) = paste0("set_", LETTERS[seq_len(ncol(object))])
    }
    assertthat::assert_that(class(object) == "data.frame")
    return(object)
}

# ggplot version of vennDiagram from limma package
ggVenn = function(object, group_names = NULL, counts_txt_size = 5, counts_as_labels = F, show_outside_count = F, 
    lwd = 3, circle.col = NULL, fill_circles = T, fill_alpha = ifelse(fill_circles, 0.5, 0), counts.col = NULL) {
    library(ggplot2)
    object = input2membership_matrix(object)
    # assert_that(is.matrix(object) | is.data.frame(object) | class(object) == 'DataFrame' |
    # is.list(object))
    all(apply(object, 2, class) == "logical")
    object <- limma::vennCounts(object)
    set_counts <- object[, "Counts"]
    nsets <- ncol(object) - 1
    if (nsets > 3) 
        stop("Can't plot Venn diagram for more than 3 sets")
    if (is.null(group_names)) 
        group_names <- factor(colnames(object)[1:nsets], levels = colnames(object)[1:nsets])
    if (is.null(circle.col)) {
        circle.col = suppressWarnings(RColorBrewer::brewer.pal(8, "Dark2"))[seq_len(nsets)]
        # fill.col = circle.col
    } else {
        not_hex = substr(circle.col, 0, 1) != "#"
        circle.col[not_hex] = col2hex(circle.col[not_hex])
    }
    if (length(circle.col) < nsets) 
        circle.col <- rep(circle.col, length.out = nsets)
    if (is.null(counts.col)) 
        counts.col <- par("col")
    col_scale = circle.col
    names(col_scale) = group_names
    ahex = substr(rgb(red = 1, blue = 1, green = 1, alpha = fill_alpha), start = 8, stop = 9)
    fill_scale = paste0(col_scale, ahex)
    names(fill_scale) = names(col_scale)
    
    p = ggplot() + labs(fill = "", color = "") + scale_color_manual(values = col_scale) + scale_size_identity() + 
        scale_fill_manual(values = fill_scale) + theme_minimal() + theme(plot.background = element_blank(), 
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        legend.position = "top") + guides(fill = guide_legend(override.aes = list(shape = 21, size = 10)))
    p = p + coord_fixed()
    xcentres <- switch(nsets, 0, c(-1, 1), c(-1, 1, 0))
    ycentres <- switch(nsets, 0, c(0, 0), c(1, 1, -2)/sqrt(3))
    r <- 1.5
    xtext <- switch(nsets, -1.2, c(-1.2, 1.2), c(-1.2, 1.2, 0))
    ytext <- switch(nsets, 1.8, c(1.8, 1.8), c(2.4, 2.4, -3))
    
    df_circles = data.frame(xcentres, ycentres, r, size = lwd, col = circle.col, group = group_names)
    
    
    p = p + ggforce::geom_circle(data = df_circles, mapping = aes(x0 = xcentres, y0 = ycentres, r = r, 
        size = size, col = group, fill = group))
    showCounts <- switch(nsets, function(counts, text_size, col) {
        df_text = data.frame(x = 0, y = 0, label = counts[2], size = counts_txt_size, col = col)
        if (show_outside_count) {
            df_text = rbind(df_text, data.frame(x = 2.3, y = -2.1, label = counts[1], size = counts_txt_size, 
                col = col))
        }
        df_text
    }, function(counts, text_size, col) {
        df_text = data.frame(x = c(1.5, -1.5, 0), y = c(0.1, 0.1, 0.1), label = counts[2:4], size = counts_txt_size, 
            col = col)
        if (show_outside_count) {
            df_text = rbind(df_text, data.frame(x = 2.3, y = -2.1, label = counts[1], size = counts_txt_size, 
                col = col))
        }
        df_text
    }, function(counts, text_size, col) {
        df_text = data.frame(x = c(0, 1.5, 0.75, -1.5, -0.75, 0, 0), y = c(-1.7, 1, -0.35, 1, -0.35, 
            0.9, 0), label = counts[2:8], size = counts_txt_size, col = col)
        if (show_outside_count) {
            df_text = rbind(df_text, data.frame(x = 2.5, y = -3, label = counts[1], size = counts_txt_size, 
                col = col))
        }
        df_text
    })
    
    counts_method = ifelse(counts_as_labels, geom_label, geom_text)
    df_text = showCounts(counts = set_counts, text_size = counts_txt_size, col = counts.col[1])
    
    p = p + counts_method(data = df_text, mapping = aes(x = x, y = y, label = label, size = counts_txt_size))
    return(p)
}
ggBars = function(df) {
    df = input2membership_matrix(df)
    hit_counts = colSums(df)
    hit_counts_df = data.frame(count = hit_counts, group = factor(names(hit_counts), levels = names(hit_counts)))
    p <- ggplot(hit_counts_df, aes(x = group, y = count, fill = group)) + labs(x = "") + geom_bar(width = 1, 
        stat = "identity") + guides(fill = "none") + theme_bw() + theme(axis.text.x = element_text(angle = 90, 
        hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    p = p + annotate("text", y = hit_counts/2, x = 1:length(hit_counts), label = hit_counts)
    return(p)
}

ggPie = function(df) {
    df = input2membership_matrix(df)
    hit_counts = colSums(df)
    hit_counts_df = data.frame(count = hit_counts, group = factor(names(hit_counts), levels = names(hit_counts)))
    
    p <- ggplot(hit_counts_df, aes(x = "", y = count, fill = (group))) + labs(x = "") + geom_bar(width = 1, 
        stat = "identity") + guides(x = "none") + theme(axis.text.x = element_blank(), panel.background = element_blank(), 
        axis.ticks = element_blank()) + coord_polar("y", start = 0) + scale_y_reverse()
    hc = rev(hit_counts)/sum(hit_counts)
    hc = c(0, cumsum(hc)[-length(hc)]) + hc/2
    p = p + annotate(geom = "text", y = hc * sum(hit_counts), x = 1.1, label = rev(hit_counts_df$count))
    return(p)
}

# from https://gist.github.com/trinker/31edc08d0a4ec4c73935a23040c2f6cb p_load(dplyr, venneuler,
# ggforce, textshape)
#' Try to load a bed-like file and convert it to a GRanges object
#'
#' @param memb A membership table
#' @return ggplot of venneuler results
#' @import venneuler
#' @examples
ggEuler = function(memb, line_width = 2) {
    memb = input2membership_matrix(memb)
    cn = colnames(memb)
    todo = expand.grid(lapply(1:ncol(memb), function(x) 0:1))
    grp_names = apply(todo, 1, function(x) {
        x = as.logical(x)
        nam = paste(cn[x], collapse = "&")
        nam
    })
    grp_counts = apply(todo, 1, function(x) {
        x = as.logical(x)
        count = sum(apply(memb, 1, function(y) {
            all(x == y)
        }))
        count
    })
    names(grp_counts) = grp_names
    
    eul = list(placeholder = grp_counts[-1])
    
    ve = venneuler::venneuler(grp_counts[-1])
    df = data.frame(ve$centers, diameters = ve$diameters, labels = ve$labels, stringsAsFactors = FALSE)
    df$r = df$diameters/2
    col_scale = RColorBrewer::brewer.pal(nrow(df), "Dark2")
    blank = rep(NA, length(df$labels))
    
    add_fill = T
    if (add_fill) {
        p = ggplot(df, aes(x0 = x, y0 = y, r = r, fill = labels, col = labels, size = 2, alpha = 0.08))
    } else {
        p = ggplot(df, aes(x0 = x, y0 = y, r = r, col = labels, size = 2))
    }
    p = p + ggforce::geom_circle() + labs(fill = "", color = "") + scale_size_identity() + scale_shape_identity() + 
        scale_alpha_identity() + scale_fill_manual(values = col_scale) + scale_color_manual(values = col_scale) + 
        theme_minimal() + theme(plot.background = element_blank(), axis.title = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "top") + guides(fill = guide_legend(override.aes = list(shape = 21, 
        size = 10))) + coord_fixed()  #+
    p
}

ggMembershipHeatmap = function(df, raster_approximation = T, raster_width_min = 1000, raster_height_min = 1000) {
    df = input2membership_matrix(df)
    mat = ifelse(as.matrix(df), 0, 1)
    for (i in rev(1:ncol(mat))) {
        mat = mat[order(mat[, i], decreasing = T), , drop = F]
    }
    hdt = data.table::as.data.table(cbind(mat, row = 1:nrow(mat)))
    hdt = data.table::melt(hdt, id.vars = "row", variable.name = "groups")
    
    if (raster_approximation) {
        dmat = as.matrix(data.table::dcast(hdt, value.var = "value", formula = rev(row) ~ groups))[, 
            1 + seq_len(ncol(df)), drop = F]
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
        require(grid)
        png_grob = grid::rasterGrob(png::readPNG("tmp.png"), width = unit(1, "npc"), height = unit(1, 
            "npc"), interpolate = F)
        file.remove("tmp.png")
        p = ggplot(hdt) + geom_tile(aes(x = groups, y = row, fill = NA, col = NULL)) + scale_fill_manual(values = c(`TRUE` = "black", 
            `FALSE` = "#EFEFEF")) + scale_alpha(0) + labs(fill = "binding", y = "", title = "Marks per region clustering") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank(), 
                panel.grid.minor.x = element_blank(), plot.background = element_blank(), panel.background = element_blank(), 
                axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
                axis.text = element_text(size = rel(1.5)), legend.key.size = unit(0.12, units = "npc"), 
                legend.text = element_text(size = rel(2)), legend.title = element_text(size = rel(3)), 
                axis.title = element_blank()) + annotation_custom(png_grob, xmin = 0.5, xmax = ncol(dmat) + 
            0.5, ymin = 0.5, ymax = nrow(dmat) + 0.5)
    } else {
        hdt[, `:=`(bool, value == 1)]
        p = ggplot(hdt) + geom_tile(aes(x = groups, y = row, fill = bool, col = NULL)) + scale_fill_manual(values = c(`FALSE` = "black", 
            `TRUE` = "#EFEFEF")) + labs(fill = "binding", y = "", title = "Marks per region clustering") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank(), 
                panel.grid.minor.x = element_blank(), plot.background = element_blank(), panel.background = element_blank(), 
                axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
                axis.text = element_text(size = rel(1.5)), legend.key.size = unit(0.12, units = "npc"), 
                legend.text = element_text(size = rel(2)), legend.title = element_text(size = rel(3)), 
                axis.title = element_blank())
    }
    return(p)
}
