
ggCircle <- function(r, xc, yc, color="black", fill=NA, ...) {
  x <- xc + r*cos(seq(0, pi, length.out=100))
  ymax <- yc + r*sin(seq(0, pi, length.out=100))
  ymin <- yc + r*sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}

#ggplot version of vennDiagram from limma package
ggVenn = function (object, names = NULL,
                   labels_size = 10,
                   counts_size = 10,
                   show_outside_count = F,
                   lwd = 1,
                   circle.col = NULL,
                   fill_circles = T,
                   counts.col = NULL){
  library(ggplot2)
  if(class(object) == "GRanges"){
    object = elementMetadata(object)
  }
  assert_that(is.matrix(object) | is.data.frame(object) | class(object) == "DataFrame")
  all(apply(object, 2, class) == "logical")
  object <- vennCounts(object)
  z <- object[, "Counts"]
  nsets <- ncol(object) - 1
  if (nsets > 3)
    stop("Can't plot Venn diagram for more than 3 sets")
  if (is.null(names))
    names <- factor(colnames(object)[1:nsets], levels = colnames(object)[1:nsets])
  if (is.null(circle.col)) {
    circle.col = suppressWarnings(RColorBrewer::brewer.pal(nsets, "Dark2"))
    # fill.col = circle.col
  }
  if (length(circle.col) < nsets)
    circle.col <- rep(circle.col, length.out = nsets)
  if (is.null(counts.col))
    counts.col <- par("col")
  col_scale = circle.col
  names(col_scale) = names
  blank = rep(NA, length(names))
  p = ggplot(data.frame(blank), aes(fill = names, col = names, x = NA, y = NA, size = -1)) +
    labs(fill = "", color = "") +
    scale_size_identity() +
    scale_shape_identity() +
    scale_fill_manual(values = col_scale) +
    scale_color_manual(values = col_scale) +
    geom_point()+
    theme_minimal() +
    theme(plot.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "top") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 10 )))
  if (nsets <= 3) {
    p = p + coord_fixed()
    xcentres <- switch(nsets, 0, c(-1, 1), c(-1, 1, 0))
    ycentres <- switch(nsets, 0, c(0, 0), c(1, 1, -2)/sqrt(3))
    r <- 1.5
    xtext <- switch(nsets,
                    -1.2,
                    c(-1.2, 1.2),
                    c(-1.2, 1.2, 0))
    ytext <- switch(nsets,
                    1.8,
                    c(1.8, 1.8),
                    c(2.4, 2.4, -3))
    for (circle in 1:nsets) {
      if (!(fill_circles)){
        p = p + ggCircle(r = r, xc = xcentres[circle],
                         yc = ycentres[circle],
                         col = circle.col[circle],
                         size = lwd)
      }else{
        RGB <- col2rgb(circle.col[circle])/255
        ALPHA <- 0.06
        RGB.ALP <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1],
                       alpha = ALPHA)
        p = p + ggCircle(r = r, xc = xcentres[circle],
                         yc = ycentres[circle],
                         fill = alpha(circle.col[circle], ALPHA),
                         col = circle.col[circle],
                         size = lwd)
      }
      # p = p + annotate("text", x = xtext[circle], y = ytext[circle], label = names[circle], size = labels_size)
    }

    # switch(nsets, rect(-3, -2.5, 3, 2.5), rect(-3, -2.5,
    #                                            3, 2.5), rect(-3, -3.5, 3, 3.3))
    showCounts <- switch(nsets,
                         function(p, counts, text_size, adj, col) {
                           p = p + annotate("text", x = 0, y = 0, label = counts[2], size = counts_size, col = col)
                           if(show_outside_count){
                             p = p + annotate("text", x = 2.3, y = -2.1, label = counts[1], size = counts_size, col = col)

                           }
                           p
                         }, function(p, counts, text_size, adj, col) {
                           p = p + annotate("text",
                                            x = c(1.5, -1.5, 0),
                                            y = c(0.1, 0.1, 0.1),
                                            label = counts[2:4],
                                            size = counts_size, col = col)
                           if(show_outside_count){
                             p = p + annotate("text", x = 2.3, y = -2.1, label = counts[1], size = counts_size, col = col)
                           }
                           p
                         }, function(p, counts, text_size, adj, col) {
                           p = p + annotate("text",
                                            x = c(0, 1.5, .75, -1.5, -.75, 0, 0),
                                            y = c(-1.7, 1, -.35, 1, -.35, .9, 0),
                                            label = counts[2:8],
                                            size = counts_size, col = col)
                           if(show_outside_count){
                             p = p + annotate("text", x = 2.5, y = -3, label = counts[1], size = counts_size, col = col)
                           }
                           p
                         })
    p = showCounts(p, counts = z, text_size = counts_size, adj = c(.5,.5), col = counts.col[1])

    #only up to 3 sets supported so far
  }
  p = p + scale_x_discrete(expand = c(.05))+ scale_y_discrete(expand = c(.05))
  return(p)
}

ggBars = function(df){
  df = as.data.frame(df)
  hit_counts = colSums(df)
  hit_counts_df = data.frame(count = hit_counts,
                             group = factor(names(hit_counts), levels = names(hit_counts)))
  p<- ggplot(hit_counts_df, aes(x=group, y=count, fill=group))+
    labs(x = "") +
    geom_bar(width = 1, stat = "identity") +
    guides(fill = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  p = p + annotate("text", y = hit_counts / 2, x = 1:length(hit_counts), label = hit_counts)
  return(p)
}

ggPie = function(df){
  df = as.data.frame(df)
  hit_counts = colSums(df)
  hit_counts_df = data.frame(count = hit_counts,
                             group = factor(names(hit_counts), levels = names(hit_counts)))

  p<- ggplot(hit_counts_df, aes(x="", y=count, fill=(group)))+
    labs(x = "") +
    geom_bar(width = 1, stat = "identity") +
    guides(x = "none") +
    theme(axis.text.x = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank()) +
    coord_polar("y", start = 0) + scale_y_reverse()
  hc = rev(hit_counts) / sum(hit_counts)
  hc = c(0, cumsum(hc)[-length(hc)]) + hc / 2
  p = p + annotate(geom = "text", y = hc * sum(hit_counts), x = 1.1, label = rev(hit_counts_df$count))
  return(p)
}

#from https://gist.github.com/trinker/31edc08d0a4ec4c73935a23040c2f6cb
# p_load(dplyr, venneuler, ggforce, textshape)
#' Try to load a bed-like file and convert it to a GRanges object
#'
#' @param memb A membership table
#' @return ggplot of venneuler results
#' @import venneuler
#' @examples
ggEuler = function(memb, line_width = 2){
  cn = colnames(memb)
  # cn = LETTERS[1:4]
  todo = expand.grid(lapply(1:ncol(memb), function(x)0:1))
  grp_names = apply(todo, 1, function(x){
    x = as.logical(x)
    nam = paste(cn[x], collapse = "&")
    nam
  })
  grp_counts = apply(todo, 1, function(x){
    x = as.logical(x)
    count = sum(apply(memb, 1, function(y){
      all(x == y)
    }))
    count
  })
  names(grp_counts) = grp_names

  # eul = list("cat" = grp_counts[-1], "dog" = grp_counts[-1])
  eul = list("placeholder" = grp_counts[-1])

  ve = venneuler::venneuler(grp_counts[-1])
  df = data.frame(ve$centers, diameters = ve$diameters, labels = ve$labels, stringsAsFactors = FALSE)
  df$r = df$diameters / 2
  col_scale = RColorBrewer::brewer.pal(nrow(df), "Dark2")
  blank = rep(NA, length(df$labels))
  p = ggplot(data.frame(blank), aes(fill = df$labels, col = df$labels, x = NA, y = NA, size = -1)) +
    labs(fill = "", color = "") +
    scale_size_identity() +
    scale_shape_identity() +
    scale_fill_manual(values = col_scale) +
    scale_color_manual(values = col_scale) +
    geom_point()+
    theme_minimal() +
    theme(plot.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "top") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 10 ))) +
    coord_fixed() +
    scale_x_discrete(expand = c(.05))+ scale_y_discrete(expand = c(.05))
  for(i in seq_len(nrow(df))){
    circ = df[i,]
    p = p + ggCircle(xc = circ$x+.5, yc = circ$y+.5, r = circ$r, color = col_scale[i], size = 3)
  }

  add_fill = T
  if(add_fill){
    q = ggplot(df, aes(x0 = x, y0 = y, r = r, fill = labels, col = labels, size = 2, alpha = .08))
  }else{
    q = ggplot(df, aes(x0 = x, y0 = y, r = r, col = labels, size = 2))
  }
  q = q +
    ggforce::geom_circle() +
    labs(fill = "", color = "") +
    scale_size_identity() +
    scale_shape_identity() +
    scale_alpha_identity() +
    scale_fill_manual(values = col_scale) +
    scale_color_manual(values = col_scale) +
    theme_minimal() +
    theme(plot.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "top") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 10 ))) +
    coord_fixed() #+
  q
}
