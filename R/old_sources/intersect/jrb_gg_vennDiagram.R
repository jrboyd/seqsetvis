#tweaked vennDiagram from limma package
vennCounts = limma::vennCounts

gg_circle <- function(r, xc, yc, color="black", fill=NA, ...) {
  x <- xc + r*cos(seq(0, pi, length.out=100))
  ymax <- yc + r*sin(seq(0, pi, length.out=100))
  ymin <- yc + r*sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}

gg_venn = function (object, names = NULL, 
                    labels_size = 10,
                    counts_size = 10,
                    show_outside_count = F,
                    lwd = 1, 
                    circle.col = NULL, 
                    counts.col = NULL){
  require(ggplot2)
  object <- vennCounts(object)
  z <- object[, "Counts"]
  nsets <- ncol(object) - 1
  if (nsets > 3) 
    stop("Can't plot Venn diagram for more than 3 sets")
  if (is.null(names)) 
    names <- colnames(object)[1:nsets]
  FILL.COL <- TRUE
  if (is.null(circle.col)) {
    circle.col <- par("col")
    FILL.COL <- FALSE
  }
  if (length(circle.col) < nsets) 
    circle.col <- rep(circle.col, length.out = nsets)
  if (is.null(counts.col)) 
    counts.col <- par("col")
  col_scale = circle.col
  names(col_scale) = names
  blank = rep(NA, length(names))
  p = ggplot(data.frame(blank), aes(fill = names, col = names, x = NA, y = NA, size = -1)) +
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
    p = p + coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4))
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
      if (!FILL.COL) 
        p = p + gg_circle(r = r, xc = xcentres[circle], 
                          yc = ycentres[circle], 
                          col = circle.col[circle], 
                          size = lwd)
      if (FILL.COL) {
        RGB <- col2rgb(circle.col[circle])/255
        ALPHA <- 0.06
        RGB.ALP <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1], 
                       alpha = ALPHA)
        p = p + gg_circle(r = r, xc = xcentres[circle], 
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
  return(p)
}

