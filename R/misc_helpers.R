#' converts a valid r color name ("black", "red", "white", etc.) to a hex value
#' @export
#' @param color_name character. one or more r color names.
#' @importFrom grDevices col2rgb rgb
#' @return hex value of colors coded by colors()
#' @examples
#' col2hex(c("red", "green", "blue"))
#' col2hex(c("lightgray", "gray", "darkgray"))
col2hex = function(color_name) {
  grDevices::rgb(t(grDevices::col2rgb(color_name))/255)
}

#' convert a list of sets, each list item should be a character vector
#' denoting items in sets
#' @param set_list a list of character vectors.  default names will be added if
#' missing
#' @return converts list of characters/numeric to membership table matrix
set_list2memb = function(set_list) {
  if (is.null(names(set_list))) {
    names(set_list) = paste0("set_", LETTERS[seq_along(set_list)])
  }
  rn = unique(unlist(set_list))
  cn = names(set_list)
  memb = matrix(FALSE, nrow = length(rn), ncol = length(cn))
  rownames(memb) = rn
  colnames(memb) = cn
  for (column in cn) {
    memb[set_list[[column]], column] = TRUE
  }
  return(memb)
}

#' allows RColorBrew to handle n values less than 3 and greater than 8 without
#' warnings and return expected number of colors.
#' @export
#' @param n integer value of number of colors to make palette for
#' @param pal palette recognized by RColorBrewer
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @return a character vector of hex coded colors o flength n from the color
#' brewer palette pal
#' @examples
#' plot(1:2, rep(0, 2),  col = safeBrew(2, "dark2"), pch = 16, cex = 6)
#' plot(1:12, rep(0, 12),  col = safeBrew(12, "set1"), pch = 16, cex = 6)
#' plot(1:12, rep(0, 12),  col = safeBrew(12, "set2"), pch = 16, cex = 6)
#' plot(1:12, rep(0, 12),  col = safeBrew(12, "set3"), pch = 16, cex = 6)
safeBrew = function(n, pal = "Dark2"){
    if(n < 1) stop("n must be at least 1")
    pal_info = RColorBrewer::brewer.pal.info
    pal_info$brewName = rownames(pal_info)
    rownames(pal_info) = tolower(rownames(pal_info))
    pal = tolower(pal)
    if(!any(pal == rownames(pal_info)))
        stop(paste("Palette", pal, "not a valid RColorBrewer palette, see RColorBrewer::brewer.pal.info"))
    maxColors = pal_info[pal,]$maxcolors
    nbrew = max(n, 3) %>% min(., maxColors)
    RColorBrewer::brewer.pal(n = nbrew, name = pal_info[pal,]$brewName)[(seq_len(n)-1) %% maxColors + 1]
}
