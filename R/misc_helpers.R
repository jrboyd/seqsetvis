#' converts a valid r color name ("black", "red", "white", etc.) to a hex value
#' @export
#' @param color_name character. one or more r color names.
#' @importFrom grDevices col2rgb rgb
#' @return hex value of colors coded by colors()
#' @examples
#' col2hex(c("red", "green", "blue"))
#' col2hex(c("lightgray", "gray", "darkgray"))
col2hex = function(color_name) {
    stopifnot(is.character(color_name))
    grDevices::rgb(t(grDevices::col2rgb(color_name))/255)
}

#' convert a list of sets, each list item should be a character vector
#' denoting items in sets
#' @param set_list a list of character vectors.  default names will be added if
#' missing
#' @return converts list of characters/numeric to membership table matrix
set_list2memb = function(set_list) {
    stopifnot(is.list(set_list))
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
    stopifnot(is.numeric(n))
    stopifnot(is.character(pal))
    if(n < 1) stop("n must be at least 1")
    pal_info = RColorBrewer::brewer.pal.info
    pal_info$brewName = rownames(pal_info)
    rownames(pal_info) = tolower(rownames(pal_info))
    pal = tolower(pal)
    if(!any(pal == rownames(pal_info)))
        stop("Palette", pal, "not a valid RColorBrewer palette, see RColorBrewer::brewer.pal.info")
    maxColors = pal_info[pal,]$maxcolors
    nbrew = min(max(n, 3), maxColors)
    RColorBrewer::brewer.pal(n = nbrew, name = pal_info[pal,]$brewName)[(seq_len(n)-1) %% maxColors + 1]
}

# x: the vector n: the number of samples centered: if FALSE, then average current sample and previous (n-1) samples if TRUE, then average symmetrically in past and
# future. (If n is even, use one more sample from future.)
# from http://www.cookbook-r.com/Manipulating_data/Calculating_a_moving_average/
movingAverage <- function(x, n = 1, centered = T) {

    if (centered) {
        before <- floor((n - 1)/2)
        after <- ceiling((n - 1)/2)
    } else {
        before <- n - 1
        after <- 0
    }

    # Track the sum and count of number of non-NA items
    s <- rep(0, length(x))
    count <- rep(0, length(x))

    # Add the centered data
    new <- x
    # Add to count list wherever there isn't a
    count <- count + (!is.na(new))
    # Now replace NA_s with 0_s and add to total
    new[is.na(new)] <- 0
    s <- s + new

    # Add the data from before
    i <- 1
    while (i <= before) {
        # This is the vector with offset values to add
        new <- c(rep(NA, i), x[1:(length(x) - i)])

        count <- count + (!is.na(new))
        new[is.na(new)] <- 0
        s <- s + new

        i <- i + 1
    }

    # Add the data from after
    i <- 1
    while (i <= after) {
        # This is the vector with offset values to add
        new <- c(x[(i + 1):length(x)], rep(NA, i))

        count <- count + (!is.na(new))
        new[is.na(new)] <- 0
        s <- s + new

        i <- i + 1
    }

    # return sum divided by count
    s/count
}
