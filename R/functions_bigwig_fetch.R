

#' Fetch values from a bigwig appropriate for heatmaps etc.
#'
#' \code{fetchWindowedBigwig} Gets values for each region of the query GRanges (\code{qgr}).
#' Values correspond to the center of each window of size \code{win_size}.  A tidy formatted data.table
#' object is returned suitable for plotting using ggplots.
#' @export
#' @param bw_file The character vector path to bigwig files to read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @return A tidy formatted data.table containing fetched values.
#'
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved from \code{bw_file}
#'
fetchWindowedBigwig = function(bw_file, qgr, win_size = 50) {
    queryHits = id = x = NULL
    if (!all(width(qgr)%%win_size == 0)) {
        stop(paste("all widths of qgr are not evenly divisible by win_size,", win_size))
    }
    if (!all(width(qgr) == width(qgr)[1])) {
        warning(paste("Widths vary between GRanges.  Before aggregating or co-plotting, recommend one of:
                  1) modifying input qgr such that all GRanges have same width
                  2) filtering output so that the same values of x are present for every region"))
    }
    # suppressWarnings({
    bw_gr = rtracklayer::import.bw(bw_file, which = qgr)
    # })
    windows = slidingWindows(qgr, width = win_size, step = win_size)
    if (is.null(qgr$id)) {
        if (!is.null(names(qgr))) {
            qgr$id = names(qgr)
        } else {
            qgr$id = paste0("region_", seq_along(qgr))
        }
    }
    names(windows) = qgr$id
    # print(object.size(windows), units = 'GB')
    windows = unlist(windows)
    windows$id = names(windows)
    # mid_gr = function(gr) {
    #     start(gr) + floor((width(gr) - 1)/2)
    # }
    # mids = mid_gr(windows)
    # start(windows) = mids
    # end(windows) = mids
    windows = resize(windows, width = 1, fix = "center")
    olaps = suppressWarnings(data.table::as.data.table(findOverlaps(query = windows, subject = bw_gr)))
    # patch up missing/out of bound data with 0
    missing_idx = setdiff(1:length(windows), olaps$queryHits)
    if (length(missing_idx) > 0) {
        olaps = rbind(olaps, data.table::data.table(queryHits = missing_idx, subjectHits = length(bw_gr) + 1))[order(queryHits)]
        bw_gr = c(bw_gr, GRanges(seqnames(bw_gr)[length(bw_gr)], IRanges::IRanges(1, 1), score = 0))
    }
    # set y and output windows = windows[olaps$queryHits]
    windows$y = bw_gr[olaps$subjectHits]$score
    bw_dt = data.table::as.data.table(windows)
    bw_dt[, `:=`(x, start - min(start) + win_size/2), by = id]
    bw_dt[, `:=`(x, x - round(mean(x))), by = id]
    shift = round(win_size/2)
    bw_dt[, `:=`(start, start - shift + 1)]
    bw_dt[, `:=`(end, end + win_size - shift)]
    return(bw_dt)
}

#' Iterates a character vector (ideally named) and calls \code{fetchWindowedBigwig}
#' on each.  Appends grouping variable to each resulting data.table and uses rbindlist to
#' efficiently combine results
#'
#' \code{fetchWindowedBigwigList} iteratively calls \code{fetchWindowedBigwig}.
#' See \code{\link{fetchWindowedBigwig}} for more info.
#' @export
#' @param bw_files The character vector paths to bigwig files to read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param bw_names names to use in final data.table to designate source bigwig
#' @param bw_variable_name The column name where bw_names are stored.
#' Default is 'sample'
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @return A tidy formatted data.table containing fetched values.
#'
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved from \code{bw_file}
fetchWindowedBigwigList = function(bw_files, qgr, bw_names = names(bw_files), bw_variable_name = "sample", win_size = 50) {
    if (is.null(bw_names)) {
        bw_names = basename(bw_files)
    }
    if (any(duplicated(bw_names))) {
        stop(paste0("some bw_names are duplicated:\n", paste(collapse = "\n", unique(bw_names[duplicated(bw_names)]))))
    }

    load_bw = function(nam) {
        print(paste0("loading ", nam, " ..."))
        f = bw_files[nam]
        dt = fetchWindowedBigwig(bw_file = f, win_size = win_size, qgr = qgr)
        dt[[bw_variable_name]] = nam
        print(paste0("finished loading ", nam, "."))
        dt
    }

    bw_list = lapply(names(bw_files), load_bw)
    all_bw_dt = data.table::rbindlist(bw_list)
    return(all_bw_dt)
}



