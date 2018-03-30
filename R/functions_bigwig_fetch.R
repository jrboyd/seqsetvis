# windowViewGRanges_dt = function(score_gr, query_gr, window_size){
#     x = id = NULL
#     windows = slidingWindows(query_gr, width = window_size, step = window_size)
#     if (is.null(query_gr$id)) {
#         if (!is.null(names(query_gr))) {
#             query_gr$id = names(query_gr)
#         } else {
#             query_gr$id = paste0("region_", seq_along(query_gr))
#         }
#     }
#     names(windows) = query_gr$id
#     windows = unlist(windows)
#     windows$id = names(windows)
#     windows = resize(windows, width = 1, fix = "center")
#     olaps = suppressWarnings(data.table::as.data.table(findOverlaps(query = windows, subject = score_gr)))
#     # patch up missing/out of bound data with 0
#     missing_idx = setdiff(seq_along(windows), olaps$queryHits)
#     if (length(missing_idx) > 0) {
#         olaps = rbind(olaps, data.table::data.table(queryHits = missing_idx, subjectHits = length(score_gr) + 1))[order(queryHits)]
#         score_gr = c(score_gr, GRanges(seqnames(score_gr)[length(score_gr)], IRanges::IRanges(1, 1), score = 0))
#     }
#     # set y and output windows = windows[olaps$queryHits]
#     windows$y = score_gr[olaps$subjectHits]$score
#     score_dt = data.table::as.data.table(windows)
#     score_dt[, `:=`(x, start - min(start) + window_size/2), by = id]
#     score_dt[, `:=`(x, x - round(mean(x))), by = id]
#     shift = round(window_size/2)
#     score_dt[, `:=`(start, start - shift + 1)]
#     score_dt[, `:=`(end, end + window_size - shift)]
# }


#' Fetch values from a bigwig appropriate for heatmaps etc.
#'
#' \code{fetchWindowedBigwig_dt} Gets values for each region of the query GRanges (\code{qgr}).
#' Values correspond to the center of each window of size \code{win_size}.  A tidy formatted data.table
#' object is returned suitable for plotting using ggplots.
#' @export
#' @param bw_file The character vector path to bigwig files to read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @return A tidy formatted data.table containing fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved from \code{bw_file}
#' @examples
#' if(Sys.info()['sysname'] != "Windows"){
#' library(GenomicRanges)
#' bw_f = system.file("extdata/test_bigwigs/test_loading.bw",
#'     package = "seqsetvis", mustWork = TRUE)
#' qgr = GRanges("chrTest", IRanges(1, 30))
#' bw_dt = fetchWindowedBigwig_dt(bw_f, qgr, win_size = 10)
#' }
fetchWindowedBigwig_dt = function(bw_file, qgr, win_size = 50) {
    queryHits = id = x = NULL
    stopifnot(is.character(bw_file))
    stopifnot(class(qgr) == "GRanges")
    stopifnot(is.numeric(win_size))
    if (!all(width(qgr)%%win_size == 0)) {
        stop("all widths of qgr are not evenly divisible by win_size, ", win_size)
    }
    if (!all(width(qgr) == width(qgr)[1])) {
        warning("Widths vary between GRanges.  Before aggregating or co-plotting, recommend one of:
            1) modifying input qgr such that all GRanges have same width
            2) filtering output so that the same values of x are present for every region")
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
    windows = unlist(windows)
    windows$id = names(windows)
    windows = resize(windows, width = 1, fix = "center")
    olaps = suppressWarnings(data.table::as.data.table(findOverlaps(query = windows, subject = bw_gr)))
    # patch up missing/out of bound data with 0
    missing_idx = setdiff(seq_along(windows), olaps$queryHits)
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
#' @return A GRanges containing fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved from \code{bw_file}
#' @examples
#' if(Sys.info()['sysname'] != "Windows"){
#' library(GenomicRanges)
#' bw_f = system.file("extdata/test_bigwigs/test_loading.bw",
#'     package = "seqsetvis", mustWork = TRUE)
#' qgr = GRanges("chrTest", IRanges(1, 30))
#' bw_gr = fetchWindowedBigwig(bw_f, qgr, win_size = 10)
#' }
fetchWindowedBigwig = function(bw_file, qgr, win_size = 50) {
    GRanges(fetchWindowedBigwig_dt(bw_file, qgr, win_size))
}

#' Iterates a character vector (ideally named) and calls \code{fetchWindowedBigwig}
#' on each.  Appends grouping variable to each resulting data.table and uses rbindlist to
#' efficiently combine results
#'
#' \code{fetchWindowedBigwigList} iteratively calls \code{fetchWindowedBigwig}.
#' See \code{\link{fetchWindowedBigwig}} for more info.
#' @export
#' @param bw_files The character vector or list of paths to bigwig files to
#'  read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param bw_names names to use in final data.table to designate source bigwig
#' @param bw_variable_name The column name where bw_names are stored.
#' Default is 'sample'
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @return A GRanges containing fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be
#' retrieved from \code{bw_file}
#' @examples
#' if(Sys.info()['sysname'] != "Windows"){
#' library(GenomicRanges)
#' bw_f = system.file("extdata/test_bigwigs/test_loading.bw",
#'     package = "seqsetvis", mustWork = TRUE)
#' bw_files = c("a" = bw_f, "b" = bw_f)
#' qgr = GRanges("chrTest", IRanges(1, 30))
#' bw_gr = fetchWindowedBigwigList(bw_files, qgr, win_size = 10)
#' bw_gr2 = fetchWindowedBigwigList(as.list(bw_files), qgr, win_size = 10)
#' }
fetchWindowedBigwigList = function(bw_files, qgr, bw_names = names(bw_files),
                                   bw_variable_name = "sample", win_size = 50) {
    all_bw_gr = GRanges(fetchWindowedBigwigList_dt(bw_files, qgr, bw_names, bw_variable_name, win_size))
    return(all_bw_gr)
}


#' Iterates a character vector (ideally named) and calls \code{fetchWindowedBigwig}
#' on each.  Appends grouping variable to each resulting data.table and uses rbindlist to
#' efficiently combine results
#'
#' \code{fetchWindowedBigwigList} iteratively calls \code{fetchWindowedBigwig}.
#' See \code{\link{fetchWindowedBigwig}} for more info.
#' @export
#' @param bw_files The character vector or list of paths to bigwig files to
#'  read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param bw_names names to use in final data.table to designate source bigwig
#' @param bw_variable_name The column name where bw_names are stored.
#' Default is 'sample'
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @return A tidy formatted data.table containing fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be
#' retrieved from \code{bw_file}
#' @examples
#' if(Sys.info()['sysname'] != "Windows"){
#' library(GenomicRanges)
#' bw_f = system.file("extdata/test_bigwigs/test_loading.bw",
#'     package = "seqsetvis", mustWork = TRUE)
#' bw_files = c("a" = bw_f, "b" = bw_f)
#' qgr = GRanges("chrTest", IRanges(1, 30))
#' bw_dt = fetchWindowedBigwigList(bw_files, qgr, win_size = 10)
#' bw_dt2 = fetchWindowedBigwigList(as.list(bw_files), qgr, win_size = 10)
#' }
fetchWindowedBigwigList_dt = function(bw_files, qgr, bw_names = names(bw_files),
                                   bw_variable_name = "sample", win_size = 50) {
    if(is.list(bw_files)){
        bw_files = unlist(bw_files)
    }
    if (is.null(bw_names)) {
        bw_names = basename(bw_files)
    }
    names(bw_files) = bw_names
    stopifnot(is.character(bw_files))
    stopifnot(class(qgr) == "GRanges")
    stopifnot(is.character(bw_names))
    stopifnot(is.character(bw_variable_name))
    stopifnot(is.numeric(win_size))

    if (any(duplicated(bw_names))) {
        stop("some bw_names are duplicated:\n",
                    paste(collapse = "\n", unique(bw_names[duplicated(bw_names)])))
    }

    load_bw = function(nam) {
        message("loading ", nam, " ...")
        f = bw_files[nam]
        dt = fetchWindowedBigwig_dt(bw_file = f, win_size = win_size, qgr = qgr)
        dt[[bw_variable_name]] = nam
        message("finished loading ", nam, ".")
        dt
    }

    bw_list = lapply(names(bw_files), load_bw)
    all_bw_dt = data.table::rbindlist(bw_list)
    return(all_bw_dt)
}


