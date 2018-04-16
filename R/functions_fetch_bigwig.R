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
#' @param return_data.table logical. If TRUE the internal data.table is
#' returned instead of GRanges.  Default is FALSE.
#' @return A GRanges (or data.table if specified) containing fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved from \code{bw_file}
#' @examples
#' if(Sys.info()['sysname'] != "Windows"){
#' library(GenomicRanges)
#' bw_f = system.file("extdata/test_loading.bw",
#'     package = "seqsetvis", mustWork = TRUE)
#' qgr = GRanges("chrTest", IRanges(1, 30))
#' bw_gr = fetchWindowedBigwig(bw_f, qgr, win_size = 10)
#'
#' bw_dt = fetchWindowedBigwig(bw_f, qgr, win_size = 10,
#'     return_data.table = TRUE)
#' }
fetchWindowedBigwig = function(bw_file,
                               qgr,
                               win_size = 50,
                               return_data.table = FALSE) {
    stopifnot(is.character(bw_file))
    stopifnot(class(qgr) == "GRanges")
    stopifnot(is.numeric(win_size))
    qgr = prepare_fetch_GRanges(qgr = qgr, win_size = win_size, target_size = NULL)
    score_gr = rtracklayer::import.bw(bw_file, which = qgr)
    out = viewGRangesWindowed_dt(score_gr, qgr, win_size)
    if(!return_data.table){
        out = GRanges(out)
    }
    return(out)
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
#' @param return_data.table logical. If TRUE the internal data.table is
#' returned instead of GRanges.  Default is FALSE.
#' @return A tidy formatted GRanges (or data.table if specified) containing
#' fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be
#' retrieved from \code{bw_file}
#' @examples
#' if(Sys.info()['sysname'] != "Windows"){
#' library(GenomicRanges)
#' bw_f = system.file("extdata/test_loading.bw",
#'     package = "seqsetvis", mustWork = TRUE)
#' bw_files = c("a" = bw_f, "b" = bw_f)
#' qgr = GRanges("chrTest", IRanges(1, 30))
#' bw_gr = fetchWindowedBigwigList(bw_files, qgr, win_size = 10)
#' bw_gr2 = fetchWindowedBigwigList(as.list(bw_files), qgr, win_size = 10)
#'
#' bw_dt = fetchWindowedBigwigList(bw_files, qgr, win_size = 10,
#'     return_data.table = TRUE)
#' }
fetchWindowedBigwigList = function(bw_files,
                                   qgr,
                                   bw_names = names(bw_files),
                                   bw_variable_name = "sample",
                                   win_size = 50,
                                   return_data.table = FALSE) {
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
    qgr = prepare_fetch_GRanges(qgr = qgr, win_size = win_size, target_size = NULL)
    load_bw = function(nam) {
        message("loading ", nam, " ...")
        f = bw_files[nam]
        dt = fetchWindowedBigwig(bw_file = f,
                                 win_size = win_size,
                                 qgr = qgr,
                                 return_data.table = TRUE)
        dt[[bw_variable_name]] = nam
        message("finished loading ", nam, ".")
        dt
    }
    bw_list = lapply(names(bw_files), load_bw)
    out = data.table::rbindlist(bw_list)
    if(!return_data.table){
        out = GRanges(out)
    }
    return(out)
}



