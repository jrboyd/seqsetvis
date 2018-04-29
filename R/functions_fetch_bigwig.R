#' Fetch values from a bigwig appropriate for heatmaps etc.
#'
#' \code{fetchWindowedBigwig} Gets values for each region of the query
#' GRanges (\code{qgr}). Values correspond to the center of each window of
#' size \code{win_size}.  A tidy formatted data.table
#' object is returned suitable for plotting using ggplots.
#' @export
#' @param bw_file The character vector path to bigwig files to read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @param win_method character.  one of c("sample", "summary").  Determines
#' if \code{\link{viewGRangesWinSample_dt}} or
#' \code{\link{viewGRangesWinSummary_dt}} is used to represent each region in
#' qgr.
#' @param summary_FUN function.  only relevant if win_method is "summary".
#' passed to \code{\link{viewGRangesWinSummary_dt}}.
#' @param x0 character, one of c("center", "center_unstranded",
#' "left", "left_unstranded")
#' @param return_data.table logical. If TRUE the internal data.table is
#' returned instead of GRanges.  Default is FALSE.
#' @return A GRanges (or data.table if specified) containing fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved
#' from \code{bw_file}
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
                               win_method = c("sample", "summary")[1],
                               summary_FUN = stats::weighted.mean,
                               x0 = c("left", "left_unstranded", "center",
                                      "center_unstranded")[3],
                               return_data.table = FALSE) {
    stopifnot(is.character(bw_file))
    stopifnot(class(qgr) == "GRanges")
    stopifnot(is.numeric(win_size))

    switch (win_method,
            sample = {
                qgr = prepare_fetch_GRanges(qgr, win_size)
                score_gr = rtracklayer::import.bw(bw_file, which = qgr)
                out = viewGRangesWinSample_dt(score_gr, qgr, win_size, x0 = x0)
            },
            summary = {
                score_gr = rtracklayer::import.bw(bw_file, which = qgr)
                out = viewGRangesWinSummary_dt(score_gr, qgr, win_size,
                                               summary_FUN = summary_FUN,
                                               x0 = x0)
            }
    )

    if(!return_data.table){
        out = GRanges(out)
    }
    return(out)
}


#' Iterates a character vector (ideally named) and calls
#' \code{fetchWindowedBigwig} on each.  Appends grouping variable to each
#' resulting data.table and uses rbindlist to efficiently combine results.
#'
#' \code{fetchWindowedBigwigList} iteratively calls \code{fetchWindowedBigwig}.
#' See \code{\link{fetchWindowedBigwig}} for more info.
#' @export
#' @param file_paths The character vector or list of paths to bigwig files to
#'  read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param unique_names names to use in final data.table to designate source
#' bigwig.
#' @param names_variable The column name where unique_names are stored.
#' Default is 'sample'
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @param win_method character.  one of c("sample", "summary").  Determines
#' if \code{\link{viewGRangesWinSample_dt}} or
#' \code{\link{viewGRangesWinSummary_dt}} is used to represent each region in
#' qgr.
#' @param summary_FUN function.  only relevant if win_method is "summary".
#' passed to \code{\link{viewGRangesWinSummary_dt}}.
#' @param x0 character, one of c("center", "center_unstranded",
#' "left", "left_unstranded")
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
fetchWindowedBigwigList = function(file_paths,
                                   qgr,
                                   unique_names = names(file_paths),
                                   names_variable = "sample",
                                   win_size = 50,
                                   win_method = c("sample", "summary")[1],
                                   summary_FUN = stats::weighted.mean,
                                   x0 = c("left", "left_unstranded", "center",
                                          "center_unstranded")[3],
                                   return_data.table = FALSE) {

    load_bw = function(f, nam, qgr) {
        message("loading ", f, " ...")
        dt = fetchWindowedBigwig(bw_file = f,
                                 qgr = qgr,
                                 win_size = win_size,
                                 win_method = win_method,
                                 summary_FUN = summary_FUN,
                                 x0 = x0,
                                 return_data.table = TRUE)
        dt[[names_variable]] = nam
        message("finished loading ", nam, ".")
        dt
    }

    fetchWindowedSignalList(file_paths = file_paths,
                            qgr = qgr,
                            load_signal = load_bw,
                            unique_names = unique_names,
                            names_variable = names_variable,
                            win_size = win_size,
                            return_data.table = return_data.table)
}



