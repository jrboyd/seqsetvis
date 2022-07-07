#' Iterates a character vector (ideally named) and calls
#' \code{ssvFetchBigwig.single} on each.  Appends grouping variable to each
#' resulting data.table and uses rbindlist to efficiently combine results.
#'
#' \code{ssvFetchBigwig} iteratively calls \code{fetchWindowedBigwig.single}.
#' See \code{\link{ssvFetchBigwig.single}} for more info.
#' @export
#' @param file_paths character vector of file_paths to load from. Alternatively,
#' file_paths can be a data.frame or data.table whose first column is a
#' character vector of paths and additial columns will be used as metadata.
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
#' @param fragLens never used by ssvFetchBigwig. Ignore.
#' @param anchor character, one of c("center", "center_unstranded",
#' "left", "left_unstranded")
#' @param return_data.table logical. If TRUE the internal data.table is
#' returned instead of GRanges.  Default is FALSE.
#' @param n_cores integer number of cores to use.
#' Uses mc.cores option if not supplied.
#' @param n_region_splits integer number of splits to apply to qgr. The query
#'   GRanges will be split into this many roughly equal parts for increased
#'   parallelization. Default is 1, no split.
#' @param force_skip_centerFix boolean, if TRUE all query ranges will be
#' used "as is".  This is already the case by default if win_method == "summary"
#' but may have applications where win_method == "sample".
#' @return A tidy formatted GRanges (or data.table if specified) containing
#' fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
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
#' bw_gr = ssvFetchBigwig(bw_files, qgr, win_size = 10)
#' bw_gr2 = ssvFetchBigwig(as.list(bw_files), qgr, win_size = 10)
#'
#' bw_dt = ssvFetchBigwig(bw_files, qgr, win_size = 10,
#'     return_data.table = TRUE)
#' }
ssvFetchBigwig = function(file_paths,
                          qgr,
                          unique_names = NULL,
                          names_variable = "sample",
                          win_size = 50,
                          win_method = c("sample", "summary")[1],
                          summary_FUN = stats::weighted.mean,
                          fragLens = "not_used",
                          anchor = c("left", "left_unstranded", "center",
                                     "center_unstranded")[3],
                          return_data.table = FALSE,
                          n_cores = getOption("mc.cores", 1),
                          n_region_splits = 1,
                          force_skip_centerFix = FALSE) {

    load_bw = function(f, nam, qgr) {
        message("loading ", f, " ...")
        dt = ssvFetchBigwig.single(bw_file = f,
                                   qgr = qgr,
                                   win_size = win_size,
                                   win_method = win_method,
                                   summary_FUN = summary_FUN,
                                   anchor = anchor,
                                   return_data.table = TRUE,
                                   force_skip_centerFix = force_skip_centerFix)
        # dt[[names_variable]] = nam
        message("finished loading ", nam, ".")
        dt
    }

    ssvFetchSignal(file_paths = file_paths,
                   qgr = qgr,
                   unique_names = unique_names,
                   names_variable = names_variable,
                   win_size = win_size,
                   win_method = win_method,
                   return_data.table = return_data.table,
                   load_signal = load_bw,
                   n_cores = n_cores,
                   n_region_splits = n_region_splits,
                   force_skip_centerFix = force_skip_centerFix)
}

#' Fetch values from a bigwig appropriate for heatmaps etc.
#'
#' \code{ssvFetchBigwig.single} Gets values for each region of the query
#' GRanges (\code{qgr}). Values correspond to the center of each window of
#' size \code{win_size}.  A tidy formatted data.table
#' object is returned suitable for plotting using ggplots.
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
#' @param anchor character, one of c("center", "center_unstranded",
#' "left", "left_unstranded")
#' @param return_data.table logical. If TRUE the internal data.table is
#' returned instead of GRanges.  Default is FALSE.
#' @param force_skip_centerFix boolean, if TRUE all query ranges will be
#' used "as is".  This is already the case by default if win_method == "summary"
#' but may have applications where win_method == "sample".
#' @return A GRanges (or data.table if specified) containing fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved
#' from \code{bw_file}
ssvFetchBigwig.single = function(bw_file,
                                 qgr,
                                 win_size = 50,
                                 win_method = c("sample", "summary")[1],
                                 summary_FUN = stats::weighted.mean,
                                 anchor = c("left", "left_unstranded", "center",
                                            "center_unstranded")[3],
                                 return_data.table = FALSE,
                                 force_skip_centerFix = FALSE) {
    stopifnot(is.character(bw_file))
    stopifnot(is(qgr, "GRanges"))
    stopifnot(is.numeric(win_size))
    switch (win_method,
            sample = {
                qgr = prepare_fetch_GRanges_width(qgr, win_size, skip_centerFix = force_skip_centerFix)
                fetch_gr = qgr
                names(fetch_gr) = NULL
                score_gr = rtracklayer::import.bw(bw_file, which = fetch_gr)
                out = viewGRangesWinSample_dt(score_gr, qgr, win_size, anchor = anchor)[]
            },
            summary = {
                fetch_gr = qgr
                names(fetch_gr) = NULL
                score_gr = rtracklayer::import.bw(bw_file, which = fetch_gr)
                out = viewGRangesWinSummary_dt(score_gr, qgr, win_size,
                                               summary_FUN = summary_FUN,
                                               anchor = anchor)[]
            }
    )

    if(!return_data.table){
        out = GRanges(out)
    }
    return(out)
}






