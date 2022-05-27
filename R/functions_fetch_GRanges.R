#' Fetch coverage values for a list of GRanges.
#'
#' \code{ssvFetchGRanges} Gets coverage values for each region of the query
#' GRanges (\code{qgr}). Values correspond to the center of each window of size
#' \code{win_size}.  A tidy formatted data.table object is returned suitable for
#' plotting using ggplots.
#'
#' @param grs a list of GRanges for which to calculate coverage.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#'   interval should be identical and evenly divisible by \code{win_size}.
#' @param file_attribs data.frame (1 row per item in grs) containing attributes
#'   to append to results.
#' @param unique_names The column name where unique_names are stored. Default is
#'   'sample'
#' @param names_variable The column name where unique_names are stored. Default
#'   is 'sample'
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @param win_method character.  one of c("sample", "summary").  Determines if
#'   \code{\link{viewGRangesWinSample_dt}} or
#'   \code{\link{viewGRangesWinSummary_dt}} is used to represent each region in
#'   qgr.
#' @param summary_FUN function.  only relevant if win_method is "summary".
#'   passed to \code{\link{viewGRangesWinSummary_dt}}.
#' @param target_strand character. if one of "+" or "-", reads are filtered to
#'   match. ignored if any other value.
#' @param use_coverage boolean or NULL, if TRUE, query regions are scored by the
#'   number of intervals overlapping.  Default of NULL checks if attrib_var is
#'   "score" and uses coverage if so.
#' @param attrib_var character, column in mcols of GRanges to pull values from.
#'   Default of "score" is compatible with internal coverage calculation or
#'   bedgraph-like files.
#' @param fill_value numeric or character value to use where queried regions are
#'   empty.  Default is 0 and appropriate for both calculated coverage and
#'   bedgraph/bigwig like files.  Will automatically switch to "MISSING" if data
#'   is guessed to be qualitative.
#' @param anchor character, one of c("center", "center_unstranded", "left",
#'   "left_unstranded")
#' @param return_data.table logical. If TRUE the internal data.table is returned
#'   instead of GRanges.  Default is FALSE.
#' @param n_cores integer number of cores to use. Uses mc.cores option if not
#'   supplied.
#' @param force_skip_centerFix boolean, if TRUE all query ranges will be
#' used "as is".  This is already the case by default if win_method == "summary"
#' but may have applications where win_method == "sample".
#' @return A tidy formatted GRanges (or data.table if specified) containing
#'   fetched values.
#' @export
#' @examples
#' ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, CTCF_in_10a_overlaps_gr, win_size = 200)
ssvFetchGRanges = function(grs,
                           qgr,
                           file_attribs = data.frame(matrix(0, nrow = length(grs), ncol = 0)),
                           unique_names = names(grs),
                           names_variable = "sample",
                           win_size = 50,
                           win_method = c("sample", "summary")[1],
                           summary_FUN = function(x, w)max(x),
                           target_strand = c("*", "+", "-", "both")[1],
                           use_coverage = NULL,
                           attrib_var = "score",
                           fill_value = 0,
                           anchor = c("left", "left_unstranded", "center",
                                      "center_unstranded")[3],
                           return_data.table = FALSE,
                           n_cores = getOption("mc.cores", 1),
                           force_skip_centerFix = FALSE){
    if(!is.list(grs)){
        if(is(grs, "GRangesList")){
            grs = as.list(grs)
        }else{
            grs = list(grs)
        }
    }
    all_gr = all(vapply(grs, function(x){
        "GRanges" %in% class(x)
    }, FUN.VALUE = TRUE))
    stopifnot(all_gr)
    if(is.null(unique_names)){
        unique_names = paste("regions", LETTERS[seq_along(grs)])
    }
    if(is.null(file_attribs[[names_variable]])){
        file_attribs[[names_variable]] = unique_names
    }
    if(is.null(use_coverage)){
        use_coverage = attrib_var == "score"
    }
    # gr_dt = lapply(grs, function(x){
    gr_dt = parallel::mclapply(grs, mc.cores = n_cores, function(x){
        if(use_coverage){
            fetchGRanges_as_coverage(x, qgr, win_size, win_method,
                                     summary_FUN, target_strand,
                                     anchor, fill_value = fill_value,
                                     force_skip_centerFix = force_skip_centerFix)
        }else{
            fetchGRanges_by_attrib_var(x, qgr, win_size, win_method,
                                       summary_FUN, target_strand,
                                       anchor, fill_value = fill_value,
                                       attrib_var = attrib_var,
                                       force_skip_centerFix = force_skip_centerFix)
        }
    })
    for(i in seq_along(gr_dt)){
        for(attrib in colnames(file_attribs)){
            gr_dt[[i]][[attrib]] = file_attribs[[attrib]][i]
        }

    }
    gr_dt = rbindlist(gr_dt)#, use.names = TRUE, idcol = names_variable)
    if(!return_data.table){
        gr_dt = GRanges(gr_dt)
    }
    gr_dt
}

# ssvFetchGRanges = function(grs,
#                            qgr,
#                            file_attribs = data.frame(matrix(0, nrow = length(grs), ncol = 0)),
#                            unique_names = names(grs),
#                            names_variable = "sample",
#                            win_size = 50,
#                            win_method = c("sample", "summary")[1],
#                            summary_FUN = function(x, w)max(x),
#                            target_strand = c("*", "+", "-", "both")[1],
#                            use_coverage = NULL,
#                            attrib_var = "score",
#                            anchor = c("left", "left_unstranded", "center",
#                                       "center_unstranded")[3],
#                            return_data.table = FALSE,
#                            n_cores = getOption("mc.cores", 1)){
fetchGRanges_as_coverage = function(x, qgr, win_size, win_method, summary_FUN,
                                    target_strand, anchor,
                                    fill_value = 0, force_skip_centerFix = FALSE){
    switch(win_method,
           sample = {
               qgr = prepare_fetch_GRanges_width(qgr,
                                           win_size,
                                           skip_centerFix = force_skip_centerFix)
               if(target_strand == "both"){
                   pos_gr = GRanges(coverage(subset(x, strand == "+")))
                   neg_gr = GRanges(coverage(subset(x, strand == "-")))
                   pos_dt = viewGRangesWinSample_dt(
                       pos_gr, qgr,
                       win_size,
                       attrib_var = "score",
                       anchor = anchor)
                   neg_dt = viewGRangesWinSample_dt(
                       neg_gr, qgr,
                       win_size,
                       attrib_var = "score",
                       anchor = anchor)
                   pos_dt[, strand := "+"]
                   neg_dt[, strand := "-"]
                   out = rbind(
                       pos_dt,
                       neg_dt
                   )
               }else{
                   if(target_strand == "*"){
                       score_gr = GRanges(coverage(x))
                   }else{
                       score_gr = GRanges(coverage(subset(x, strand == target_strand)))
                   }

                   out = viewGRangesWinSample_dt(
                       score_gr,
                       qgr,
                       win_size,
                       attrib_var = "score",
                       anchor = anchor)
                   out[, strand := target_strand]
               }
               out
           },
           summary = {
               if(target_strand == "both"){
                   pos_gr = GRanges(coverage(subset(x, strand == "+")))
                   neg_gr = GRanges(coverage(subset(x, strand == "-")))
                   pos_dt = viewGRangesWinSummary_dt(
                       pos_gr,
                       qgr,
                       win_size,
                       summary_FUN = summary_FUN,
                       attrib_var = "score",
                       fill_value = fill_value,
                       anchor = anchor)

                   neg_dt = viewGRangesWinSummary_dt(
                       neg_gr,
                       qgr,
                       win_size,
                       attrib_var = "score",
                       fill_value = fill_value,
                       summary_FUN = summary_FUN,
                       anchor = anchor)
                   pos_dt[, strand := "+"]
                   neg_dt[, strand := "-"]
                   out = rbind(
                       pos_dt,
                       neg_dt
                   )
               }else{
                   if(target_strand == "*"){
                       score_gr = GRanges(coverage(x))
                   }else{
                       score_gr = GRanges(coverage(subset(x, strand == target_strand)))
                   }
                   out = viewGRangesWinSummary_dt(
                       score_gr,
                       qgr,
                       win_size,
                       attrib_var = "score",
                       fill_value = fill_value,
                       summary_FUN = summary_FUN,
                       anchor = anchor)
                   out[, strand := target_strand]
               }
               out
           })
}

fetchGRanges_by_attrib_var = function(x, qgr, win_size, win_method, summary_FUN,
                                      target_strand, anchor,
                                      fill_value = "MISSING", attrib_var = "name",
                                      force_skip_centerFix = FALSE){
    switch(win_method,
           sample = {
               qgr = prepare_fetch_GRanges_width(qgr, win_size,
                                           skip_centerFix = force_skip_centerFix)
               if(target_strand == "both"){
                   pos_gr = subset(x, strand == "+")
                   neg_gr = subset(x, strand == "-")
                   pos_dt = viewGRangesWinSample_dt(
                       pos_gr,
                       qgr,
                       win_size,
                       attrib_var = attrib_var,
                       fill_value = fill_value,
                       anchor = anchor
                   )
                   neg_dt = viewGRangesWinSample_dt(
                       neg_gr,
                       qgr,
                       win_size,
                       attrib_var = attrib_var,
                       fill_value = fill_value,
                       anchor = anchor
                   )
                   pos_dt[, strand := "+"]
                   neg_dt[, strand := "-"]
                   out = rbind(
                       pos_dt,
                       neg_dt
                   )
               }else{
                   if(target_strand != "*"){
                       x = subset(x, strand == target_strand)
                   }

                   out = viewGRangesWinSample_dt(
                       x,
                       qgr,
                       win_size,
                       attrib_var = attrib_var,
                       fill_value = fill_value,
                       anchor = anchor)
                   out[, strand := target_strand]
               }
               out
           },
           summary = {
               if(target_strand == "both"){
                   pos_gr = subset(x, strand == "+")
                   neg_gr = subset(x, strand == "-")
                   pos_dt = viewGRangesWinSummary_dt(
                       pos_gr,
                       qgr,
                       win_size,
                       attrib_var = attrib_var,
                       fill_value = fill_value,
                       summary_FUN = summary_FUN,
                       anchor = anchor
                   )

                   neg_dt = viewGRangesWinSummary_dt(
                       neg_gr,
                       qgr,
                       win_size,
                       attrib_var = attrib_var,
                       fill_value = fill_value,
                       summary_FUN = summary_FUN,
                       anchor = anchor
                   )
                   pos_dt[, strand := "+"]
                   neg_dt[, strand := "-"]
                   out = rbind(
                       pos_dt,
                       neg_dt
                   )
               }else{
                   if(target_strand != "*"){
                       x = subset(x, strand == target_strand)
                   }
                   out = viewGRangesWinSummary_dt(
                       x,
                       qgr,
                       win_size,
                       attrib_var = attrib_var,
                       fill_value = fill_value,
                       summary_FUN = summary_FUN,
                       anchor = anchor
                   )
                   out[, strand := target_strand]
               }
               out
           })
}
