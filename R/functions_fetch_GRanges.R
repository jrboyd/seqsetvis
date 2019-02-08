#' Fetch coverage values for a list of GRanges.
#'
#' \code{ssvFetchGRanges} Gets coverage values for each region of the query
#' GRanges (\code{qgr}). Values correspond to the center of each window of
#' size \code{win_size}.  A tidy formatted data.table
#' object is returned suitable for plotting using ggplots.
#'
#' @param grs a list of GRanges for which to calculate coverage.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#'   interval should be identical and evenly divisible by \code{win_size}.
#' @param unique_names The column name where unique_names are stored.
#' Default is 'sample'
#' @param names_variable The column name where unique_names are stored.
#' Default is 'sample'
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
#' @param names_variable
#' @param return_data.table logical. If TRUE the internal data.table is
#' returned instead of GRanges.  Default is FALSE.
#' @param n_cores integer number of cores to use.
#' Uses mc.cores option if not supplied.
#' @return A tidy formatted GRanges (or data.table if specified) containing
#' fetched values.
#' @export
#' @examples
#' ssvFetchGRanges(CTCF_in_10a_narrowPeak_grs, CTCF_in_10a_overlaps_gr, win_size = 200)
ssvFetchGRanges = function(grs, qgr,
                           unique_names = names(file_paths),
                           names_variable = "sample",
                           win_size = 50,
                           win_method = c("sample", "summary")[1],
                           summary_FUN = function(x, w)max(x),
                           target_strand = c("*", "+", "-", "both")[1],
                           anchor = c("left", "left_unstranded", "center",
                                      "center_unstranded")[3],
                           return_data.table = FALSE,
                           n_cores = getOption("mc.cores", 1)){

    # gr_dt = lapply(grs, function(x){
    gr_dt = parallel::mclapply(grs, mc.cores = n_cores, function(x){
        switch(win_method,
               sample = {
                   qgr = prepare_fetch_GRanges(qgr, win_size)
                   if(target_strand == "both"){
                       pos_gr = GRanges(coverage(subset(x, strand == "+")))
                       neg_gr = GRanges(coverage(subset(x, strand == "-")))
                       pos_dt = viewGRangesWinSample_dt(pos_gr, qgr,
                                                        win_size, anchor = anchor)
                       neg_dt = viewGRangesWinSample_dt(neg_gr, qgr,
                                                        win_size, anchor = anchor)
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

                       out = viewGRangesWinSample_dt(score_gr, qgr,
                                                     win_size, anchor = anchor)
                       out[, strand := target_strand]
                   }
                   out
               },
               summary = {
                   if(target_strand == "both"){
                       pos_gr = GRanges(coverage(subset(x, strand == "+")))
                       neg_gr = GRanges(coverage(subset(x, strand == "-")))
                       pos_dt = viewGRangesWinSummary_dt(pos_gr, qgr, win_size,
                                                         summary_FUN = summary_FUN,
                                                         anchor = anchor)

                       neg_dt = viewGRangesWinSummary_dt(neg_gr, qgr, win_size,
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
                       out = viewGRangesWinSummary_dt(score_gr, qgr, win_size,
                                                      summary_FUN = summary_FUN,
                                                      anchor = anchor)
                       out[, strand := target_strand]
                   }
                   out
               })
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

