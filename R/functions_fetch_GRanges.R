#' Title
#'
#' @param grs a list of GRanges for which to calculate coverage.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#'   interval should be identical and evenly divisible by \code{win_size}.
#' @param unique_names
#' @param win_size
#' @param win_method
#' @param summary_FUN
#' @param target_strand
#' @param anchor
#' @param names_variable
#' @param return_data.table
#'
#' @return
#' @export
#'
#' @examples
ssvFetchGRanges = function(grs, qgr,
                           unique_names = names(file_paths),
                           win_size = 50,
                           win_method = c("sample", "summary")[1],
                           summary_FUN = function(x, w)max(x),
                           target_strand = c("*", "+", "-", "both")[1],
                           anchor = c("left", "left_unstranded", "center",
                                      "center_unstranded")[3],
                           names_variable = "sample",
                           return_data.table = FALSE){

    gr_dt = lapply(grs, function(x){
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
                       score_gr = GRanges(coverage(subset(x, strand == target_strand)))
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
                       score_gr = GRanges(coverage(subset(x, strand == target_strand)))
                       out = viewGRangesWinSummary_dt(score_gr, qgr, win_size,
                                                      summary_FUN = summary_FUN,
                                                      anchor = anchor)
                       out[, strand := target_strand]
                   }
                   out
               })
    })
    gr_dt= rbindlist(gr_dt, use.names = TRUE, idcol = names_variable)
    if(!return_data.table){
        gr_dt = GRanges(gr_dt)
    }
    gr_dt
}

