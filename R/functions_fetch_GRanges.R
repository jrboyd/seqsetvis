#' Title
#'
#' @param grs
#' @param qgr
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
                           summary_FUN = stats::weighted.mean,
                           target_strand = c("*", "+", "-", "both")[1],
                           anchor = c("left", "left_unstranded", "center",
                                      "center_unstranded")[3],
                           names_variable = "sample",
                           return_data.table = FALSE){

    gr_dt = lapply(grs, function(x){
        switch(win_method,
               sample = {
                   viewGRangesWinSample_dt(GRanges(coverage(x)),
                                           qgr = qgr,
                                           anchor = anchor,
                                           window_size = win_size)
               },
               summary = {
                   viewGRangesWinSummary_dt(GRanges(coverage(x)),
                                            qgr = qgr,
                                            anchor = anchor,
                                            n_tiles = win_size)
               })

    })
    rbindlist(gr_dt, use.names = TRUE, idcol = "sample")
}



ssvFetchGRanges.single = function(grs, qgr,
                                  unique_names = names(grs),
                                  win_size = 50,
                                  win_method = c("sample", "summary")[1],
                                  summary_FUN = stats::weighted.mean,
                                  target_strand = c("*", "+", "-", "both")[1],
                                  anchor = c("left", "left_unstranded", "center",
                                             "center_unstranded")[3],
                                  names_variable = "sample",
                                  return_data.table = FALSE){
    stopifnot(is.character(win_method))
    stopifnot(length(win_method) == 1)
    stopifnot(class(qgr) == "GRanges")
    stopifnot(win_method %in% c("sample", "summary"))
    stopifnot(is.function(summary_FUN))
    stopifnot(target_strand %in% c("*", "+", "-", "both"))
    stopifnot(anchor %in% c("left", "left_unstranded", "center",
                            "center_unstranded"))
    stopifnot(splice_strategy %in% c("none", "ignore", "add", "only", "splice_count"))
    if(splice_strategy == "splice_count"){
        return(fetchBam(bam_f, qgr, NA, "*",
                        max_dupes, splice_strategy))
    }
    switch (
        win_method,
        sample = {
            qgr = prepare_fetch_GRanges(qgr, win_size)
            if(target_strand == "both"){
                pos_gr = fetchBam(bam_f, qgr, fragLen, "+",
                                  max_dupes, splice_strategy)
                neg_gr = fetchBam(bam_f, qgr, fragLen, "-",
                                  max_dupes, splice_strategy)
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
                score_gr = fetchBam(bam_f, qgr, fragLen, target_strand,
                                    max_dupes, splice_strategy)
                out = viewGRangesWinSample_dt(score_gr, qgr,
                                              win_size, anchor = anchor)
                out[, strand := target_strand]
            }


        },
        summary = {
            if(target_strand == "both"){
                pos_gr = fetchBam(bam_f, qgr, fragLen, "+",
                                  max_dupes, splice_strategy)
                neg_gr = fetchBam(bam_f, qgr, fragLen, "-",
                                  max_dupes, splice_strategy)
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
                score_gr = fetchBam(bam_f, qgr, fragLen, target_strand,
                                    max_dupes, splice_strategy)
                out = viewGRangesWinSummary_dt(score_gr, qgr, win_size,
                                               summary_FUN = summary_FUN,
                                               anchor = anchor)
                out[, strand := target_strand]
            }
        }
    )
    if(!return_data.table){
        out = GRanges(out)
    }
    return(out)
}
