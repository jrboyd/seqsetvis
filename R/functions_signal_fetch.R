# functions useful for fetching signal data regardless of source:
# bam, bigwig, etc.

#' get a windowed sampling of score_gr
#'
#' Summarizes score_gr by grabbing value of "score" every window_size bp.
#' Columns in output data.table are:
#' standard GRanges columns: seqnames, start, end, width, strand
#' id - matched to names(score_gr). if names(score_gr) is missing,
#' added as 1:length(score_gr)
#' y - value of score from score_gr
#' x - relative bp position
#'
#' @param score_gr GRanges with a "score" metadata columns.
#' @param qgr regions to view by window.
#' @param window_size qgr will be represented by value from score_gr every
#' window_size bp.
#' @param x0 character. controls how x value is derived from position for
#' each region in qgr.  0 may be the left side or center.  If not unstranded,
#' x coordinates are flipped for (-) strand.
#' @return data.table that is GRanges compatible
#' @export
#' @examples
#' bam_file = system.file("extdata/test.bam",
#'     package = "seqsetvis")
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' qgr = GenomicRanges::resize(qgr, width = 500, fix = "center")
#' bam_gr = fetchBam(bam_file, qgr)
#' bam_dt = viewGRangesWindowed_dt(bam_gr, qgr, 50)
#'
#' bw_file = system.file("extdata/MCF10A_CTCF_FE_random100.bw",
#'     package = "seqsetvis")
#' bw_gr = rtracklayer::import.bw(bw_file, which = qgr)
#' bw_dt = viewGRangesWindowed_dt(bw_gr, qgr, 50)
viewGRangesWindowed_dt = function(score_gr, qgr, window_size,
                                  x0 = c("center", "center_unstranded",
                                         "left", "left_unstranded")[1]){
    x = id = NULL
    stopifnot(class(score_gr) == "GRanges")
    stopifnot(!is.null(score_gr$score))
    stopifnot(class(qgr) == "GRanges")
    stopifnot(is.numeric(window_size))
    stopifnot(window_size >= 1)
    stopifnot(window_size %% 1 == 0)
    stopifnot(x0 %in% c("center", "center_unstranded", "left", "left_unstranded"))
    windows = slidingWindows(qgr, width = window_size, step = window_size)
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
    olaps = suppressWarnings(data.table::as.data.table(findOverlaps(query = windows, subject = score_gr)))
    # patch up missing/out of bound data with 0
    missing_idx = setdiff(seq_along(windows), olaps$queryHits)
    if (length(missing_idx) > 0) {
        olaps = rbind(olaps, data.table::data.table(queryHits = missing_idx, subjectHits = length(score_gr) + 1))[order(queryHits)]
        score_gr = c(score_gr, GRanges(seqnames(score_gr)[length(score_gr)], IRanges::IRanges(1, 1), score = 0))
    }
    # set y and output windows = windows[olaps$queryHits]
    windows$y = score_gr[olaps$subjectHits]$score
    score_dt = data.table::as.data.table(windows)

    shift = round(window_size/2)
    switch(x0,
           center = {
               score_dt[, `:=`(x, start - min(start) + shift), by = id]
               score_dt[, `:=`(x, x - round(mean(x))), by = id]
               score_dt[strand == "-", x := -1*x]
           },
           center_unstranded = {
               score_dt[, `:=`(x, start - min(start) + shift), by = id]
               score_dt[, `:=`(x, x - round(mean(x))), by = id]
           },
           left = {
               score_dt[, x := -1]
               score_dt[strand != "-", `:=`(x, start - min(start) + shift), by = id]
               #flip negative
               score_dt[strand == "-", `:=`(x, -1*(end - max(end) - shift)), by = id]
           },
           left_unstranded = {
               score_dt[, `:=`(x, start - min(start) + shift), by = id]
           }
    )

    score_dt[, `:=`(start, start - shift + 1)]
    score_dt[, `:=`(end, end + window_size - shift)]
    if(x0 == "center"){

    }
    score_dt
}

#' prepares GRanges for windowed fetching.
#'
#' output GRanges parallels input with consistent width evenly divisible by
#' win_size.  Has warning if GRanges needed resizing, otherwise no warning
#' and input GRanges is returned unchanged.
#'
#' @param qgr GRanges to prepare
#' @param win_size numeric window size for fetch
#' @param target_size numeric final width of qgr
prepare_fetch_GRanges = function(qgr, win_size, target_size = NULL){
    if(length(unique(width(qgr))) > 1 || width(qgr)[1] %% win_size != 0 ){
        if(is.null(target_size)){
            target_size = quantile(width(qgr), .75)
            target_size = round(target_size / win_size) * win_size
        }
        if(target_size %% win_size != 0){
            stop("target_size: ", target_size,
                 " not evenly divisible by win_size: ", win_size)
        }

        qgr = centerFixedSizeGRanges(qgr, fixed_size = target_size)
        warning("fetchWindowedBigwigList_dt requires widths of qgr be ",
                "identical and evenly divisible by win_size.",
                "\nA fixed width of ",
                target_size, " was applied based on the data provided.")
    }
    return(qgr)
}

#' Quantile width determination strategy
#'
#' Returns the lowest multiple of win_size greater than
#' min_quantile quantile of width(qgr)
quantileGRangesWidth = function(qgr,
                                min_quantile = .75,
                                win_size = 1){

    stopifnot(class(qgr) == "GRanges")
    stopifnot(is.numeric(min_quantile), is.numeric(win_size))
    stopifnot(min_quantile >= 0 && min_quantile <= 1)
    stopifnot(length(min_quantile) == 1 && length(win_size) == 1)
    stopifnot(win_size%%1==0)
    stopifnot(win_size >= 1)
    qwidth = quantile(width(qgr), min_quantile)
    fwidth = ceiling(qwidth / win_size) * win_size
    return(fwidth)
}

#' Derive a new GRanges of consistent width based on quantile.
#' Width is selected by rounding up to the lowest multiple of win_size greater than
#' min_quantile quantile of widths.
#'
#' @param qgr GRanges. To be resized.
#' @param min_quantile numeric [0,1]. The quantile level final width must be
#' greater than. default is 0.75
#' @param win_size integer > 0.  final width must be a multiple of win_size.
#' @return a GRanges derived from qgr (length and order match).  All ranges
#' are of same width and centered on old.  Width is at least minimum quantile
#' and a multiple of win_size
fixGRangesWidth = function(qgr,
                           min_quantile = .75,
                           win_size = 1,
                           anchor = c("center", "start")[1]){
    stopifnot(class(qgr) == "GRanges")


    fwidth = quantileGRangesWidth(qgr, min_quantile, win_size)
    qgr = setGRangesWidth(qgr = qgr, fwidth = fwidth, anchor = anchor)

}

#' Return GRanges with single width
#'
#' Essentially works like GenomicRanges::resize() but repeated applications
#' of center do not cause rounding induced drift.
setGRangesWidth = function(qgr, fwidth, anchor = c("center", "start")[1]){
    stopifnot(class(qgr) == "GRanges")
    stopifnot(is.numeric(fwidth))
    stopifnot(anchor %in% c("center", "start"))
    switch(anchor,
           center = {
               centerFixedSizeGRanges(qgr, fixed_size = fwidth)
           },
           start = {
               resize(qgr, width = fwidth, fix = "start")
           }
    )
}

