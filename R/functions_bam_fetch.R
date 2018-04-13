#' parse fragLen from MACS2 output
#' @param macs2xls_file character.  an xls file output by MACS2 to parse
#' frag length from
#' @return numeric fragment length
#' @export
#' @examples
#' xls_file = system.file("extdata/test_peaks.xls",
#'     package = "seqsetvis")
#' fragLen_fromMacs2Xls(xls_file)
fragLen_fromMacs2Xls = function(macs2xls_file){
    stopifnot(is.character(macs2xls_file))
    stopifnot(file.exists(macs2xls_file))
    str = utils::read.table(macs2xls_file, nrows = 30, comment.char = "", sep = "\n", stringsAsFactors = FALSE)$V1
    d = str[grepl(" d = ", str)]

    ts = str[grepl(" tag size is determined as ", str)]
    ts = sub("# tag size is determined as ", "", ts)
    ts = sub(" bps", "", ts)
    ts = as.numeric(ts)
    d = as.numeric(strsplit(d, " d = ")[[1]][2])
    return(d)
}

#' calculate fragLen from a bam file for specified regions
#' @param bam_f character or BamFile. bam file to read from.
#' .bai index file must be in same directory
#' @param qgr GRanges.  used as which for ScanBamParam. Can be NULL if it's
#' REALLY important to load the entire bam, force_no_which = TRUE also required.
#' @param ma_distance numeric (integer)  range to use for movingRange.
#' Default is 21.
#' @param n_regions numeric (integer) it's generally overkill to pull all
#' regions at this stage and will slow calculation down.  Default is 100.
#' @param force_no_which logical. if TRUE and qgr is also NULL, the entire
#' bam will be read.
#' @param include_plot_in_output if TRUE ouptut is a list of fragLen and a
#' ggplot showing values considered by calculation. Default is FALSE.
#' @param max_fragLen numeric.  The maximum fragLen to calculate for. Calculation
#' time is directly proportional to this number. Default
#' is 300.
#' @param ... passed to Rsamtools::ScanBamParam, can't be which or what.
#' @return numeric fragment length
#' @import Rsamtools
#' @export
#' @examples
#' bam_file = system.file("extdata/test.bam",
#'     package = "seqsetvis")
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' fragLen_calcStranded(bam_file, qgr)
#' #if plot is included, a list is returned, item 2 is the plot
#' fragLen_calcStranded(bam_file, qgr,
#'   include_plot_in_output = TRUE)[[2]]
fragLen_calcStranded = function(bam_f,
                                qgr,
                                ma_distance = 21,
                                n_regions = 100,
                                force_no_which = FALSE,
                                include_plot_in_output = FALSE,
                                max_fragLen = 300, ...){
    x = y = N = NULL #reserve bindings for data.table
    if(is.null(qgr)){
        if(force_no_which){
            sbParam = Rsamtools::ScanBamParam(what = c("rname", "strand", "pos", "qwidth"), ...)
        }else{
            stop("No qgr was set for ScanBamParam which arg.  ",
                 "This will probably be very slow and uneccessary.  ",
                 "Recall with :\nforce_no_which = TRUE\n if you're certain.")
        }
    }else{
        set.seed(0)
        sbParam = Rsamtools::ScanBamParam(which = sample(qgr, min(length(qgr), n_regions)), what = c("rname", "strand", "pos", "qwidth"), ...)
    }
    bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname, strand = x$strand, start = x$pos, width = x$qwidth)
    })
    bam_dt = data.table::rbindlist(bam_dt, use.names = TRUE, idcol = "which_label")
    bam_dt[, end := start + width - 1L]
    bam_dt[, pos := start]
    bam_dt[strand == "-", pos := end]

    xs = 0:max_fragLen
    perc = vapply(xs, function(x){
        bam_dt[strand == "+", pos := pos + 1L]
        tmp = bam_dt[, .N, by = list(seqnames, pos, strand) ]
        tmp = tmp[, .N, by = list(seqnames, pos)]
        (tmp[N > 1, .N] / nrow(tmp)) * 100
    }, 1)
    ma_perc21 = movingAverage(perc, n = ma_distance)
    fragLenMa = which.max(ma_perc21)
    if(!include_plot_in_output){
        return(fragLenMa)
    }else{
        pdt = data.table(x = xs, raw = perc, moving_average = ma_perc21)

        pdt = data.table::melt(pdt, id.vars = "x", variable.name = "transform", value.name = "y")
        p = ggplot(pdt) +
            labs(x = "Fragment Length", y = "% strand match",
                 title = paste(basename(bam_f)),
                 subtitle = paste("Fragment Length determined by strand match maximization",
                                  paste("Moving averge window", ma_distance, "applied"),
                                  sep = "\n")) +
            geom_line(aes(x = x, y = y, color = transform)) +
            scale_color_manual(values = c(raw = "black", moving_average = "red")) +
            annotate("line", x = rep(fragLenMa, 2), y = range(perc), color = "green") +
            annotate("label", x = fragLenMa, y = mean(range(perc)), label = fragLenMa, color = "black")
        return(list(fragLenMa, p))
    }
}

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

#' fetch a bam file pileup with the ability to consider cross strand correlation
#' @param bam_f character or BamFile to load
#' @param qgr GRanges regions to fetchs
#' @param fragLen numeric, NULL, or NA.  if numeric, supplied value is used.
#' if NULL, value is calculated with fragLen_calcStranded
#' if NA, raw bam pileup with no cross strand shift is returned.
#' @param target_strand character. if one of "+" or "-", reads are filtered
#' to match. ignored if any other value.
#' @param ... passed to ScanBamParam(), can't be which or what.
#' @return GRanges containing tag pileup values in score meta column.  tags are
#' optionally extended to fragment length (fragLen) prior to pile up.
#' @export
#' @examples
#' bam_file = system.file("extdata/test.bam", package = "seqsetvis")
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' fetchBam(bam_file, qgr)
#' fetchBam(bam_file, qgr, fragLen = 180, target_strand = "+")
fetchBam = function(bam_f, qgr, fragLen = NULL, target_strand = c("*", "+", "-")[1], ...){
    if(is.null(fragLen)){
        fragLen = fragLen_calcStranded(bam_f, qgr)
        message("fragLen was calculated as: ", fragLen)
    }
    if(!is.na(fragLen)){
        stopifnot(is.numeric(fragLen))
        stopifnot(fragLen %% 1 == 0)
        stopifnot(fragLen >= 1)
    }
    sbParam = Rsamtools::ScanBamParam(which = qgr, what = c("rname", "strand", "pos", "qwidth"), ...)
    bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname, strand = x$strand, start = x$pos, width = x$qwidth)
    })
    bam_dt = data.table::rbindlist(bam_dt, use.names = TRUE, idcol = "which_label")
    bam_dt[, end := start + width - 1L]

    ext_dt = copy(bam_dt)

    if(target_strand == "+"){
        ext_dt = ext_dt[strand == "+"]
    }
    if(target_strand == "-"){
        ext_dt = ext_dt[strand == "-"]
    }

    if(!is.na(fragLen)){#extension to fragLen
        ext_dt[strand == "+", end := start + as.integer(fragLen) - 1L]
        ext_dt[strand == "-", start := end - as.integer(fragLen) + 1L]
    }

    ext_cov = coverage(split(GRanges(ext_dt), ext_dt$which_label))
    score_gr = GRanges(ext_cov)
    if(target_strand == "+"){
        strand(score_gr) = "+"
    }
    if(target_strand == "-"){
        strand(score_gr) = "-"
    }
    return(score_gr)
}

#' fetch a windowed version of a bam file, returns data.table
#'
#' @param bam_f character or BamFile to load
#' @param qgr GRanges regions to fetch
#' @param fragLen numeric, NULL, or NA.  if numeric, supplied value is used.
#' if NULL, value is calculated with fragLen_calcStranded
#' if NA, raw bam pileup with no cross strand shift is returned.
#' @param win_size numeric >=1.  pileup grabbed every win_size bp
#' @param target_strand character. if one of "+" or "-", reads are filtered
#' @return tidy data.table with GRanges compatible columns.  pileup is
#' calculated only every win_size bp.
#' @export
#' @examples
#' bam_file = system.file("extdata/test.bam",
#'     package = "seqsetvis")
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' bam_dt = fetchWindowedBam_dt(bam_file, qgr)
#' bam_dt = fetchWindowedBam_dt(bam_file, qgr, fragLen = 180,
#'     win_size = 10, target_strand = "+")
fetchWindowedBam_dt = function(bam_f, qgr, fragLen = NULL, win_size = 50, target_strand = c("*", "+", "-")[1]) {
    score_gr = fetchBam(bam_f, qgr, fragLen, target_strand)
    viewGRangesWindowed_dt(score_gr, qgr, win_size)
}

#' fetch a windowed version of a bam file, returns GRanges
#'
#' @param bam_f character or BamFile to load
#' @param qgr GRanges regions to fetchs
#' @param fragLen numeric, NULL, or NA.  if numeric, supplied value is used.
#' if NULL, value is calculated with fragLen_calcStranded
#' if NA, raw bam pileup with no cross strand shift is returned.
#' @param win_size numeric >=1.  pileup grabbed every win_size bp
#' @param target_strand character. if one of "+" or "-", reads are filtered
#' accordingly. ignored if any other value.
#' @return tidy GRanges with pileups from bam file.  pileup is
#' calculated only every win_size bp.
#' @export
#' @examples
#' bam_file = system.file("extdata/test.bam",
#'     package = "seqsetvis")
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' bam_gr = fetchWindowedBam(bam_file, qgr)
#' bam_gr = fetchWindowedBam(bam_file, qgr, fragLen = 180,
#'     win_size = 10, target_strand = "+")
fetchWindowedBam = function(bam_f, qgr, fragLen = NULL, win_size = 50, target_strand = c("*", "+", "-")[1]) {
    GRanges(fetchWindowedBam_dt(bam_f, qgr, fragLen, win_size, target_strand))
}
