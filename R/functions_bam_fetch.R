#' parse fragLen from MACS2 output
#' @export
fragLen_fromMacs2Xls = function(macs2xls_file){
    stopifnot(is.character(macs2xls_file))
    stopifnot(file.exists(macs2xls_file))
    str = read.table(macs2xls_file, nrows = 30, comment.char = "", sep = "\n", stringsAsFactors = F)$V1
    d = str[grepl(" d = ", str)]

    ts = str[grepl(" tag size is determined as ", str)]
    ts = sub("# tag size is determined as ", "", ts)
    ts = sub(" bps", "", ts)
    ts = as.numeric(ts)
    d = as.numeric(strsplit(d, " d = ")[[1]][2])
    return(d)
}

#' calculate fragLen from a bam file for specified regions
#' @import Rsamtools
#' @param ... passed to Rsamtools::ScanBamParam, can't be which or what.
#' @export
fragLen_calcStranded = function(bam_f,
                                qgr,
                                n_regions = 100,
                                force_no_which = FALSE,
                                return_fragLen_only = TRUE,
                                max_fragLen = 300, ...){
    if(is.null(qgr)){
        if(force_no_which){
            sbParam = Rsamtools::ScanBamParam(what = c("rname", "strand", "pos", "qwidth"), ...)
        }else{
            stop("No qgr was set for ScanBamParam which arg.  ",
                 "This will probably be very slow and uneccessary.  ",
                 "Recall with :\nforce_no_which = TRUE\n if you're certain.")
        }
    }else{
        if(is.character(qgr)){
            peak_gr = easyLoad_narrowPeak(qgr)[[1]]
            qgr = peak_gr[order(peak_gr$signalValue, decreasing = T)][1:n_regions]
            qgr = GenomicRanges::resize(qgr, width = 2500, fix = "center")
        }
        sbParam = Rsamtools::ScanBamParam(which = sample(qgr, min(length(qgr), n_regions)), what = c("rname", "strand", "pos", "qwidth"), ...)
    }
    bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname, strand = x$strand, start = x$pos, width = x$qwidth)
    })
    bam_dt = data.table::rbindlist(bam_dt, use.names = T, idcol = "which_label")
    bam_dt[, end := start + width - 1L]
    bam_dt[, pos := start]
    bam_dt[strand == "-", pos := end]

    xs = 0:max_fragLen
    perc = pbapply::pbsapply(xs, function(x){
        bam_dt[strand == "+", pos := pos + 1L]
        tmp = bam_dt[, .N, by = .(seqnames, pos, strand) ]
        tmp = tmp[, .N, by = .(seqnames, pos)]
        (tmp[N > 1, .N] / nrow(tmp)) * 100
    })
    ma_perc21 = movingAverage(perc, n = 21)
    fragLenMa = which.max(ma_perc21)
    if(return_fragLen_only){
        return(fragLenMa)
    }else{
        pdt = data.table(x = xs, raw = perc, MA_21 = ma_perc21)
        pdt = data.table::melt(pdt, id.vars = "x", variable.name = "transform", value.name = "y")
        p = ggplot(pdt) +
            labs(x = "Fragment Length", y = "% strand match",
                 title = paste(basename(bam_f)),
                 subtitle = paste("Fragment Length determined by strand match maximization",
                                  "Moving averge window 21 applied",
                                  sep = "\n")) +
            geom_line(aes(x = x, y = y, color = transform)) +
            scale_color_manual(values = c(raw = "black", MA_21 = "red")) +
            annotate("line", x = rep(fragLenMa, 2), y = range(perc), color = "green") +
            annotate("label", x = fragLenMa, y = mean(range(perc)), label = fragLenMa, color = "black")
        return(list(fragLenMa, p))
    }
}

#' get a windowed sampling of score_gr
#' @export
windowViewGRanges_dt = function(score_gr, query_gr, window_size){
    x = id = NULL
    windows = slidingWindows(query_gr, width = window_size, step = window_size)
    if (is.null(query_gr$id)) {
        if (!is.null(names(query_gr))) {
            query_gr$id = names(query_gr)
        } else {
            query_gr$id = paste0("region_", seq_along(query_gr))
        }
    }
    names(windows) = query_gr$id
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
    score_dt[, `:=`(x, start - min(start) + window_size/2), by = id]
    score_dt[, `:=`(x, x - round(mean(x))), by = id]
    shift = round(window_size/2)
    score_dt[, `:=`(start, start - shift + 1)]
    score_dt[, `:=`(end, end + window_size - shift)]
    score_dt
}

#' fetch a bam file pileup with the ability to consider cross strand correlation
#' @param target_strand character. if one of "+" or "-", reads are filtered
#' @param fragLen numeric, NULL, or NA.  if numeric, supplied value is used.
#' @param ... passed to ScanBamParam(), can't be which or what.
#' if NULL, value is calculated with fragLen_calcStranded
#' if NA, raw bam pileup with no cross strand shift is returned.
#' accordingly. ignored if any other value.
#' @export
fetchBam = function(bam_f, qgr, fragLen = NULL, win_size = 50, target_strand = c("*", "+", "-")[1], ...){
    if(is.null(fragLen)){
        fragLen = fragLen_calcStranded(bam_f, qgr)
        message("fragLen was calculated as: ", fragLen)
    }
    if(!is.na(fragLen)){
        stopifnot(is.numeric(fragLen))
        stopifnot(fragLen %% 1 == 0)
        stopifnot(fragLen >= 1)
    }

    if(is.character(qgr)){
        qgr = easyLoad_narrowPeak(qgr)[[1]]
    }
    if(length(unique(width(qgr))) > 1 || width(qgr)[1] %% win_size != 0){
        qgr = fixGRangesWidth(qgr, min_quantile = .75, win_size = win_size)
        target_size = width(qgr)[1]
        message("A fixed width of ",
                target_size, " was applied based on the data provided.")
    }
    sbParam = Rsamtools::ScanBamParam(which = qgr, what = c("rname", "strand", "pos", "qwidth"), ...)
    bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname, strand = x$strand, start = x$pos, width = x$qwidth)
    })
    bam_dt = data.table::rbindlist(bam_dt, use.names = T, idcol = "which_label")
    bam_dt[, end := start + width - 1L]

    ext_dt = copy(bam_dt)

    if(target_strand == "+"){
        ext_dt = ext_dt[strand == "+"]
    }
    if(target_strand == "-"){
        ext_dt = ext_dt[strand == "-"]
    }

    if(!is.na(fragLen)){#extension to fragLen
        ext_dt[strand == "+", end := start + fragLen - 1L]
        ext_dt[strand == "-", start := end - fragLen + 1L]
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
#' @param target_strand character. if one of "+" or "-", reads are filtered
#' @param fragLen numeric, NULL, or NA.  if numeric, supplied value is used.
#' if NULL, value is calculated with fragLen_calcStranded
#' if NA, raw bam pileup with no cross strand shift is returned.
#' accordingly. ignored if any other value.
#' @export
fetchWindowedBam_dt = function(bam_f, qgr, fragLen = NULL, win_size = 50, target_strand = c("*", "+", "-")[1]) {
    score_gr = fetchBam(bam_f, qgr, fragLen, win_size, target_strand)
    windowViewGRanges_dt(score_gr, qgr, win_size)
}

#' fetch a windowed version of a bam file, returns GRanges
#' @param target_strand character. if one of "+" or "-", reads are filtered
#' @param fragLen numeric, NULL, or NA.  if numeric, supplied value is used.
#' if NULL, value is calculated with fragLen_calcStranded
#' if NA, raw bam pileup with no cross strand shift is returned.
#' accordingly. ignored if any other value.
#' @export
fetchWindowedBam = function(bam_f, qgr, fragLen = NULL, win_size = 50, target_strand = c("*", "+", "-")[1]) {
    GRanges(fetchWindowedBam_dt(bam_f, qgr, fragLen, win_size, target_strand))
}
