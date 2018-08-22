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
    str = utils::read.table(macs2xls_file,
                            nrows = 30,
                            comment.char = "",
                            sep = "\n",
                            stringsAsFactors = FALSE)$V1
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
#' @param max_fragLen numeric.  The maximum fragLen to calculate for.
#' Calculation time is directly proportional to this number. Default
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
            sbParam = Rsamtools::ScanBamParam(
                what = c("rname", "strand", "pos", "qwidth"), ...)
        }else{
            stop("No qgr was set for ScanBamParam which arg.  ",
                 "This will probably be very slow and uneccessary.  ",
                 "Recall with :\nforce_no_which = TRUE\n if you're certain.")
        }
    }else{
        set.seed(0)
        sbParam = Rsamtools::ScanBamParam(
            which =  sample(qgr, min(length(qgr), n_regions)),
            what = c("rname", "strand", "pos", "qwidth"),
            ...)
    }
    bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname,
                   strand = x$strand,
                   start = x$pos,
                   width = x$qwidth)
    })
    bam_dt = data.table::rbindlist(bam_dt,
                                   use.names = TRUE,
                                   idcol = "which_label")
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

        pdt = data.table::melt(pdt, id.vars = "x",
                               variable.name = "transform", value.name = "y")
        p = ggplot(pdt) +
            labs(x = "Fragment Length", y = "% strand match",
                 title = paste(basename(bam_f)),
                 subtitle = paste("Fragment Length determined by",
                                  "strand match maximization",
                                  paste("Moving averge window",
                                        ma_distance,
                                        "applied"),
                                  sep = "\n")) +
            geom_line(aes(x = x, y = y, color = transform)) +
            scale_color_manual(values = c(raw = "black",
                                          moving_average = "red")) +
            annotate("line",
                     x = rep(fragLenMa, 2),
                     y = range(perc),
                     color = "green") +
            annotate("label",
                     x = fragLenMa,
                     y = mean(range(perc)),
                     label = fragLenMa,
                     color = "black")
        return(list(fragLenMa, p))
    }
}

#' Remove duplicate reads based on stranded start position.  This is an
#' over-simplification.  For better duplicate handling, duplicates must be
#' marked in bam and flag passed to fetchBam() ... for ScanBamParam
#'
#' flag = scanBamFlag(isDuplicate = FALSE)
#'
#' @param reads_dt data.table of reads as loaded by fetchBam
#' @param max_dupes maximum allowed positional duplicates
#'
#' @return reads_dt with duplicated reads over max_dupes removed
#'
#' @examples
.rm_dupes = function(reads_dt, max_dupes){
    N = which_label = NULL
    reads_dt[, N := 1L]
    reads_dt[strand == "+", N := seq_len(.N)[order(width, decreasing = TRUE)], by = list(which_label, start)]
    reads_dt[strand == "-", N := seq_len(.N)[order(width, decreasing = TRUE)], by = list(which_label, end)]
    reads_dt = reads_dt[N <= max_dupes]
    reads_dt$N = NULL
    reads_dt
}

#' fetch a bam file pileup with the ability to consider read extension to fragment size (fragLen)
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
fetchBam = function(bam_f,
                    qgr,
                    fragLen = NULL,
                    target_strand = c("*", "+", "-")[1],
                    max_dupes = Inf,
                    ...){
    stopifnot(is.numeric(max_dupes))
    stopifnot(max_dupes >= 1)
    if(is.null(fragLen)){
        fragLen = fragLen_calcStranded(bam_f, qgr)
        message("fragLen was calculated as: ", fragLen)
    }
    if(!is.na(fragLen)){
        stopifnot(is.numeric(fragLen))
        stopifnot(fragLen %% 1 == 0)
        stopifnot(fragLen >= 1)
    }
    sbgr = qgr
    strand(sbgr) = "*"
    sbParam = Rsamtools::ScanBamParam(
        which = sbgr,
        what = c("rname", "strand", "pos", "qwidth"), ...)
    bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname, strand = x$strand,
                   start = x$pos, width = x$qwidth)
    })
    bam_dt = data.table::rbindlist(bam_dt,
                                   use.names = TRUE,
                                   idcol = "which_label")
    bam_dt[, end := start + width - 1L]

    if(target_strand == "+"){
        bam_dt = bam_dt[strand == "+"]
    }
    if(target_strand == "-"){
        bam_dt = bam_dt[strand == "-"]
    }

    if(max_dupes < Inf){
        bam_dt = .rm_dupes(bam_dt, max_dupes)
    }

    #
    # max_dupes = 2
    # bam_dt[, N := 1L]
    # bam_dt[strand == "+", N := seq_len(.N)[order(width, decreasing = TRUE)], by = .(which_label, start)]
    # bam_dt[strand == "-", N := seq_len(.N)[order(width, decreasing = TRUE)], by = .(which_label, end)]
    # bam_dt = bam_dt[N <= max_dupes]

    if(!is.na(fragLen)){#extension to fragLen
        bam_dt[strand == "+", end := start + as.integer(fragLen) - 1L]
        bam_dt[strand == "-", start := end - as.integer(fragLen) + 1L]
    }


    ext_cov = coverage(split(GRanges(bam_dt), bam_dt$which_label))
    score_gr = GRanges(ext_cov)
    if(target_strand == "+"){
        strand(score_gr) = "+"
    }
    if(target_strand == "-"){
        strand(score_gr) = "-"
    }
    return(score_gr)
}

#' fetch a windowed version of a bam file, returns GRanges
#'
#' @param bam_f character or BamFile to load
#' @param qgr GRanges regions to fetchs
#' @param win_size numeric >=1.  pileup grabbed every win_size bp for win_method
#' sample.  If win_method is summary, this is the number of windows used
#' (confusing, sorry).
#' @param win_method character.  one of c("sample", "summary").  Determines
#' if \code{\link{viewGRangesWinSample_dt}} or
#' \code{\link{viewGRangesWinSummary_dt}} is used to represent each region in
#' qgr.
#' @param summary_FUN function.  only relevant if win_method is "summary".
#' passed to \code{\link{viewGRangesWinSummary_dt}}.
#' @param fragLen numeric, NULL, or NA.  if numeric, supplied value is used.
#' if NULL, value is calculated with fragLen_calcStranded
#' if NA, raw bam pileup with no cross strand shift is returned.
#' @param target_strand character. if one of "+" or "-", reads are filtered
#' accordingly. ignored if any other value.
#' @param anchor character, one of c("center", "center_unstranded",
#' "left", "left_unstranded")
#' @param return_data.table logical. If TRUE the internal data.table is
#' returned instead of GRanges.  Default is FALSE.
#' @param max_dupes numeric >= 1.  duplicate reads by strandd start position
#' over this number are removed, Default is Inf.
#' @return tidy GRanges (or data.table if specified) with pileups from bam
#' file.  pileup is calculated only every win_size bp.
#' @export
#' @importFrom stats weighted.mean
#' @examples
#' bam_file = system.file("extdata/test.bam",
#'     package = "seqsetvis")
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' bam_gr = ssvFetchBam.single(bam_file, qgr)
#' bam_gr = ssvFetchBam.single(bam_file, qgr, fragLen = 180,
#'     win_size = 10, target_strand = "+")
#'
#' bam_dt = ssvFetchBam.single(bam_file, qgr,
#'     return_data.table = TRUE)
ssvFetchBam.single = function(bam_f,
                              qgr,
                              win_size = 50,
                              win_method = c("sample", "summary")[1],
                              summary_FUN = stats::weighted.mean,
                              fragLen = NULL,
                              target_strand = c("*", "+", "-")[1],
                              anchor = c("left", "left_unstranded", "center",
                                         "center_unstranded")[3],
                              return_data.table = FALSE,
                              max_dupes = Inf) {
    stopifnot(is.character(win_method))
    stopifnot(length(win_method) == 1)
    stopifnot(class(qgr) == "GRanges")
    stopifnot(win_method %in% c("sample", "summary"))
    stopifnot(is.function(summary_FUN))
    switch (win_method,
            sample = {
                qgr = prepare_fetch_GRanges(qgr, win_size)
                score_gr = fetchBam(bam_f, qgr, fragLen, target_strand, max_dupes)
                out = viewGRangesWinSample_dt(score_gr, qgr, win_size, anchor = anchor)
            },
            summary = {
                score_gr = fetchBam(bam_f, qgr, fragLen, target_strand, max_dupes)
                out = viewGRangesWinSummary_dt(score_gr, qgr, win_size,
                                               summary_FUN = summary_FUN,
                                               anchor = anchor)
            }
    )


    if(!return_data.table){
        out = GRanges(out)
    }
    return(out)
}

#' Iterates a character vector (ideally named) and calls \code{ssvFetchBam.single}
#' on each.  Appends grouping variable to each resulting data.table and uses
#' rbindlist to efficiently combine results
#'
#' \code{ssvFetchBam} iteratively calls \code{fetchWindowedBam.single}.
#' See \code{\link{ssvFetchBam.single}} for more info.
#' @export
#' @param file_paths The character vector or list of paths to bigwig files to
#'  read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param unique_names names to use in final data.table to designate source
#' bigwig. Default is 'sample'
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @param win_method character.  one of c("sample", "summary").  Determines
#' if \code{\link{viewGRangesWinSample_dt}} or
#' \code{\link{viewGRangesWinSummary_dt}} is used to represent each region in
#' qgr.
#' @param summary_FUN function.  only relevant if win_method is "summary".
#' passed to \code{\link{viewGRangesWinSummary_dt}}.
#' @param fragLens numeric. The fragment length to use to extend reads.  The
#' default value "auto" causes an automatic calculation from 100 regions in
#' qgr.  NA causes no extension of reads to fragment size.
#' @param target_strand character. One of c("*", "+", "-"). Controls
#' filtering of reads by strand.  Default of "*" combines both strands.
#' @param anchor character, one of c("center", "center_unstranded",
#' "left", "left_unstranded")
#' @param names_variable The column name where unique_names are stored.
#' @param return_data.table logical. If TRUE the internal data.table is
#' returned instead of GRanges.  Default is FALSE.
#' @param max_dupes numeric >= 1.  duplicate reads by strandd start position
#' over this number are removed, Default is Inf.
#' @return A tidy formatted GRanges (or data.table if specified) containing
#' fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be
#' retrieved from \code{bw_file}
#' @importFrom stats weighted.mean
#' @examples
#' if(Sys.info()['sysname'] != "Windows"){
#' library(GenomicRanges)
#' bam_f = system.file("extdata/test.bam",
#'     package = "seqsetvis", mustWork = TRUE)
#' bam_files = c("a" = bam_f, "b" = bam_f)
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' bw_gr = ssvFetchBam(bam_files, qgr, win_size = 10)
#' bw_gr2 = ssvFetchBam(as.list(bam_files), qgr, win_size = 10)
#'
#' bw_dt = ssvFetchBam(bam_files, qgr, win_size = 10,
#'     return_data.table = TRUE)
#' }
ssvFetchBam = function(file_paths,
                       qgr,
                       unique_names = names(file_paths),
                       win_size = 50,
                       win_method = c("sample", "summary")[1],
                       summary_FUN = stats::weighted.mean,
                       fragLens = "auto",
                       target_strand = c("*", "+", "-")[1],
                       anchor = c("left", "left_unstranded", "center",
                                  "center_unstranded")[3],
                       names_variable = "sample",
                       return_data.table = FALSE,
                       max_dupes = Inf){
    stopifnot(all(is.character(fragLens) | is.numeric(fragLens) | is.na(fragLens)))
    stopifnot(length(fragLens) == 1 || length(fragLens) == length(file_paths))
    if(length(fragLens == 1)){
        fragLens = rep(fragLens[1], length(file_paths))
    }
    names(fragLens) = file_paths

    load_bam = function(f, nam, qgr) {
        message("loading ", f, " ...")
        fl = fragLens[f]
        if(!is.na(fl))
            if(fl == "auto"){
                fl = NULL
            }
        dt = ssvFetchBam.single(bam_f = f,
                                qgr = qgr,
                                win_size = win_size,
                                win_method = win_method,
                                summary_FUN = summary_FUN,
                                fragLen = fl,
                                target_strand = target_strand,
                                anchor = anchor,
                                return_data.table = TRUE,
                                max_dupes = max_dupes)
        dt[[names_variable]] = nam
        message("finished loading ", nam, ".")
        dt
    }

    ssvFetchSignal(file_paths = file_paths,
                   qgr = qgr,
                   load_signal = load_bam,
                   unique_names = unique_names,
                   names_variable = names_variable,
                   win_size = win_size,
                   win_method = win_method,
                   return_data.table = return_data.table)

}
