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

#' Calculate cross correlation by using shiftApply on read coverage Rle
#'
#' @param bam_file character. Path to .bam file, must have index at .bam.bai.
#' @param query_gr GRanges.  Regions to calculate cross correlation for.
#' @param max_dupes integer.  Duplicate reads above this value will be removed.
#' @param fragment_sizes integer.  fragment size range to search for maximum
#'   correlation.
#' @param read_length integer. Any values outside fragment_range that must be
#'   searched.  If not supplied will be determined from bam_file.  Set as NA
#'   to disable this behavior.
#' @param ... arguments passed to ScanBamParam
#' @return named list of results
#' @export
#' @import GenomicRanges GenomicAlignments pbapply Rsamtools S4Vectors
#' @examples
#' bam_f = system.file("extdata/test.bam",
#'     package = "seqsetvis", mustWork = TRUE)
#' query_gr = CTCF_in_10a_overlaps_gr[1:2]
#' crossCorrByRle(bam_f, query_gr[1:2], fragment_sizes = seq(50, 300, 50))
crossCorrByRle = function(bam_file,
                          query_gr,
                          max_dupes = 1,
                          fragment_sizes = 50:300,
                          read_length = NULL,
                          ...){
    rn = NULL # reserve for data.table
    if(is.null(query_gr$name)){
        if(is.null(names(query_gr))){
            query_gr$name = paste0("peak_", seq_along(query_gr))
        }else{
            query_gr$name = names(query_gr)
        }
    }else{
        if(is.null(names(query_gr))){
            names(query_gr) = query_gr$name
        }else{
            #both names() and $name are set, leave it alone
        }

    }
    q_widths = apply(cbind(width(query_gr), max(fragment_sizes)*3), 1, max)
    query_gr = resize(query_gr, q_widths, fix = "center")
    names(query_gr) = query_gr$name
    # query_gr = resize(query_gr, 500, fix = "center")

    query_gr = harmonize_seqlengths(query_gr, bam_file)

    Param <- ScanBamParam(
        which=query_gr,
        what=c("flag","mapq"), ...)
    temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
    dt = as.data.table(temp)
    if(is.null(read_length)){
        read_length = getReadLength(bam_file, query_gr)
    }
    if(is.na(read_length)){
        read_length = numeric()
    }
    fragment_sizes = sort(union(read_length, fragment_sizes))

    PosCoverage <- coverage(GenomicRanges::shift(GRanges(temp[strand(temp)=="+"])), -read_length)
    PosCoverage = PosCoverage[query_gr]
    names(PosCoverage) = query_gr$name

    NegCoverage <- coverage(GRanges(temp[strand(temp)=="-"]))
    NegCoverage = NegCoverage[query_gr]
    names(NegCoverage) = query_gr$name
    # ShiftMatCor = vapply(as.list(seq_along(query_gr)), FUN.VALUE = as.list(seq(fragment_sizes)),
    #                      function(i){
    #     ShiftsCorTemp <- S4Vectors::shiftApply(fragment_sizes,
    #                                            PosCoverage[[i]],
    #                                            NegCoverage[[i]],
    #                                            cor, simplify = FALSE,
    #                                            verbose = FALSE)
    # })

    ShiftMatCor = vapply((seq_along(query_gr)), FUN.VALUE = numeric(length(fragment_sizes)),
                         function(i){
                             ShiftsCorTemp <- S4Vectors::shiftApply(fragment_sizes,
                                                                    PosCoverage[[i]],
                                                                    NegCoverage[[i]],
                                                                    cor, simplify = FALSE,
                                                                    verbose = FALSE)
                             unlist(ShiftsCorTemp)
                         })
    #necessary due to singleton query_gr or shift not resulting in matrix
    # ShiftMatCor = matrix(unlist(ShiftMatCor),
    #                      byrow = FALSE,
    #                      nrow = length(fragment_sizes),
    #                      ncol = length(query_gr))
    ShiftMatCor[is.nan(ShiftMatCor)] = 0

    colnames(ShiftMatCor) = query_gr$name
    rownames(ShiftMatCor) = fragment_sizes
    shift_dt = as.data.table(ShiftMatCor, keep.rownames = TRUE)
    shift_dt[, shift := as.numeric(rn)]
    shift_dt$rn = NULL
    shift_dt = melt(shift_dt, id.vars = "shift",
                    variable.name = "id", value.name = "correlation")
    return(shift_dt)
}

#' harmonize_seqlengths
#'
#' ensures compatibility between seqlength of gr and bam_file based on header
#'
#' @param gr GRanges, object to harmonize seqlengths for
#' @param bam_file character, a path to a valid bam file
#'
#' @return gr with seqlengths matching bam_file
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools scanBamHeader
#' @export
#' @examples
#' library(GenomicRanges)
#' gr = GRanges("chr1", IRanges(1, 100))
#' #seqlengths has not been set
#' seqlengths(gr)
#' bam = system.file("extdata/test.bam", package = "seqsetvis")
#' gr2 = harmonize_seqlengths(gr, bam)
#' #seqlengths now set
#' seqlengths(gr2)
harmonize_seqlengths = function(gr, bam_file){
    chr_lengths = Rsamtools::scanBamHeader(bam_file)[[1]]$targets
    GenomeInfoDb::seqlengths(gr) =
        chr_lengths[names(GenomeInfoDb::seqlengths(gr))]
    too_long = end(gr) >
        GenomeInfoDb::seqlengths(gr)[as.character(seqnames(gr))]
    if(any(too_long)){
        message(sum(too_long),
                " region shifted for extending beyond seqlengths")
        fix_gr = gr[too_long]
        shift_by = -(end(fix_gr) - GenomeInfoDb::seqlengths(fix_gr)[
            as.character(seqnames(fix_gr))])
        gr[too_long] = GenomicRanges::shift(fix_gr, shift_by)
    }
    too_short = start(gr) < 1
    if(any(too_short)){
        message(sum(too_short),
                " region shifted for starting before seqlengths")
        fix_gr = gr[too_short]
        shift_by = 1 - start(fix_gr)
        gr[too_short] = GenomicRanges::shift(fix_gr, shift_by)
    }
    gr
}

getReadLength = function(bam_file,
                         query_gr){
    Param <- Rsamtools::ScanBamParam(which=sample(query_gr,
                                                  min(10, length(query_gr))),
                                     what=c("flag","mapq"))
    temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
    readlength=as.numeric(names(sort(table(width(temp)), decreasing = TRUE))[1])
    readlength
}

#' calculate fragLen from a bam file for specified regions
#' @param bam_f character or BamFile. bam file to read from.
#' .bai index file must be in same directory
#' @param qgr GRanges.  used as which for ScanBamParam. Can be NULL if it's
#' REALLY important to load the entire bam, force_no_which = TRUE also required.
#' @param n_regions numeric (integer) it's generally overkill to pull all
#' regions at this stage and will slow calculation down.  Default is 100.
#' @param include_plot_in_output if TRUE ouptut is a list of fragLen and a
#' ggplot showing values considered by calculation. Default is FALSE.
#' @param test_fragLen numeric.  The set of fragment lenghts to gather
#' strand cross correlation for.
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
                                n_regions = 100,
                                include_plot_in_output = FALSE,
                                test_fragLen = seq(100, 400, 5),
                                ...){

    x = y = N = correlation = id = NULL #reserve bindings for data.table
    if(n_regions < length(qgr)){
        qgr = sample(qgr, n_regions)
    }
    cc_res = crossCorrByRle(bam_f, qgr, max_dupes = 1,
                            fragment_sizes = test_fragLen, ...)
    fragLen = cc_res[, list(shift = shift[which.max(correlation)]),
                     by = list(id)][, mean(shift)]
    fragLen = round(fragLen)
    p_res = cc_res[, list(correlation = mean(correlation)), by = list(shift)]
    if(!include_plot_in_output){
        return(fragLen)
    }else{
        pdt = data.table(x = p_res$shift, raw = p_res$correlation,
                         moving_average = movingAverage(p_res$correlation))
        fragLen = pdt$x[which.max(pdt$moving_average)]
        pdt = data.table::melt(pdt, id.vars = "x",
                               variable.name = "transform", value.name = "y")
        p = ggplot(pdt) +
            labs(x = "Fragment Length", y = "strand cross correlation",
                 title = paste(basename(bam_f)),
                 subtitle = paste("Fragment Length determined by",
                                  "strand cross correlation maximization",
                                  sep = "\n")) +
            geom_line(aes(x = x, y = y, color = transform)) +
            # scale_color_manual(values = c(raw = "black",
            #                               moving_average = "red")) +
            annotate("line",
                     x = rep(fragLen, 2),
                     y = range(pdt$y),
                     color = "green") +
            annotate("label",
                     x = fragLen,
                     y = mean(range(pdt$y)),
                     label = fragLen,
                     color = "black")
        return(list(fragLen, p))
    }
}

#' Expand intermediate bam fetch by cigar codes
#'
#' see \href{https://samtools.github.io/hts-specs/SAMv1.pdf}{sam specs} for
#' cigar details
#'
#'
#' @param cigar_dt data.table with 5 required named columns in any order.
#'   c("which_label", "seqnames", "strand", "start", "cigar")
#'
#' @return data.table with cigar entries expanded
.expand_cigar_dt_recursive = function(cigar_dt){
    cigar_w = cigar = cigar_type = tmp = which_label = NULL #dt bindings
    stopifnot(all(c("which_label", "seqnames",
                    "strand", "start",
                    "cigar") %in% colnames(cigar_dt)))
    reg1 = regexpr("[MIDNSHP=X]", cigar_dt$cigar)
    cigar_dt[, cigar_w := as.integer(substr(cigar, 1, reg1-1)) - 1L]
    cigar_dt[, cigar_type :=
                 substr(cigar, reg1, reg1+attr(reg1, "match.length")-1L)]
    cigar_dt[, tmp := substring(cigar, reg1 + 1)]
    # these do not consume reference
    cigar_dt[cigar_type %in% c("I", "S", "H", "P"), cigar_w := 0L]


    ass1 = cigar_dt[, list(which_label, seqnames, strand,
                           start, end  = start + cigar_w, cigar_type)]


    next_dt = cigar_dt[tmp != "", list(which_label,
                                       seqnames,
                                       strand,
                                       start = end + cigar_w + 1L,
                                       end = end + cigar_w + 1L,
                                       cigar = tmp)]
    if(nrow(next_dt) > 0){
        return(rbind(ass1,
                     .expand_cigar_dt_recursive(next_dt)))
    }else{
        return(ass1)
    }
}



#' Expand intermediate bam fetch by cigar codes
#'
#' see \href{https://samtools.github.io/hts-specs/SAMv1.pdf}{sam specs} for
#' cigar details
#'
#'
#' @param cigar_dt data.table with 5 required named columns in any order.
#'   c("which_label", "seqnames", "strand", "start", "cigar")
#' @param op_2count Cigar codes to count. Default is alignment (M), deletion
#'   (D), match (=), and mismatch (X).  Other useful codes may be skipped
#'   regions for RNA splicing (N).  The locations of any insterions (I) or
#'   clipping/padding (S, H, or P) will be a single bp immediately before the
#'   interval.
#'
#' @return data.table with cigar entries expanded
.expand_cigar_dt = function(cigar_dt, op_2count = c("M", "D", "=", "X")){
    cigar_type = NULL
    cigar_dt = copy(cigar_dt)
    cigar_dt[, end := start]
    exp_dt = .expand_cigar_dt_recursive(cigar_dt)
    exp_dt[cigar_type %in% op_2count]
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
.rm_dupes = function(reads_dt, max_dupes){
    ndupe = which_label = NULL
    reads_dt[, ndupe := 1L]
    reads_dt[strand == "+",
             ndupe := seq_len(.N)[order(width, decreasing = TRUE)],
             by = list(which_label, start)]
    reads_dt[strand == "-",
             ndupe := seq_len(.N)[order(width, decreasing = TRUE)],
             by = list(which_label, end)]
    reads_dt = reads_dt[ndupe <= max_dupes]
    reads_dt$ndupe = NULL
    reads_dt
}

#' fetch a bam file pileup with the ability to consider read extension to
#' fragment size (fragLen)
#' @param bam_f character or BamFile to load
#' @param qgr GRanges regions to fetchs
#' @param fragLen numeric, NULL, or NA.  if numeric, supplied value is used. if
#'   NULL, value is calculated with fragLen_calcStranded (default) if NA, raw
#'   bam pileup with no cross strand shift is returned.
#' @param target_strand character. if one of "+" or "-", reads are filtered to
#'   match. ignored if any other value.
#' @param max_dupes numeric >= 1.  duplicate reads by strandd start position
#'   over this number are removed, Default is Inf.
#' @param splice_strategy character, one of c("none", "ignore", "add", "only").
#'   Default is "none" and split read alignments are asssumed not present.
#'   fragLen must be NA for any other value to be valid.  "ignore" will not
#'   count spliced regions.  "add" counts spliced regions along with others,
#'   "only" will only count spliced regions and ignore others.
#' @param ... passed to ScanBamParam(), can't be which or what.
#' @return GRanges containing tag pileup values in score meta column.  tags are
#'   optionally extended to fragment length (fragLen) prior to pile up.
fetchBam = function(bam_f,
                    qgr,
                    fragLen = NULL,
                    target_strand = c("*", "+", "-")[1],
                    max_dupes = Inf,
                    splice_strategy = c("none", "ignore", "add",
                                        "only", "splice_count")[1],
                    ...){
    which_label = NULL #reserve binding
    stopifnot(is.numeric(max_dupes))
    stopifnot(max_dupes >= 1)
    if(!is.na(fragLen) && splice_strategy != "none"){
        stop("fragLen must be NA if splice_strategy is not 'none'.")
    }
    if( ! splice_strategy %in% c("none", "ignore", "add",
                                 "only", "splice_count")){
        stop('splice_strategy must be one of: "none", "ignore", "add", "only"')
    }
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
        what = c("rname", "strand", "pos", "qwidth", "cigar"), ...)
    bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname, strand = x$strand,
                   start = x$pos, width = x$qwidth, cigar = x$cigar)
    })
    bam_dt = data.table::rbindlist(bam_dt,
                                   use.names = TRUE,
                                   idcol = "which_label")
    toflip = sub(":-", "", as.character(subset(qgr, strand == "-")))
    if(target_strand %in% c("+", "-")){
        bam_dt = bam_dt[strand == target_strand & !(which_label %in% toflip) |
                            strand != target_strand & which_label %in% toflip]
    }

    bam_dt[, end := start + width - 1L]

    if(max_dupes < Inf){
        bam_dt = .rm_dupes(bam_dt, max_dupes)
    }
    if(is.na(fragLen)){
        bam_dt = switch(
            splice_strategy,
            none = {bam_dt},
            ignore = {.expand_cigar_dt(bam_dt)},
            add = {.expand_cigar_dt(bam_dt,
                                    op_2count = c("M", "D", "=", "X", "N"))},
            only = {.expand_cigar_dt(bam_dt, op_2count = c("N"))},
            splice_count = {.expand_cigar_dt(bam_dt, op_2count = c("N"))}
        )
    }else{
        # bam_dt[, end := start + width - 1L]
        bam_dt[strand == "+", end := start + as.integer(fragLen) - 1L]
        bam_dt[strand == "-", start := end - as.integer(fragLen) + 1L]
    }
    if(splice_strategy == "splice_count"){
        return(bam_dt[, .N,
                      by = list(which_label, seqnames, start, end, strand)])
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
#'   sample.  If win_method is summary, this is the number of windows used
#'   (confusing, sorry).
#' @param win_method character.  one of c("sample", "summary").  Determines if
#'   \code{\link{viewGRangesWinSample_dt}} or
#'   \code{\link{viewGRangesWinSummary_dt}} is used to represent each region in
#'   qgr.
#' @param summary_FUN function.  only relevant if win_method is "summary".
#'   passed to \code{\link{viewGRangesWinSummary_dt}}.
#' @param fragLen numeric, NULL, or NA.  if numeric, supplied value is used. if
#'   NULL, value is calculated with fragLen_calcStranded if NA, raw bam pileup
#'   with no cross strand shift is returned.
#' @param target_strand character. if one of "+" or "-", reads are filtered
#'   accordingly. ignored if any other value.
#' @param anchor character, one of c("center", "center_unstranded", "left",
#'   "left_unstranded")
#' @param return_data.table logical. If TRUE the internal data.table is returned
#'   instead of GRanges.  Default is FALSE.
#' @param max_dupes numeric >= 1.  duplicate reads by strandd start position
#'   over this number are removed, Default is Inf.
#' @param splice_strategy character, one of c("none", "ignore", "add", "only",
#'   "splice_count"). Default is "none" and spliced alignment are asssumed not
#'   present. fragLen must be NA for any other value to be valid.  "ignore" will
#'   not count spliced regions.  add" counts spliced regions along with others,
#'   "only" will only count spliced regions and ignore others.
#' @param ... passed to Rsamtools::ScanBamParam()
#' @return tidy GRanges (or data.table if specified) with pileups from bam file.
#'   pileup is calculated only every win_size bp.
#' @importFrom stats weighted.mean
ssvFetchBam.single = function(bam_f,
                              qgr,
                              win_size = 50,
                              win_method = c("sample", "summary")[1],
                              summary_FUN = stats::weighted.mean,
                              fragLen = NULL,
                              target_strand = c("*", "+", "-", "both")[1],
                              anchor = c("left", "left_unstranded", "center",
                                         "center_unstranded")[3],
                              return_data.table = FALSE,
                              max_dupes = Inf,
                              splice_strategy = c("none", "ignore", "add",
                                                  "only", "splice_count")[1],
                              ...) {
    stopifnot(is.character(win_method))
    stopifnot(length(win_method) == 1)
    stopifnot(class(qgr) == "GRanges")
    stopifnot(win_method %in% c("sample", "summary"))
    stopifnot(is.function(summary_FUN))
    stopifnot(target_strand %in% c("*", "+", "-", "both"))
    stopifnot(anchor %in% c("left", "left_unstranded", "center",
                            "center_unstranded"))
    stopifnot(splice_strategy %in% c("none", "ignore", "add",
                                     "only", "splice_count"))
    if(splice_strategy == "splice_count"){
        return(fetchBam(bam_f, qgr, NA, "*",
                        max_dupes, splice_strategy, ...))
    }
    switch (
        win_method,
        sample = {
            qgr = prepare_fetch_GRanges(qgr, win_size)
            if(target_strand == "both"){
                pos_gr = fetchBam(bam_f, qgr, fragLen, "+",
                                  max_dupes, splice_strategy, ...)
                neg_gr = fetchBam(bam_f, qgr, fragLen, "-",
                                  max_dupes, splice_strategy, ...)
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
                                    max_dupes, splice_strategy, ...)
                out = viewGRangesWinSample_dt(score_gr, qgr,
                                              win_size, anchor = anchor)
                out[, strand := target_strand]
            }


        },
        summary = {
            if(target_strand == "both"){
                pos_gr = fetchBam(bam_f, qgr, fragLen, "+",
                                  max_dupes, splice_strategy, ...)
                neg_gr = fetchBam(bam_f, qgr, fragLen, "-",
                                  max_dupes, splice_strategy, ...)
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
                                    max_dupes, splice_strategy, ...)
                out = viewGRangesWinSummary_dt(score_gr, qgr, win_size,
                                               summary_FUN = summary_FUN,
                                               anchor = anchor)
                out[, strand := target_strand]
            }
        }
    )
    # if(any(strand(qgr) == "-")){
    #     toflip = names(subset(qgr, strand == "-"))
    #     out[id %in% toflip & strand != "*", strand := ifelse(strand == "+", "-", "+")]
    # }
    out = out[order(x)][order(id)][order(strand)]

    if(!return_data.table){
        out = GRanges(out)
    }
    return(out)
}

#' Iterates a character vector (ideally named) and calls
#' \code{ssvFetchBam.single} on each.  Appends grouping variable to each
#' resulting data.table and uses rbindlist to efficiently combine results
#'
#' \code{ssvFetchBam} iteratively calls \code{fetchWindowedBam.single}. See
#' \code{\link{ssvFetchBam.single}} for more info.
#' @export
#' @param file_paths character vector of file_paths to load from. Alternatively,
#' file_paths can be a data.frame or data.table whose first column is a
#' character vector of paths and additial columns will be used as metadata.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#'   interval should be identical and evenly divisible by \code{win_size}.
#' @param unique_names names to use in final data.table to designate source
#'   bigwig. Default is 'sample'
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @param win_method character.  one of c("sample", "summary").  Determines if
#'   \code{\link{viewGRangesWinSample_dt}} or
#'   \code{\link{viewGRangesWinSummary_dt}} is used to represent each region in
#'   qgr.
#' @param summary_FUN function.  only relevant if win_method is "summary".
#'   passed to \code{\link{viewGRangesWinSummary_dt}}.
#' @param fragLens numeric. The fragment length to use to extend reads.  The
#'   default value "auto" causes an automatic calculation from 100 regions in
#'   qgr.  NA causes no extension of reads to fragment size.
#' @param target_strand character. One of c("*", "+", "-"). Controls filtering
#'   of reads by strand.  Default of "*" combines both strands.
#' @param flip_strand boolean. if TRUE strands are flipped.
#' @param anchor character, one of c("center", "center_unstranded", "left",
#'   "left_unstranded")
#' @param names_variable The column name where unique_names are stored.
#' @param return_data.table logical. If TRUE the internal data.table is returned
#'   instead of GRanges.  Default is FALSE.
#' @param max_dupes numeric >= 1.  duplicate reads by strandd start position
#'   over this number are removed, Default is Inf.
#' @param splice_strategy character, one of c("none", "ignore", "add", "only",
#'   "splice_count"). Default is "none" and spliced alignment are asssumed not
#'   present. fragLen must be NA for any other value to be valid.  "ignore" will
#'   not count spliced regions.  add" counts spliced regions along with others,
#'   "only" will only count spliced regions and ignore others.
#' @param n_cores integer number of cores to use.
#' @param ... passed to Rsamtools::ScanBamParam()
#' Uses mc.cores option if not supplied.
#' @return A tidy formatted GRanges (or data.table if specified) containing
#'   fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#'   10, values from positions chr1 5,15,25...85, and 95 will be retrieved from
#'   \code{bw_file}
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
                       target_strand = c("*", "+", "-", "both")[1],
                       flip_strand = FALSE,
                       anchor = c("left", "left_unstranded", "center",
                                  "center_unstranded")[3],
                       names_variable = "sample",
                       return_data.table = FALSE,
                       max_dupes = Inf,
                       splice_strategy = c("none", "ignore", "add",
                                           "only", "splice_count")[1],
                       n_cores = getOption("mc.cores", 1),
                       ...){
    stopifnot(all(is.character(fragLens) |
                      is.numeric(fragLens) |
                      is.na(fragLens)))
    stopifnot(length(fragLens) == 1 || length(fragLens) == length(file_paths))
    if(length(fragLens == 1)){
        fragLens = rep(fragLens[1], length(file_paths))
    }
    names(fragLens) = file_paths

    load_bam = function(f, nam, qgr) {
        message("loading ", f, " ...")
        if(!file.exists(paste0(f, ".bai"))){
            warning("creating index for ", f)
            Rsamtools::indexBam(f)
        }
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
                                max_dupes = max_dupes,
                                splice_strategy = splice_strategy,
                                ...)
        # dt[[names_variable]] = rep(nam, nrow(dt))
        message("finished loading ", nam, ".")
        dt
    }

    bdt = ssvFetchSignal(file_paths = file_paths,
                         qgr = qgr,
                         load_signal = load_bam,
                         unique_names = unique_names,
                         names_variable = names_variable,
                         win_size = win_size,
                         win_method = win_method,
                         return_data.table = TRUE,
                         n_cores = n_cores)
    if(flip_strand){
        if(target_strand == "*"){
            warning('flip_strand is not compatible with target_strand="*", ignored.')
        }else{
            bdt[, strand := ifelse(strand == "+", "-", "+")]
        }
    }
    if(!return_data.table){
        bdt = GRanges(bdt)
    }
    bdt

}

#' ssvFetchBamPE
#'
#' wrapper to handle standard FR paired end data
#'
#' @param file_paths character vector of file_paths to load from. Alternatively,
#' file_paths can be a data.frame or data.table whose first column is a
#' character vector of paths and additial columns will be used as metadata.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#'   interval should be identical and evenly divisible by \code{win_size}.
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @param target_strand character. if one of "+" or "-", reads are filtered to
#'   match. ignored if any other value.
#' @param splice_strategy  character, one of c("none", "ignore", "add", "only",
#'   "splice_count"). Default is "none" and spliced alignment are asssumed not
#'   present. fragLen must be NA for any other value to be valid.  "ignore" will
#'   not count spliced regions.  add" counts spliced regions along with others,
#'   "only" will only count spliced regions and ignore others.
#' @param return_data.table logical. If TRUE the internal data.table is returned
#'   instead of GRanges.  Default is FALSE.
#' @return a GRanges (or data.table if return_data.table == TRUE)
#' @export
#'
#' @examples
#' if(Sys.info()['sysname'] != "Windows"){
#' library(GenomicRanges)
#' bam_f = system.file("extdata/test.bam",
#'     package = "seqsetvis", mustWork = TRUE)
#' bam_files = c("a" = bam_f, "b" = bam_f)
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' bw_gr = ssvFetchBamPE(bam_files, qgr, win_size = 10)
#' }
ssvFetchBamPE = function(file_paths, qgr, win_size = 50, target_strand = "both", splice_strategy = "ignore",
                         return_data.table = FALSE){
    strand(qgr) = "*"
    bam_r1 = ssvFetchBam(file_paths = file_paths, qgr = qgr, target_strand = target_strand, splice_strategy = splice_strategy,
                         return_data.table = TRUE, fragLens = NA, win_size = win_size,
                         flag = scanBamFlag(isFirstMateRead = TRUE))
    cn = colnames(bam_r1)
    bam_r1$read = "r1"
    bam_r2 = ssvFetchBam(file_paths = file_paths, qgr = qgr, target_strand = target_strand,
                         return_data.table = TRUE, fragLens = NA, win_size = win_size,
                         flag = scanBamFlag(isSecondMateRead = TRUE), flip_strand = TRUE)
    bam_r2$read = "r2"

    # ggplot(rbind(bam_r1, bam_r2), aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("read")
    # bam_r1[, strand := ifelse(strand == "+", "-", "+") ]
    # ggplot(rbind(bam_r1, bam_r2), aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("read")
    bam_dt = rbind(bam_r1, bam_r2)[, cn, with = FALSE]
    bam_dt = bam_dt[, .(y = sum(y)), by = c(cn[cn != "y"])][, cn, with = FALSE]
    # if(as.character(strand(qgr)) == "+"){
    #     bam_dt[, strand := ifelse(strand == "+", "-", "+") ]
    # }
    # ggplot(bam_dt, aes(x = x, y = y, color = strand)) + geom_path()
    if(!return_data.table){
        bam_dt = GRanges(bam_dt)
    }
    bam_dt
}
