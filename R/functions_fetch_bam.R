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
#' @param file_attribs optional data.frame/data.table with one row per item in
#' file paths.  Each column will be a variable added to final tidy output.
#' @param return_data.table logical. If TRUE the internal data.table is returned
#'   instead of GRanges.  Default is FALSE.
#' @param max_dupes numeric >= 1.  duplicate reads by strandd start position
#'   over this number are removed, Default is Inf.
#' @param splice_strategy character, one of c("none", "ignore", "add", "only",
#'   "splice_count"). Default is "none" and spliced alignment are asssumed not
#'   present. fragLen will be forced to be NA for any other value.  "ignore" will
#'   not count spliced regions.  add" counts spliced regions along with others,
#'   "only" will only count spliced regions and ignore others.
#' @param n_cores integer number of cores to use.
#' Uses mc.cores option if not supplied.
#' @param n_region_splits integer number of splits to apply to qgr. The query
#'   GRanges will be split into this many roughly equal parts for increased
#'   parallelization. Default is 1, no split.
#' @param return_unprocessed boolean. if TRUE returns read alignment in data.table. Default is FALSE.
#' @param force_skip_centerFix boolean, if TRUE all query ranges will be
#' used "as is".  This is already the case by default if win_method == "summary"
#' but may have applications where win_method == "sample".
#' @param ... passed to Rsamtools::ScanBamParam()
#' @return A tidy formatted GRanges (or data.table if specified) containing
#'   fetched values.
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#'   10, values from positions chr1 5,15,25...85, and 95 will be retrieved from
#'   \code{bw_file}
#' @importFrom stats weighted.mean
#' @importFrom pbapply pblapply
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
                       unique_names = NULL,
                       names_variable = "sample",
                       file_attribs = NULL,
                       win_size = 50,
                       win_method = c("sample", "summary")[1],
                       summary_FUN = stats::weighted.mean,
                       fragLens = "auto",
                       target_strand = c("*", "+", "-", "both")[1],
                       flip_strand = FALSE,
                       anchor = c("left", "left_unstranded", "center",
                                  "center_unstranded")[3],
                       return_data.table = FALSE,
                       max_dupes = Inf,
                       splice_strategy = c("none", "ignore", "add",
                                           "only", "splice_count")[1],
                       n_cores = getOption("mc.cores", 1),
                       n_region_splits = 1,
                       return_unprocessed = FALSE,
                       force_skip_centerFix = FALSE,
                       ...){
    stopifnot(all(is.character(fragLens) |
                      is.numeric(fragLens) |
                      is.na(fragLens)))
    exp_fragLen = ifelse(is.data.frame(file_paths) || is.data.table(file_paths), nrow(file_paths), length(file_paths))
    stopifnot(length(fragLens) == 1 || length(fragLens) == exp_fragLen)
    if(splice_strategy != "none") fragLens = NA
    if(return_unprocessed) fragLens = NA
    if(length(fragLens) == 1){
        if (is.data.frame(file_paths) || is.data.table(file_paths)) {
            fragLens = rep(fragLens[1], nrow(file_paths))
        }else{
            fragLens = rep(fragLens[1], length(file_paths))
        }
    }
    tmp = .get_file_attribs(file_paths, file_attribs)
    file_paths = tmp$file_paths
    file_attribs = tmp$file_attribs
    remove(tmp)
    unique_names = .get_unique_names(unique_names, file_paths, file_attribs, names_variable)

    names(fragLens) = unique_names
    names(file_paths) = unique_names

    if(any(fragLens == 'auto', na.rm = TRUE)){
        message("Calculating average fragment length for ", sum(fragLens == "auto", na.rm = TRUE), " bam file(s).")
        calc_fl = ssv_mclapply(names(fragLens), function(nam){
            if(!is.na(fragLens[nam]) && fragLens[nam] == "auto"){
                fl = fragLen_calcStranded(file_paths[nam], qgr, flip_strand = flip_strand)
                message("fragLen for ", basename(nam), " was calculated as: ", fl)
            }else{
                fl = fragLens[nam]
            }
            fl
        })

        calc_fl = unlist(calc_fl)
        names(calc_fl) = names(fragLens)
        fragLens = calc_fl
    }

    fragLens = as.numeric(fragLens)
    names(fragLens) = unique_names

    load_bam = function(f, nam, qgr) {
        message("loading ", nam, " ...")
        if(!file.exists(paste0(f, ".bai"))){
            warning("creating index for ", f)
            Rsamtools::indexBam(f)
        }
        fl = fragLens[nam]
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
                                flip_strand = flip_strand,
                                return_unprocessed = return_unprocessed,
                                force_skip_centerFix = force_skip_centerFix,
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
                         file_attribs = file_attribs,
                         win_size = win_size,
                         win_method = win_method,
                         return_data.table = TRUE,
                         n_cores = n_cores,
                         n_region_splits = n_region_splits,
                         force_skip_centerFix = force_skip_centerFix)

    if(!return_data.table & !return_unprocessed){
        bdt = GRanges(bdt)
    }
    bdt
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
#' @param flip_strand if TRUE, strand alignment is flipped prior to fragLen
#'   extension. Default is FALSE.
#' @param return_unprocessed boolean. if TRUE returns read alignment in data.table. Default is FALSE.
#' @param force_skip_centerFix boolean, if TRUE all query ranges will be
#' used "as is".  This is already the case by default if win_method == "summary"
#' but may have applications where win_method == "sample".
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
                              flip_strand = FALSE,
                              return_unprocessed = FALSE,
                              force_skip_centerFix = FALSE,
                              ...) {
    x = id = y = NULL #binding for data.table
    stopifnot(is.character(win_method))
    stopifnot(length(win_method) == 1)
    stopifnot(is(qgr, "GRanges"))
    stopifnot(win_method %in% c("sample", "summary"))
    stopifnot(is.function(summary_FUN))
    stopifnot(target_strand %in% c("*", "+", "-", "both"))
    stopifnot(anchor %in% c("left", "left_unstranded", "center",
                            "center_unstranded"))
    stopifnot(splice_strategy %in% c("none", "ignore", "add",
                                     "only", "splice_count"))
    if(return_unprocessed){
        score_gr = fetchBam(bam_f,
                            qgr,
                            fragLen = NA,
                            target_strand = "*",
                            max_dupes,
                            splice_strategy = "none",
                            return_unprocessed = return_unprocessed,
                            ...)
        tmp = qgr
        strand(tmp) = "*"
        score_gr = merge(score_gr, data.table(which_label = as.character(tmp), id = names(qgr)), by = "which_label")
        return(score_gr)
    }
    if(splice_strategy == "splice_count"){
        return(fetchBam(bam_f,
                        qgr,
                        fragLen = NA,
                        target_strand = target_strand,
                        max_dupes, splice_strategy,
                        flip_strand = flip_strand,
                        ...))
    }

    #splitting by qgr strand is necessary for strand sensitive fetch when
    #features on opposing strands overlap
    strands_TODO = unique(as.character(strand(qgr)))
    out_l = list()
    for(s in strands_TODO){
        strand_qgr = subset(qgr, strand == s)
        switch (
            win_method,
            sample = {
                strand_qgr = prepare_fetch_GRanges_width(strand_qgr, win_size, skip_centerFix = force_skip_centerFix)
                if(target_strand == "both"){
                    pos_gr = fetchBam(bam_f, strand_qgr, fragLen, "+",
                                      max_dupes, splice_strategy,
                                      flip_strand = flip_strand, ...)
                    neg_gr = fetchBam(bam_f, strand_qgr, fragLen, "-",
                                      max_dupes, splice_strategy,
                                      flip_strand = flip_strand, ...)
                    pos_dt = viewGRangesWinSample_dt(pos_gr, strand_qgr,
                                                     win_size, anchor = anchor)
                    neg_dt = viewGRangesWinSample_dt(neg_gr, strand_qgr,
                                                     win_size, anchor = anchor)
                    pos_dt[, strand := "+"]
                    neg_dt[, strand := "-"]
                    out = rbind(
                        pos_dt,
                        neg_dt
                    )
                }else{
                    score_gr = fetchBam(bam_f, strand_qgr, fragLen, target_strand,
                                        max_dupes, splice_strategy,
                                        flip_strand = flip_strand, ...)
                    out = viewGRangesWinSample_dt(score_gr, strand_qgr,
                                                  win_size, anchor = anchor)
                    out[, strand := target_strand]
                }


            },
            summary = {
                if(target_strand == "both"){
                    pos_gr = fetchBam(bam_f, strand_qgr, fragLen, "+",
                                      max_dupes, splice_strategy,
                                      flip_strand = flip_strand, ...)
                    neg_gr = fetchBam(bam_f, strand_qgr, fragLen, "-",
                                      max_dupes, splice_strategy,
                                      flip_strand = flip_strand, ...)
                    pos_dt = viewGRangesWinSummary_dt(pos_gr, strand_qgr, win_size,
                                                      summary_FUN = summary_FUN,
                                                      anchor = anchor)

                    neg_dt = viewGRangesWinSummary_dt(neg_gr, strand_qgr, win_size,
                                                      summary_FUN = summary_FUN,
                                                      anchor = anchor)
                    pos_dt[, strand := "+"]
                    neg_dt[, strand := "-"]
                    out = rbind(
                        pos_dt,
                        neg_dt
                    )
                }else{
                    score_gr = fetchBam(bam_f, strand_qgr, fragLen, target_strand,
                                        max_dupes, splice_strategy,
                                        flip_strand = flip_strand, ...)
                    out = viewGRangesWinSummary_dt(score_gr, strand_qgr, win_size,
                                                   summary_FUN = summary_FUN,
                                                   anchor = anchor)
                    out[, strand := target_strand]
                }
            }
        )
        out_l[[s]] = out
    }
    out = rbindlist(out_l)
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
#' @param flip_strand if TRUE, strand alignment is flipped prior to fragLen
#'   extension. Default is FALSE.
#' @param return_unprocessed boolean. if TRUE returns read alignment in data.table. Default is FALSE.
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
                    flip_strand = FALSE,
                    return_unprocessed = FALSE,
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
    qgr = harmonize_seqlengths(qgr, bam_f)
    if(is.null(fragLen)){
        fragLen = fragLen_calcStranded(bam_f, qgr, flip_strand = flip_strand)
        message("fragLen was calculated as: ", fragLen)
    }
    if(!is.na(fragLen)){
        stopifnot(is.numeric(fragLen))
        stopifnot(fragLen %% 1 == 0)
        stopifnot(fragLen >= 1)
    }
    sbgr = qgr
    strand(sbgr) = "*"

    if(return_unprocessed){
        sbParam = Rsamtools::ScanBamParam(
            which = sbgr,
            what = Rsamtools::scanBamWhat(),
            ...)
        bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
        bam_dt = lapply(bam_raw, function(x){
            data.table(seqnames = x$rname,
                       strand = x$strand,
                       start = x$pos,
                       width = x$qwidth,
                       cigar = x$cigar,
                       qname = x$qname,
                       flag = x$flag,
                       mapq = x$mapq,
                       mrnm = x$mrnm,
                       mpos = x$mpos,
                       isize = x$isize,
                       seq = as.character(x$seq),
                       qual = as.character(x$qual))
        })
        bam_dt = data.table::rbindlist(bam_dt,
                                       use.names = TRUE,
                                       idcol = "which_label")
        return(bam_dt)
    }else{
        sbParam = Rsamtools::ScanBamParam(
            which = sbgr,
            what = c("rname", "strand", "pos", "qwidth", "cigar"),
            ...)
        bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
        bam_dt = lapply(bam_raw, function(x){
            data.table(seqnames = x$rname, strand = x$strand,
                       start = x$pos, width = x$qwidth, cigar = x$cigar)
        })
        bam_dt = data.table::rbindlist(bam_dt,
                                       use.names = TRUE,
                                       idcol = "which_label")
    }


    bam_dt = bam_dt[!is.na(width)]
    #as.character for GRanges does a couple annoying things. for GRanges
    #of width 1, end is omitted.  Appending strand also not desirable here.
    # want chr1:17508106-17508106 not chr1:17508106:-
    gr_as.character = function(gr){
        if(length(gr) == 0){
            character()
        }else{
            paste0(seqnames(gr), ":", start(gr), "-", end(gr))
        }
    }

    toflip = gr_as.character(subset(qgr, strand == "-"))
    if(target_strand %in% c("+", "-")){
        if(!flip_strand){
            bam_dt = bam_dt[strand == target_strand & !(which_label %in% toflip) |
                                strand != target_strand & which_label %in% toflip]
        }else{
            bam_dt = bam_dt[strand != target_strand & !(which_label %in% toflip) |
                                strand == target_strand & which_label %in% toflip]
            bam_dt[, strand := ifelse(strand == "+", "-", "+")]
        }

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
        if(flip_strand){
            bam_dt[strand == "-", end := start + as.integer(fragLen) - 1L]
            bam_dt[strand == "+", start := end - as.integer(fragLen) + 1L]
        }else{
            bam_dt[strand == "+", end := start + as.integer(fragLen) - 1L]
            bam_dt[strand == "-", start := end - as.integer(fragLen) + 1L]
        }

    }
    if(splice_strategy == "splice_count"){
        sp_dt = bam_dt[, .N,
                       by = list(which_label, seqnames, start, end, strand)]
        if(target_strand %in% c("+", "-")){
            sp_dt = sp_dt[strand == target_strand]
        }
        return(sp_dt)
    }
    ext_cov = coverage(split(GRanges(bam_dt), bam_dt$which_label))
    score_gr = GRanges(ext_cov)
    score_gr = harmonize_seqlengths(score_gr, bam_f)

    if(length(score_gr) == 0){
        score_gr = sbgr
        mcols(score_gr) = NULL
        score_gr$score = 0
    }
    if(target_strand == "+"){
        strand(score_gr) = "+"
    }
    if(target_strand == "-"){
        strand(score_gr) = "-"
    }
    return(score_gr)
}

#' harmonize_seqlengths
#'
#' ensures compatibility between seqlength of gr and bam_file based on header
#'
#' @param query_gr GRanges, object to harmonize seqlengths for
#' @param bam_file character, a path to a valid bam file
#' @param force_fix Logical, if TRUE incompatible seqnames are removed from the query_gr.  Default is FALSE.
#'
#' @return GRanges with seqlengths matching bam_file
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Rsamtools scanBamHeader
#' @export
#' @examples
#' library(GenomicRanges)
#' query_gr = GRanges("chr1", IRanges(1, 100))
#' #seqlengths has not been set
#' seqlengths(query_gr)
#' bam = system.file("extdata/test.bam", package = "seqsetvis")
#' gr2 = harmonize_seqlengths(query_gr, bam)
#' #seqlengths now set
#' seqlengths(gr2)
harmonize_seqlengths = function(query_gr, bam_file, force_fix = FALSE){
    chr_lengths = Rsamtools::scanBamHeader(bam_file)[[1]]$targets
    if(!all(as.character(seqnames(query_gr)) %in% names(chr_lengths))){
        if(force_fix){
            len.orig = length(query_gr)
            query_gr = subset(query_gr, seqnames %in% names(chr_lengths))
            len.new = length(query_gr)
            message(len.orig - len.new, " regions were removed due to force_fix = TRUE where seqnames were missing from bam file.")
        }else{
        stop("There are chromosomes present in query GRanges (",
             paste(setdiff(as.character(seqnames(query_gr)), names(chr_lengths)), collapse = ", "),
                   ") not present in bam file (", bam_file, ")!\n",
             "If you want to proceed, run harmonize_seqlengths with force_fix = TRUE on query_gr to remove incompatible seqnames and retry.")
        }
    }
    to_remove =setdiff(names(GenomeInfoDb::seqlengths(query_gr)), unique(as.character(seqnames(query_gr))))
    if(length(to_remove) > 0){
        query_gr = GenomeInfoDb::dropSeqlevels(query_gr, to_remove)
    }
    GenomeInfoDb::seqlengths(query_gr) =
        chr_lengths[names(GenomeInfoDb::seqlengths(query_gr))]
    too_long = end(query_gr) >
        GenomeInfoDb::seqlengths(query_gr)[as.character(seqnames(query_gr))]
    if(any(too_long)){
        message(sum(too_long),
                " region shifted for extending beyond seqlengths")
        fix_gr = query_gr[too_long]
        shift_by = -(end(fix_gr) - GenomeInfoDb::seqlengths(fix_gr)[
            as.character(seqnames(fix_gr))])
        query_gr[too_long] = GenomicRanges::shift(fix_gr, shift_by)
    }
    too_short = start(query_gr) < 1
    if(any(too_short)){
        message(sum(too_short),
                " region shifted for starting before seqlengths")
        fix_gr = query_gr[too_short]
        shift_by = 1 - start(fix_gr)
        query_gr[too_short] = GenomicRanges::shift(fix_gr, shift_by)
    }
    query_gr
}

#' determine the most common read length for input bam_file.  uses 50 randomly
#' selected regions from query_gr.  If fewer than 20 reads are present, loads
#' all of query_gr.
#'
#' @param bam_file indexed bam file
#' @param query_gr GRanges to read from bam file
#' @export
#' @return numeric of most common read length.
#' @examples
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' bam_file = system.file("extdata/test.bam", package = "seqsetvis", mustWork = TRUE)
#' getReadLength(bam_file, qgr)
getReadLength = function(bam_file,
                         query_gr){
    Param <- Rsamtools::ScanBamParam(which=sample(query_gr,
                                                  min(50, length(query_gr))),
                                     what=c("flag","mapq"))
    temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
    if(length(temp) < 20){
        Param <- Rsamtools::ScanBamParam(which=sample(query_gr),
                                         what=c("flag","mapq"))
        temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
    }
    if(length(temp) == 0){
        warning("No reads could be found in query_gr to determine read length. Setting to 180.",
                "Please verify that query_gr is appropriate for bam_file: ", bam_file)
        readlength = 180
    }else{
        readlength=as.numeric(names(sort(table(width(temp)), decreasing = TRUE))[1])
    }
    readlength
}


#' Expand cigar codes to GRanges
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
#' @param return_data.table if TRUE, a data.table is returned, else a GRanges.
#' Default is FALSE.
#' @export
#' @return data.table with cigar entries expanded
#' @examples
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' bam_file = system.file("extdata/test.bam", package = "seqsetvis", mustWork = TRUE)
#' raw_dt = ssvFetchBam(bam_file, qgr, return_unprocessed = TRUE)
#' expandCigar(raw_dt)
expandCigar = function(cigar_dt, op_2count = c("M", "D", "=", "X"), return_data.table = FALSE){
    if("GRanges" %in% class(cigar_dt)){
        cigar_dt = as.data.table(cigar_dt)
    }
    stopifnot(c("which_label", "seqnames", "strand", "start", "cigar") %in% colnames(cigar_dt))
    res = .expand_cigar_dt(cigar_dt, op_2count)
    if(!return_data.table){
        res = GRanges(res)
    }
    return(res)
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

