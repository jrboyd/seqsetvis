#' ssvFetchBam for paired-end ChIP-seq files.
#'
#' Iterates a character vector (ideally named) and calls
#' \code{ssvFetchBamPE.single} on each.  Appends grouping variable to each
#' resulting data.table and uses rbindlist to efficiently combine results
#'
#' #' In contrast to ssvFetchBam, extension of reads to estimated fragment size is
#' not an issue as each read pair represents a fragment of exact size.
#'
#' \code{ssvFetchBamPE} iteratively calls \code{fetchWindowedBam.single}. See
#' \code{\link{ssvFetchBamPE.single}} for more info.
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
#' @param anchor character, one of c("center", "center_unstranded", "left",
#'   "left_unstranded")
#' @param names_variable The column name where unique_names are stored.
#' @param return_data.table logical. If TRUE the internal data.table is returned
#'   instead of GRanges.  Default is FALSE.
#' @param max_dupes numeric >= 1.  duplicate reads by strandd start position
#'   over this number are removed, Default is Inf.
#' @param n_cores integer number of cores to use.
#' @param min_isize integer. Read pairs must have an isize greater than or equal to this value.  Default is 1.
#' @param max_isize integer. Read pairs must have an isize less than or equal to this value.  Default is Inf.
#' @param return_unprocessed boolean. if TRUE returns read alignment in data.table. Default is FALSE.
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
#' bam_f = system.file("extdata/Bcell_PE.mm10.bam",
#'     package = "seqsetvis", mustWork = TRUE)
#' bam_files = c("a" = bam_f, "b" = bam_f)
#' data("Bcell_peaks")
#' qgr = Bcell_peaks
#' bw_gr = ssvFetchBamPE(bam_files, qgr, win_size = 10)
#' bw_gr2 = ssvFetchBamPE(as.list(bam_files), qgr, win_size = 10)
#'
#' bw_dt = ssvFetchBamPE(bam_files, qgr, win_size = 10,
#'     return_data.table = TRUE)
#' }
ssvFetchBamPE = function(file_paths,
                         qgr,
                         unique_names = NULL,
                         win_size = 50,
                         win_method = c("sample", "summary")[1],
                         summary_FUN = stats::weighted.mean,
                         anchor = c("left", "left_unstranded", "center",
                                    "center_unstranded")[3],
                         names_variable = "sample",
                         return_data.table = FALSE,
                         max_dupes = Inf,
                         n_cores = getOption("mc.cores", 1),
                         min_isize = 1,
                         max_isize = Inf,
                         return_unprocessed = FALSE,
                         force_skip_centerFix = FALSE,
                         ...){
    load_bamPE = function(f, nam, qgr) {
        message("loading ", f, " ...")
        if(!file.exists(paste0(f, ".bai"))){
            warning("creating index for ", f)
            Rsamtools::indexBam(f)
        }
        dt = ssvFetchBamPE.single(
            bam_f = f,
            qgr = qgr,
            win_size = win_size,
            win_method = win_method,
            summary_FUN = summary_FUN,
            anchor = anchor,
            return_data.table = TRUE,
            max_dupes = max_dupes,
            min_isize = min_isize,
            max_isize = max_isize,
            return_unprocessed = return_unprocessed,
            force_skip_centerFix = force_skip_centerFix,
            ...)
        # dt[[names_variable]] = rep(nam, nrow(dt))
        message("finished loading ", nam, ".")
        dt
    }

    bdt = ssvFetchSignal(file_paths = file_paths,
                         qgr = qgr,
                         load_signal = load_bamPE,
                         unique_names = unique_names,
                         names_variable = names_variable,
                         win_size = win_size,
                         win_method = win_method,
                         return_data.table = TRUE,
                         n_cores = n_cores,
                         force_skip_centerFix = force_skip_centerFix)


    if(!return_data.table & !return_unprocessed){
        bdt = GRanges(bdt)
    }
    bdt
}

#' fetch a windowed version of a paired-end bam file, returns GRanges
#' In contrast to ssvFetchBam, extension of reads to estimated fragment size is
#' not an issue as each read pair represents a fragment of exact size.
#'
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
#' @param anchor character, one of c("center", "center_unstranded", "left",
#'   "left_unstranded")
#' @param return_data.table logical. If TRUE the internal data.table is returned
#'   instead of GRanges.  Default is FALSE.
#' @param max_dupes numeric >= 1.  duplicate reads by strandd start position
#'   over this number are removed, Default is Inf.
#' @param min_isize integer. Read pairs must have an isize greater than or equal to this value.  Default is 1.
#' @param max_isize integer. Read pairs must have an isize less than or equal to this value.  Default is Inf.
#' @param return_unprocessed boolean. if TRUE returns read alignment in data.table. Default is FALSE.
#' @param ... passed to Rsamtools::ScanBamParam()
#' @return tidy GRanges (or data.table if specified) with pileups from bam file.
#'   pileup is calculated only every win_size bp.
#' @importFrom stats weighted.mean
ssvFetchBamPE.single = function(bam_f,
                                qgr,
                                win_size = 50,
                                win_method = c("sample", "summary")[1],
                                summary_FUN = stats::weighted.mean,
                                anchor = c("left", "left_unstranded", "center",
                                           "center_unstranded")[3],
                                return_data.table = FALSE,
                                max_dupes = Inf,
                                min_isize = 1,
                                max_isize = Inf,
                                return_unprocessed = FALSE,
                                force_skip_centerFix = FALSE,
                                ...) {
    x = id = y = NULL
    stopifnot(is.character(win_method))
    stopifnot(length(win_method) == 1)
    stopifnot(is(qgr, "GRanges"))
    stopifnot(win_method %in% c("sample", "summary"))
    stopifnot(is.function(summary_FUN))
    stopifnot(anchor %in% c("left", "left_unstranded", "center",
                            "center_unstranded"))
    switch (
        win_method,
        sample = {
            qgr = prepare_fetch_GRanges(qgr, win_size, skip_centerFix = force_skip_centerFix)
            score_gr = fetchBamPE(bam_f,
                                  qgr,
                                  max_dupes,
                                  min_isize,
                                  max_isize,
                                  return_unprocessed = return_unprocessed,
                                  ...)
            if(return_unprocessed){
                score_gr = merge(score_gr, data.table(which_label = as.character(qgr), id = qgr$name), by = "which_label")
                return(score_gr)
            }else{
                out = viewGRangesWinSample_dt(score_gr,
                                              qgr,
                                              win_size,
                                              anchor = anchor)
            }
        },
        summary = {
            score_gr = fetchBamPE(bam_f,
                                  qgr,
                                  max_dupes,
                                  min_isize,
                                  max_isize,
                                  return_unprocessed = return_unprocessed,
                                  ...)
            if(return_unprocessed){
                score_gr = merge(score_gr, data.table(which_label = as.character(qgr), id = qgr$name), by = "which_label")
                return(score_gr)
            }else{
                out = viewGRangesWinSummary_dt(score_gr, qgr, win_size,
                                               summary_FUN = summary_FUN,
                                               anchor = anchor)
            }
        }
    )

    out = out[order(x)][order(id)][order(strand)]

    if(!return_data.table){
        out = GRanges(out)
    }
    return(out)
}

fetchBamPE = function(bam_f,
                      qgr,
                      max_dupes = Inf,
                      min_isize = 1,
                      max_isize = Inf,
                      return_unprocessed = FALSE,
                      ...
                      ){
    isize = paired = qname = which_label = NULL #reserve bindings
    stopifnot(is.numeric(min_isize))
    stopifnot(is.numeric(max_isize))
    stopifnot(min_isize > 0)
    stopifnot(max_isize >= min_isize)
    sbgr = qgr
    strand(sbgr) = "*"
    sbParam = Rsamtools::ScanBamParam(
        which = sbgr,
        # what = Rsamtools::scanBamWhat())
        what = c("rname", "qname", "isize", "strand", "pos", "qwidth", "cigar"),
        ...)
    bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)

    bam_dt = lapply(bam_raw, function(x){
        as.data.table(data.frame(x))
        # data.table(seqnames = x$rname, qname = x$qname, strand = x$strand,
        #            start = x$pos, width = x$qwidth, cigar = x$cigar)
    })
    bam_dt = data.table::rbindlist(bam_dt,
                                   use.names = TRUE,
                                   idcol = "which_label")
    if(return_unprocessed){
        return(bam_dt)
    }
    # ggplot(bam_dt[isize >= 0, ], aes(x = isize)) + geom_histogram() + facet_wrap("which_label")

    bam_dt = bam_dt[!is.na(qwidth)]
    bam_dt[, paired := .N == 2, .(qname, which_label)]
    bam_dt = bam_dt[paired == TRUE]
    bam_dt = bam_dt[isize >= min_isize & isize <= max_isize]

    #prepare for GRanges conversion
    bam_dt = bam_dt[, list(
        which_label,
        seqnames = rname,
        strand = "*",
        start = pos,
        width = isize,
        end = pos + isize - 1L
    )]
    if(max_dupes < Inf){
        bam_dt = .rm_dupesPE(bam_dt, max_dupes)
    }

    ext_cov = coverage(split(GRanges(bam_dt), bam_dt$which_label))
    score_gr = GRanges(ext_cov)
    score_gr
}

#' Remove duplicate paired-end reads based on start and end position.  This is
#' an over-simplification.  For better duplicate handling, duplicates must be
#' marked in bam and flag passed to fetchBamPE() ... for ScanBamParam
#'
#' flag = scanBamFlag(isDuplicate = FALSE)
#'
#' @param reads_dt data.table of reads as loaded by fetchBamPE
#' @param max_dupes maximum allowed positional duplicates
#'
#' @return reads_dt with duplicated reads over max_dupes removed
.rm_dupesPE = function(reads_dt, max_dupes){
    ndupe = which_label = NULL
    reads_dt[, ndupe := 1L]
    reads_dt[,
             ndupe := seq_len(.N),
             by = list(which_label, start, end)]
    reads_dt = reads_dt[ndupe <= max_dupes]
    reads_dt$ndupe = NULL
    reads_dt
}
