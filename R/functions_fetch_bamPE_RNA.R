#' ssvFetchBamPE.RNA
#'
#' @param file_paths character vector of file_paths to load from. Alternatively,
#'   file_paths can be a data.frame or data.table whose first column is a
#'   character vector of paths and additial columns will be used as metadata.
#' @param qgr Set of GRanges to query. For valid results the width of each
#'   interval should be identical and evenly divisible by win_size.
#' @param unique_names names to use in final data.table to designate source
#'   bigwig. Default is 'sample'
#' @param win_size The window size that evenly divides widths in qgr.
#' @param target_strand character. if one of "+" or "-", reads are filtered to
#'   match. ignored if any other value.
#' @param absolute_strand If TRUE, strandedness of `qgr` will be ignored. This
#'   is useful when creating tracks for similar.
#' @param splice_strategy character, one of c("none", "ignore", "add", "only",
#'   "splice_count"). Default is "none" and spliced alignment are asssumed not
#'   present. fragLen must be NA for any other value to be valid. "ignore" will
#'   not count spliced regions. add" counts spliced regions along with others,
#'   "only" will only count spliced regions and ignore others.
#' @param return_data.table logical. If TRUE the internal data.table is returned
#'   instead of GRanges. Default is FALSE.
#' @param win_method 	character. one of c("sample", "summary"). "sample" selects
#'   values at intervals and "summary" applies a weight mean function to all
#'   values in window.
#' @param max_dupes numeric >= 1. duplicate reads by strandd start position over
#'   this number are removed, Default is Inf.
#' @param flip_strand logical. if TRUE strands are flipped.
#' @param sum_reads logical. If true R1 and R2 reads are added together. If
#'   FALSE they are returned separately, identified by the "read" attribute.
#' @param n_cores integer number of cores to use. Uses mc.cores option if not
#'   supplied.
#' @param force_skip_centerFix boolean, if TRUE all query ranges will be used
#'   "as is". This is already the case by default if win_method == "summary" but
#'   may have applications where win_method == "sample".
#' @param n_region_splits integer number of splits to apply to qgr. The query
#'   GRanges will be split into this many roughly equal parts for increased
#'   parallelization. Default is 1, no split.
#'
#' @return A tidy formatted GRanges (or data.table if specified) containing
#'   fetched values.
#' @export
#'
#' @examples
#' library(GenomicRanges)
#' pkg_dir = system.file(package = "seqsetvis", "extdata", mustWork = TRUE)
#' bam_files_esr1 = dir(pkg_dir, pattern = "H1.+R1.ESR1_RNA.+bam$", full.names = TRUE)
#' names(bam_files_esr1) = sub("_R.+", "", basename(bam_files_esr1))
#' query_gr = GenomicRanges::GRanges("chr6:151656691-152129619:+")
#' query_gr = GenomicRanges::GRanges("chr6:152116691-152129619:+")
#'
#' strand(query_gr) = "+"
#'
#' prof_dt = ssvFetchBamPE.RNA(bam_files_esr1, query_gr, return_data.table = TRUE, win_size = 1)
#' prof_dt
ssvFetchBamPE.RNA = function(
        file_paths,
        qgr,
        unique_names = NULL,
        win_size = 50,
        target_strand = "both",
        absolute_strand = FALSE,
        splice_strategy = "ignore",
        return_data.table = FALSE,
        win_method = "sample",
        max_dupes = Inf,
        flip_strand = FALSE,
        sum_reads = TRUE,
        n_cores = getOption("mc.cores", 1),
        force_skip_centerFix = TRUE,
        n_region_splits = 1){
    y = cn = NULL #reserve bindings
    if(absolute_strand){
        strand(qgr) = "*"
    }
    bam_r1 = ssvFetchBam(
        file_paths = file_paths,
        qgr = qgr,
        unique_names = unique_names,
        target_strand = target_strand,
        splice_strategy = splice_strategy,
        return_data.table = TRUE,
        fragLens = NA,
        anchor = "left",
        win_size = win_size,
        flag = Rsamtools::scanBamFlag(isFirstMateRead = TRUE),
        win_method = win_method,
        flip_strand = !flip_strand,
        max_dupes = max_dupes,
        n_cores = n_cores,
        force_skip_centerFix = force_skip_centerFix,
        n_region_splits = n_region_splits
    )
    cn = colnames(bam_r1)
    bam_r1$read = "r1"
    bam_r2 = ssvFetchBam(
        file_paths = file_paths,
        qgr = qgr,
        unique_names = unique_names,
        target_strand = target_strand,
        splice_strategy = splice_strategy,
        return_data.table = TRUE,
        fragLens = NA,
        anchor = "left",
        win_size = win_size,
        flag = Rsamtools::scanBamFlag(isSecondMateRead = TRUE),
        win_method = win_method,
        flip_strand = flip_strand,
        max_dupes = max_dupes,
        n_cores = n_cores,
        force_skip_centerFix = force_skip_centerFix,
        n_region_splits = n_region_splits
    )
    bam_r2$read = "r2"

    if(sum_reads){
        bam_dt = rbind(bam_r1, bam_r2)[, cn, with = FALSE]
        bam_dt = bam_dt[, list(y = sum(y)), by = c(cn[cn != "y"])][, cn, with = FALSE]
    }else{
        bam_dt = rbind(bam_r1, bam_r2)
    }

    if(!return_data.table){
        bam_dt = GRanges(bam_dt)
    }
    bam_dt
}
