#' Intersect a list of GRanges to create a single GRanges object of
#' merged ranges including metadata describing overlaps per input GRanges
#'
#' @export
#' @param grs A list of GRanges
#' @param ext An integer specifying how far to extend ranges before merging.
#' in effect, ranges withing 2*ext of one another will be joined during the
#' merge
#' @param use_first A logical.  If True, instead of merging all grs, only use
#' first and add metadata logicals for others.
#' @param preserve_mcols Controls carrying forward mcols metadata from input
#'   list of GRanges. If TRUE, all mcols will be carried forward with the item
#'   name appended. If a character vector, only those attributes will be carried
#'   and all must be present in all GRanges. The default of FALSE will carry
#'   nothing forward and only membership table will be generated.
#'   [ssvAnnotateSubjectGRanges] is used internally.
#' @param ... arguments passed to IRanges::findOverlaps, i.e. maxgap, minoverlap, type, select, invert.
#' @return GRanges with metadata columns describing overlap of input grs.
#' @examples
#' library(GenomicRanges)
#' a = GRanges("chr1", IRanges(1:7*10, 1:7*10))
#' b = GRanges("chr1", IRanges(5:10*10, 5:10*10))
#' ssvOverlapIntervalSets(list(a, b))
#' @import GenomicRanges
ssvOverlapIntervalSets = function(grs, ext = 0, use_first = FALSE, preserve_mcols = FALSE, ...){
  queryHits = NULL
  if(is(grs, "GRangesList")){
    grs = as.list(grs)
  }
  if(!is.list(grs)){
    stop("grs must be a List or GRangesList")
  }
  if(is.null(names(grs))){
    names(grs) = paste0("set_", LETTERS[seq_along(grs)])
  }
  if(use_first){
    base_gr = grs[[1]]
    mcols(base_gr) = NULL
  }else{
    grs_nometa = lapply(grs, function(x){mcols(x) = NULL; x})
    base_gr = reduce(unlist(GRangesList(grs_nometa)))
    start(base_gr) = start(base_gr) - ext
    end(base_gr) = end(base_gr) + ext
    base_gr = reduce(base_gr)
    start(base_gr) = start(base_gr) + ext
    end(base_gr) = end(base_gr) - ext
  }
  suppressWarnings({
    for(i in seq_along(grs)){
      nam = names(grs)[i]
      mcols(base_gr)[[nam]] = FALSE
      olaps = findOverlaps(base_gr, grs[[i]], ...)
      mcols(base_gr)[[nam]][queryHits(olaps)] = TRUE
    }
    if(use_first){
      base_gr$no_hit = apply(mcols(base_gr), 1, function(x)all(!x))
    }
  })
  names(base_gr) = seq_along(base_gr)
  if(preserve_mcols){
      base_gr = ssvAnnotateSubjectGRanges(grs, base_gr)
  }
  return(base_gr)
}

#' Intersect a list of GRanges to create a single GRanges object of merged
#' ranges including metadata describing overlaps per input GRanges.
#'
#' In constrast to ssvOverlapIntervalSets, only regions where a consensus of
#' input grs are present are preserved and annotated.
#'
#' @param grs A list of GRanges
#' @param ext An integer specifying how far to extend ranges before merging. in
#'   effect, ranges withing 2*ext of one another will be joined during the merge
#' @param min_number An integer number specifying the absloute minimum of input
#'   grs that must overlap for a site to be considered consensus.
#' @param min_fraction A numeric between 0 and 1 specifying the fraction of grs
#'   that must overlap to be considered consensus.
#' @param preserve_mcols Controls carrying forward mcols metadata from input
#'   list of GRanges. If TRUE, all mcols will be carried forward with the item
#'   name appended. If a character vector, only those attributes will be carried
#'   and all must be present in all GRanges. The default of FALSE will carry
#'   nothing forward and only membership table will be generated.
#'   [ssvAnnotateSubjectGRanges] is used internally.
#' @param ... arguments passed to IRanges::findOverlaps, i.e. maxgap, minoverlap, type, select, invert.
#' @details Only the most stringent of min_number or min_fraction will be
#'   applied.
#' @return GRanges with metadata columns describing consensus overlap of input
#'   grs.
#' @importFrom GenomeInfoDb seqlengths
#' @export
#'
#' @examples
#' library(GenomicRanges)
#' a = GRanges("chr1", IRanges(1:7*10, 1:7*10))
#' b = GRanges("chr1", IRanges(5:10*10, 5:10*10))
#' ssvConsensusIntervalSets(list(a, b))
ssvConsensusIntervalSets = function(grs, ext = 0, min_number = 2, min_fraction = .5, preserve_mcols = FALSE, ...){
  if(is(grs, "GRangesList")){
    grs = as.list(grs)
  }
  if(!is.list(grs)){
    stop("grs must be a List or GRangesList")
  }
  if(is.null(names(grs))){
    names(grs) = paste0("set_", LETTERS[seq_along(grs)])
  }
  stopifnot(min_fraction >= 0 & min_fraction <= 1)
  stopifnot(min_number >= 0 & min_number <= length(grs))
  min_by_fraction = ceiling(length(grs) * min_fraction)
  min_consensus = 0
  if(min_by_fraction > min_number){
    message("based on min_fraction of ", min_fraction, " consensus threshold is ", min_by_fraction, ".")
    min_consensus = min_by_fraction
  }else{
    message("min_number of ", min_number, " used for consensus threshold.")
    min_consensus = min_number
  }

  if(min_consensus > length(grs)){
    stop("Consensus threshold of ", min_consensus, " is too high for length of input grs, ", length(grs))
  }

  grs_ext = lapply(grs, function(gr){
    resize(gr, width(gr)+2*ext, fix = "center")
  })

  grs_cov = GRanges(coverage(unlist(GRangesList(grs_ext))) >= max(min_consensus, 1))
  grs_cov = subset(grs_cov, score == TRUE)
  grs_cov$score = NULL
  for(nam in names(grs)){
    gr = grs[[nam]]
    olaps = findOverlaps(grs_cov, gr, ...)
    mcols(grs_cov)[[nam]] = FALSE
    mcols(grs_cov)[[nam]][queryHits(olaps)] = TRUE
  }
  GenomeInfoDb::seqlengths(grs_cov) = NA
  if(preserve_mcols){
      grs_cov = ssvAnnotateSubjectGRanges(annotation_source = grs, subject_gr = grs_cov)
  }
  grs_cov
}

