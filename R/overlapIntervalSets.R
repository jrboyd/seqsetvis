#' Intersect a list of GRanges to create a single GRanges object of
#' merged ranges including metadata describing overlaps per input GRanges
#' @export
#' @param grs A list of GRanges
#' @param ext An integer specifying how far to extend ranges before merging.
#' in effect, ranges withing 2*ext of one another will be joined during the
#' merge
#' @param use_first A logical.  If True, instead of merging all grs, only use
#' first and add metadata logicals for others.
#' @return GRanges with metadata columns describing overlap of input grs
#' @examples
#' library(GenomicRanges)
#' a = GRanges("chr1", IRanges(1:7*10, 1:7*10))
#' b = GRanges("chr1", IRanges(5:10*10, 5:10*10))
#' ssvOverlapIntervalSets(list(a, b))
#' @import GenomicRanges
ssvOverlapIntervalSets = function(grs, ext = 0, use_first = FALSE){
  queryHits = NULL
  if(is(grs, "GRangesList")){
    grs = as.list(grs)
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
      olaps = findOverlaps(base_gr, grs[[i]])
      mcols(base_gr)[[nam]][queryHits(olaps)] = TRUE
    }
    if(use_first){
      base_gr$no_hit = apply(mcols(base_gr), 1, function(x)all(!x))
    }
    # base_gr$group = "no_hit"
    # for(i in rev(seq_along(grs))){
    #   i_gr = grs[[i]]
    #   olaps = findOverlaps(base_gr, i_gr)
    #   base_gr[queryHits(olaps)]$group = names(grs)  [i]
    # }
  })
  names(base_gr) = seq_along(base_gr)
  return(base_gr)
}


