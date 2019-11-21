#' collapse_gr
#'
#' collapse non-contiguous regions (i.e. exons) into a contiguous coordinate
#' starting at 1.  this is strand sensitive and intended for use with
#' all exons of a single gene.
#'
#' @param genome_gr a GRanges of regions on a single chromosome.  Regions
#' are intended to be non-contiguous and may even overlap.
#'
#' @return a new GRanges object with same mcols as input with all intervals
#' starting at 1 and no empty space between syntenic regions.
#' @export
#'
#' @examples
#' library(data.table)
#' library(GenomicRanges)
#' dev_dat = data.table(seqnames = "chrTest",
#'                      transcript_id = c(1, 1,  2, 2, 3, 3, 3),
#'                      start = c(5,  30,  8, 30, 2, 30, 40),
#'                      end = c(10, 35, 15, 38, 7, 35, 45),
#'                      strand = "+")
#'
#' genome_gr = GRanges(dev_dat)
#' collapse_gr(genome_gr)
#'
#' neg_gr = genome_gr
#' strand(neg_gr) = "-"
#' collapse_gr(neg_gr)
collapse_gr = function(genome_gr){
    stopifnot("GRanges" %in% class(genome_gr))
    stopifnot(length(genome_gr) > 0)
    stopifnot(length(unique(seqnames(genome_gr))) == 1)
    stopifnot(length(unique(strand(genome_gr))) == 1)
    if(as.character(strand(genome_gr[1])) == "-"){
        tmp = mcols(genome_gr)
        genome_gr = GRanges(seqnames(genome_gr),
                IRanges(max(end(genome_gr)) - end(genome_gr) + 1,
                        max(end(genome_gr))- start(genome_gr) + 1))
        mcols(genome_gr) = tmp

    }
    strand(genome_gr) = "*"
    tx_gr = genome_gr
    rm_gr = GenomicRanges::setdiff(range(genome_gr), genome_gr)
    for(i in rev(seq_along(rm_gr))){
        irm_gr = rm_gr[i]
        to_shift = start(tx_gr) > end(irm_gr)
        tx_gr[to_shift]  = GenomicRanges::shift(tx_gr[to_shift],
                                                shift = -width(irm_gr))
    }
    tx_gr = GenomicRanges::shift(tx_gr, -min(start(tx_gr)) + 1)
    return(tx_gr)
}

#' convert_collapsed_coord
#'
#' (preliminary implementation, sub-optimal)
#'
#' see \code{\link{collapse_gr}} for explanation of intended uses. this function
#' translates all values of x from original genomic coordinates to new
#' coordinate space created by \code{\link{collapse_gr}}.
#'
#'
#' @param genome_gr non-contiguous regions to collapse a la
#'   \code{\link{collapse_gr}}
#' @param x numeric, positions within genome_gr to convert to collapsed
#'   coordinates.
#'
#' @return numeric, positions of every value of x within collapse coordinates.
#' values outside of collapsed regions (an intron or outside range) will be NA.
#' @export
#'
#' @examples
#' library(data.table)
#' library(GenomicRanges)
#' dev_dat = data.table(seqnames = "chrTest",
#'                      transcript_id = c(1, 1,  2, 2, 3, 3, 3),
#'                      start = c(5,  30,  8, 30, 2, 30, 40),
#'                      end = c(10, 35, 15, 38, 7, 35, 45),
#'                      strand = "+")
#'
#' genome_gr = GRanges(dev_dat)
#' convert_collapsed_coord(genome_gr, start(genome_gr))
#' convert_collapsed_coord(genome_gr, end(genome_gr))
convert_collapsed_coord = function(genome_gr, x){
    stopifnot("GRanges" %in% class(genome_gr))
    stopifnot(length(genome_gr) > 0)
    stopifnot(length(unique(seqnames(genome_gr))) == 1)
    stopifnot(length(unique(strand(genome_gr))) == 1)
    if(as.character(strand(genome_gr[1])) == "-"){
        x = max(end(genome_gr)) - x + 1
        tmp = mcols(genome_gr)
        genome_gr = GRanges(seqnames(genome_gr),
                            IRanges(max(end(genome_gr)) - end(genome_gr) + 1,
                                    max(end(genome_gr))- start(genome_gr) + 1))
        mcols(genome_gr) = tmp
    }
    genome_gr
    red_gr = reduce(genome_gr)
    xnew = sapply(x, function(xx){
        is_over = xx > end(red_gr)
        if(any(is_over)){
            is_in = max(which(is_over)) + 1
        }else{
            is_in = 1
        }

        xx - start(red_gr)[is_in] + sum(width(red_gr[is_over])) + 1
    })


    is_miss = sapply(x, function(xx){
        !any(start(red_gr) <= xx & end(red_gr) >= xx)
    })

    xnew[is_miss] = NA

    return(xnew)
}

