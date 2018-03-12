#' Profiles for 100 randomly selected regions from overlapping CTCF peaks in 10a cell ChIP-seq
#' Results from CTCF_in_10a_overlaps_gr
#'
#' A tidy data.table at window size 100 bp within 300 bp of peak center
#'  The variables are as follows:
#'
#' \itemize{
#'   \item seqnames. chromosome for GRanges compatibility
#'   \item start. start of interval
#'   \item end. end of interval
#'   \item id. unique identifier
#'   \item FE. fold-enrichment over input.
#'   \item x. bp relative to center
#'   \item sample. name of originating sample
#' }
#'
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_profiles_dt
#' @format A data frame with 53940 rows and 10 variables
NULL

#' 100 randomly selected regions from overlapping CTCF peaks in 10a cell ChIP-seq
#'
#' macs2 peak calls for CTCF ChIP-seq in a MCF10A progression model.
#' See  GEO series GSE98551 for details.
#'
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_overlaps_gr
#' @format GenomicRanges with 3 metadata columns of membership table
NULL

#' FTP URL path for vignette data.
#'
#' FE bigWig tracks for CTCF ChIP-seq in a MCF10A progression model.
#' See GEO series GSE98551 for details.
#'
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_bigWig_urls
#' @format named character vector of length 3
NULL

#' FTP URL path for vignette data.
#'
#' macs2 peak calls for CTCF ChIP-seq in a MCF10A progression model.
#' See  GEO series GSE98551 for details.
#'
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_narrowPeak_urls
#' @format named character vector of length 3
NULL

