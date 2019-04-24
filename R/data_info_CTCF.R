#' CTCF ChIP-seq in breast cancer cell lines
#'
#' @description Vignette data for seqsetvis was downloaded directly from GEO
#' series
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98551}{GSE98551}.
#' This data is CTCF ChIP-seq from a model of breast cancer progression
#' derived from the MCF10A cell line.
#'
#' @description Data from GEO series
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98551}{GSE98551}
#' is from the publication
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/28504305}{Fritz AJ et al. 2018}
#' @details Contains: \itemize{
#' \item \code{\link{CTCF_in_10a_overlaps_gr}}
#' \item \code{\link{CTCF_in_10a_profiles_dt}}
#' \item \code{\link{CTCF_in_10a_bigWig_urls}}
#' \item \code{\link{CTCF_in_10a_narrowPeak_urls}}
#' }
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_data
NULL

#' Profiles for 100 randomly selected regions from overlapping CTCF peaks in
#' 10a cell ChIP-seq
#' Results from CTCF_in_10a_overlaps_gr
#'
#' A tidy GRanges at window size 50 bp within 350 bp of peak center
#'  The variables are as follows:
#' @details part of \code{\link{CTCF_in_10a_data}}
#' \enumerate{
#'   \item id. unique identifier
#'   \item y. fold-enrichment over input.
#'   \item x. bp relative to center
#'   \item sample. name of originating sample
#' }
#'
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_profiles_gr
#' @format A tidy GRanges of 2100 rows and 4 metadata columns
NULL

#' Profiles for 100 randomly selected regions from overlapping CTCF peaks in
#' 10a cell ChIP-seq
#' Results from fetching bigwigs with CTCF_in_10a_overlaps_gr.
#'
#' A tidy data.table at window size 50 bp within 350 bp of peak center
#'  The variables are as follows:
#' @details part of \code{\link{CTCF_in_10a_data}}
#' \enumerate{
#'   \item seqnames. chromosome for GRanges compatibility
#'   \item start. start of interval
#'   \item end. end of interval
#'   \item width. width of interval
#'   \item strand. leftover from GRanges.
#'   \item id. unique identifier
#'   \item y. fold-enrichment over input.
#'   \item x. bp relative to center
#'   \item sample. name of originating sample
#' }
#'
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_profiles_dt
#' @format A tidy data.table of 2100 rows and 9 columns
NULL

#' 100 randomly selected regions from overlapping CTCF peaks in 10a cell
#' ChIP-seq
#'
#' MACS2 narrowPeak calls on pooled biological replicates at pval 1e-5 and
#' then 0.05 IDR filtered. IDR cutoffs determined  by comparing top 150,000
#' pvalue sorted peak in replicates.
#'
#' See  GEO series GSE98551 for details.
#' @details part of \code{\link{CTCF_in_10a_data}}
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_overlaps_gr
#' @format GenomicRanges with 3 metadata columns of membership table
NULL

#' FTP URL path for vignette data.
#'
#' FE bigWig tracks for CTCF ChIP-seq in a MCF10A progression model.
#' See GEO series GSE98551 for details.
#' @details part of \code{\link{CTCF_in_10a_data}}
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_bigWig_urls
#' @format named character vector of length 3
NULL

#' FTP URL path for vignette data.
#' from
#'
#' macs2 peak calls for CTCF ChIP-seq in a MCF10A progression model.
#' See  GEO series GSE98551 for details.
#' @details part of \code{\link{CTCF_in_10a_data}}
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_narrowPeak_urls
#' @format named character vector of length 3
NULL

#' list of GRanges that results in 100 random subset when overlapped
#'
#' @details part of \code{\link{CTCF_in_10a_data}}
#' @docType data
#' @keywords datasets
#' @name CTCF_in_10a_narrowPeak_grs
#' @format named character vector of length 3
NULL

#' 4 random peaks for paired-end data
#'
#' matches \code{system.file("extdata/Bcell_PE.mm10.bam", package = "seqsetvis")}
#'
#' @details this is included only for testing ssvFetchBamPE functions.
#' @docType data
#' @keywords datasets
#' @name Bcell_peaks
#' @format GRanges length 4
NULL
