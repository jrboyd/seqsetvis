#' chromHMM state segmentation in the MCF7 cell line
#'
#' @description Vignette data for seqsetvis was downloaded directly from GEO series \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57498}{GSE57498}.
#' This data is the state segmentation by chromHMM in the MCF7 cell line.
#' chromHMM creates a hidden markov model by integrating several ChIP-seq
#' samples, in this case: \itemize{
#' \item MCF7_H3K27ac_ChIP-Seq
#' \item MCF7_H3K27me3_ChIP-Seq
#' \item MCF7_H3K4me1_ChIP-Seq
#' \item MCF7_H3K4me3_ChIP-Seq
#' \item MCF7_RNApolIIp_ChIP-Seq
#' }
#'
#' @description Data from GEO series \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57498}{GSE57498}
#' is from the publication \href{https://www.ncbi.nlm.nih.gov/pubmed/24916973}{Taberlay PC et al. 2014}
#' @details Contains: \itemize{
#' \item \code{\link{chromHMM_demo_overlaps_gr}}
#' \item \code{\link{chromHMM_demo_bw_states_gr}}
#' \item \code{\link{chromHMM_demo_state_total_widths}}
#' \item \code{\link{chromHMM_demo_state_colors}}
#' \item \code{\link{chromHMM_demo_segmentation_url}}
#' \item \code{\link{chromHMM_demo_chain_url}}
#' }
#' @docType data
#' @keywords datasets
#' @name chromHMM_demo_data
NULL



#' URL to download hg19 MCF7 chromHMM segmentation
#'
#' @details file is gzipped bed with name, score, itemRgb and thick meta columns
#' @details part of \code{\link{chromHMM_demo_data}}
#'
#' @docType data
#' @keywords datasets
#' @name chromHMM_demo_segmentation_url
#' @format a character containing a URL
NULL

#' URL to download hg19ToHg38 liftover chain from UCSC
#'
#' @details file is gzipped .txt
#' @details part of \code{\link{chromHMM_demo_data}}
#'
#' @docType data
#' @keywords datasets
#' @name chromHMM_demo_chain_url
#' @format a character containing a URL
NULL

#' state name to total width mappings, hg38
#'
#' @details part of \code{\link{chromHMM_demo_data}}
#'
#' @docType data
#' @keywords datasets
#' @name chromHMM_demo_state_total_widths
#' @format named numeric of total widths per state
NULL

#' original state name to color mappings stored in segmentation bed
#'
#' @details part of \code{\link{chromHMM_demo_data}}
#'
#' @docType data
#' @keywords datasets
#' @name chromHMM_demo_state_colors
#' @format a named character vector mapping states to hex colors
NULL

#' overlap of MCF10A CTCF with MCF7 chromHMM states, hg38.
#'
#' @details part of \code{\link{chromHMM_demo_data}}
#' @details the result of ssvOverlapIntervalSets() on MCF10A CTCF peaks and
#'          MCF7 chromHMM states with use_first = TRUE
#' @details first (the MCF10A peaks) and no_hit columns have been removed
#' each remaining column represents MCF10A peaks overlapping with a state.
#'
#' @docType data
#' @keywords datasets
#' @name chromHMM_demo_overlaps_gr
#' @format a GRanges object of length 98
#' with 10 logical metadata columns, 1 per state.
NULL

#' MCF10A CTCF profiles at 20 windows per chromHMM state, hg38.
#'
#' @details part of \code{\link{chromHMM_demo_data}}
#' @details the result of fetchWindowedBigwig() on the MCF10A_CTCF_FE.bw
#' near 20 randomly selected windows per chromHMM state.
#'
#' @docType data
#' @keywords datasets
#' @name chromHMM_demo_bw_states_gr
#' @format a GRanges object of length 4000
#' with 5 metadata columns sufficient for use with ggplot2
NULL
