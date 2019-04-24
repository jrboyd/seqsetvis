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
#' @param flip_strand boolean. if TRUE strands that reads align to are swapped.
#' This is typically only necessary if there was a mismatch between library
#' chemistry and aligner settings. Default is FALSE.
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
                                flip_strand = FALSE,
                                ...){

    x = y = N = correlation = id = NULL #reserve bindings for data.table
    if(n_regions < length(qgr)){
        qgr = sample(qgr, n_regions)
    }
    cc_res = crossCorrByRle(
        bam_f,
        qgr,
        max_dupes = 1,
        fragment_sizes = test_fragLen,
        flip_strand = flip_strand,
        ...
    )
    fragLen = cc_res[, list(shift = shift[which.max(correlation)]),
                     by = list(id)][, mean(shift)]
    fragLen = round(fragLen)
    if(!include_plot_in_output){
        return(fragLen)
    }else{
        p_res = cc_res[, list(correlation = mean(correlation)), by = list(shift)]
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
        return(list(fragLen = fragLen, plot = p))
    }
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
#' @param flip_strand boolean. if TRUE strands that reads align to are swapped.
#' This is typically only necessary if there was a mismatch between library
#' chemistry and aligner settings. Default is FALSE.
#' @param ... arguments passed to ScanBamParam
#' @return named list of results
#' @export
#' @import GenomicRanges GenomicAlignments pbapply S4Vectors
#' @importFrom Rsamtools ScanBamParam
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
                          flip_strand = FALSE,
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

    Param <- Rsamtools::ScanBamParam(
        which=query_gr,
        what=c("flag","mapq"), ...)
    if(!file.exists(paste0(bam_file, ".bai"))){
        warning("creating index for ", bam_file)
        Rsamtools::indexBam(bam_file)
    }
    temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
    if(flip_strand){
        strand(temp) = ifelse(as.character(strand(temp)) == "+", "-", "+")
    }
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
