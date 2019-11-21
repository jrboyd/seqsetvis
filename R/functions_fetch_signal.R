# functions useful for fetching signal data regardless of source:
# bam, bigwig, etc.

#' signal loading framework
#'
#' Does nothing unless load_signal is overridden to carry out reading
#' data from file_paths (likely via the appropriate ssvFetch* function,
#' ie. \code{\link{ssvFetchBigwig}} or \code{\link{ssvFetchBam}}
#'
#' @param file_paths character vector of file_paths to load from. Alternatively,
#' file_paths can be a data.frame or data.table whose first column is a
#' character vector of paths and additial columns will be used as metadata.
#' @param qgr GRanges of intervals to return from each file
#' @param unique_names unique file ids for each file in file_paths.  Default
#' is names of file_paths vector
#' @param names_variable character, variable name for column containing
#' unique_names entries.  Default is "sample"
#' @param win_size numeric/integer window size resolution to load signal at.
#' Default is 50.
#' @param win_method character.  one of c("sample", "summary").  Determines
#' if \code{\link{viewGRangesWinSample_dt}} or
#' \code{\link{viewGRangesWinSummary_dt}} is used to represent each region in
#' qgr.
#' @param return_data.table logical. If TRUE data.table is returned instead of
#' GRanges, the default.
#' @param load_signal function taking f, nam, and qgr arguments.  f is from
#' file_paths, nam is from unique_names, and qgr is qgr. See details.
#' @param n_cores integer number of cores to use.
#' Uses mc.cores option if not supplied.
#' @details load_signal is passed f, nam, and qgr and is executed in the
#' environment where load_signal is defined. See
#' \code{\link{ssvFetchBigwig}} and \code{\link{ssvFetchBam}}
#'  for examples.
#' @return A GRanges with values read from file_paths at intervals of win_size.
#' Originating file is coded by unique_names and assigned to column of name
#' names_variable.  Output is data.table is return_data.table is TRUE.
#' @export
#' @import parallel
#' @examples
#' library(GenomicRanges)
#' bam_f = system.file("extdata/test.bam",
#'     package = "seqsetvis", mustWork = TRUE)
#' bam_files = c("a" = bam_f, "b" = bam_f)
#' qgr = CTCF_in_10a_overlaps_gr[1:2]
#' qgr = resize(qgr, 500, "center")
#' load_bam = function(f, nam, qgr) {
#'     message("loading ", f, " ...")
#'     dt = seqsetvis:::ssvFetchBam.single(bam_f = f,
#'                       qgr = qgr,
#'                       win_size = 50,
#'                       fragLen = NULL,
#'                       target_strand = "*",
#'                       return_data.table = TRUE)
#'     dt[["sample"]] = nam
#'     message("finished loading ", nam, ".")
#'     dt
#' }
#' ssvFetchSignal(bam_files, qgr, load_signal = load_bam)
ssvFetchSignal = function(file_paths,
                          qgr,
                          unique_names = NULL,
                          names_variable = "sample",
                          win_size = 50,
                          win_method = c("sample", "summary")[1],
                          return_data.table = FALSE,
                          load_signal = function(f, nam, qgr) {
                              warning("nothing happened, ",
                                      "supply a function to",
                                      "load_signal parameter.")
                          },
                          n_cores = getOption("mc.cores", 1),
                          force_skip_centerFix = FALSE) {
    if (is.data.frame(file_paths) || is.data.table(file_paths)) {
        if (ncol(file_paths) == 1) {
            file_attribs = data.frame(matrix(
                0, nrow = nrow(file_paths), ncol = 0
            ))
        } else{
            file_attribs = file_paths[,-1, drop = FALSE]
        }

        file_paths = file_paths[[1]]

    } else{
        #file_paths is assumed to be a character vector
        file_attribs = data.frame(data.frame(matrix(
            0, nrow = length(file_paths), ncol = 0
        )))
        file_attribs[[names_variable]] = unique_names
    }
    if (is.list(file_paths)) {
        file_paths = unlist(file_paths)
    }
    if (is.factor(file_paths))
        file_paths = as.character(file_paths)
    if (is.null(unique_names)) {
        if(!is.null(names(file_paths))){
            unique_names = names(file_paths)
        }else if(!is.null(file_attribs[[names_variable]])){
            unique_names = file_attribs[[names_variable]]
        }else{
            unique_names = file_paths
        }
    }
    if (is.factor(unique_names))
        unique_names = as.character(unique_names)
    if (is.null(file_attribs[[names_variable]])) {
        file_attribs[[names_variable]] = unique_names
    }
    if (is.null(names(file_paths))) {
        names(file_paths) = unique_names
    }
    stopifnot(is.character(file_paths))
    stopifnot(is(qgr, "GRanges"))
    stopifnot(is.character(unique_names))
    stopifnot(is.character(names_variable))
    stopifnot(is.numeric(win_size))
    if (any(duplicated(unique_names))) {
        stop("some unique_names are duplicated:\n",
             paste(collapse = "\n",
                   unique(unique_names[duplicated(unique_names)])))
    }
    stopifnot(file.exists(file_paths))
    #if (win_method == "sample") {
    qgr = prepare_fetch_GRanges(qgr = qgr,
                                win_size = win_size,
                                target_size = NULL,
                                skip_centerFix = win_method != "sample" | force_skip_centerFix)
    #}

    nam_load_signal = function(nam) {
        f = file_paths[nam]
        load_signal(f, nam, qgr)
    }
    bw_list = parallel::mclapply(file_attribs[[names_variable]],
                                 nam_load_signal, mc.cores = n_cores)
    for (i in seq_along(bw_list)) {
        for (attrib in colnames(file_attribs)) {
            if(nrow(bw_list[[i]]) == 0){
                bw_list[[i]][[attrib]] = file_attribs[[attrib]][0]
            }else{
                bw_list[[i]][[attrib]] = file_attribs[[attrib]][i]
            }

        }

    }
    out = data.table::rbindlist(bw_list)
    if (!return_data.table) {
        out = GRanges(out)
    }
    return(out)
}

.check_qgr = function(qgr) {
    if (!is.null(qgr$name) && is.null(names(qgr))) {
        names(qgr) = qgr$name
    }
    if (is.null(qgr$id)) {
        #id not set
        if (is.null(names(qgr))) {
            #id and names not set - make names
            names(qgr) = paste0("region_", seq_along(qgr))
        }
        qgr$id = names(qgr)#names is used for id
    } else if (is.null(names(qgr))) {
        #id is used for names if id set but not names
        names(qgr) = qgr$id
    }
    qgr
}

#' get a windowed sampling of score_gr
#'
#' This method is appropriate when all GRanges in qgr are identical width and
#' when it is practical to use a window_size smaller than features in genomic
#' signal.  For instance, when retrieving signal around peaks or promoters this
#' method maintains a fixed genomic scale across regions.  This allows meaingful
#' comparison of peak widths can be made.
#'
#' Summarizes score_gr by grabbing value of "score" every window_size bp.
#' Columns in output data.table are: standard GRanges columns: seqnames, start,
#' end, width, strand id - matched to names(score_gr). if names(score_gr) is
#' missing, added as 1:length(score_gr) y - value of score from score_gr x -
#' relative bp position
#'
#' @param score_gr GRanges with a "score" metadata column.
#' @param qgr regions to view by window.
#' @param window_size qgr will be represented by value from score_gr every
#'   window_size bp.
#' @param attrib_var character name of attribute to pull data from. Default is
#'   "score", compatible with with bigWigs or bam coverage.
#' @param fill_value numeric or character value to use where queried regions are
#'   empty.  Default is 0 and appropriate for both calculated coverage and
#'   bedgraph/bigwig like files.  Will automatically switch to "MISSING" if data
#'   is guessed to be qualitative.
#' @param anchor character. controls how x value is derived from position for
#'   each region in qgr.  0 may be the left side or center.  If not unstranded,
#'   x coordinates are flipped for (-) strand. One of c("center",
#'   "center_unstranded", "left", "left_unstranded"). Default is "center".
#' @return data.table that is GRanges compatible
#' @export
#' @examples
#' bam_file = system.file("extdata/test.bam",
#'     package = "seqsetvis")
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' qgr = GenomicRanges::resize(qgr, width = 500, fix = "center")
#' bam_gr = seqsetvis:::fetchBam(bam_file, qgr)
#' bam_dt = viewGRangesWinSample_dt(bam_gr, qgr, 50)
#'
#' if(Sys.info()['sysname'] != "Windows"){
#'     bw_file = system.file("extdata/MCF10A_CTCF_FE_random100.bw",
#'         package = "seqsetvis")
#'     bw_gr = rtracklayer::import.bw(bw_file, which = qgr)
#'     bw_dt = viewGRangesWinSample_dt(bw_gr, qgr, 50)
#' }
viewGRangesWinSample_dt = function(score_gr,
                                   qgr,
                                   window_size,
                                   attrib_var = "score",
                                   fill_value = 0,
                                   anchor = c("center", "center_unstranded",
                                              "left", "left_unstranded")[1]) {
    #reserve bindings for data.table
    x = id = NULL
    stopifnot(is(score_gr, "GRanges"))
    stopifnot(!is.null(mcols(score_gr)[[attrib_var]]))
    stopifnot(is(qgr, "GRanges"))
    stopifnot(is.numeric(window_size))
    stopifnot(window_size >= 1)
    stopifnot(window_size %% 1 == 0)
    stopifnot(anchor %in% c("center", "center_unstranded",
                            "left", "left_unstranded"))
    qgr = .check_qgr(qgr)
    windows = slidingWindows(qgr, width = window_size, step = window_size)

    # names(windows) = qgr$id
    # windows's handling of names seems to have changed and now every nest
    # GRanges has parent's name
    names(windows) = NULL
    windows = unlist(windows)
    windows$id = names(windows)
    windows = resize(windows, width = 1, fix = "center")
    olaps = suppressWarnings(data.table::as.data.table(
        findOverlaps(
            query = windows,
            subject = score_gr,
            ignore.strand = TRUE
        )
    ))
    # patch up missing/out of bound data with 0
    missing_idx = setdiff(seq_along(windows), olaps$queryHits)
    if (length(missing_idx) > 0) {
        olaps = rbind(
            olaps,
            data.table::data.table(
                queryHits = missing_idx,
                subjectHits = length(score_gr) + 1
            )
        )[order(queryHits)]
        suppressWarnings({
            patch_gr = GRanges("chrPatchZero",
                               IRanges::IRanges(1, 1))
            mcols(patch_gr)[[attrib_var]] = fill_value
            score_gr = c(score_gr,
                         patch_gr)
        })
    }
    # set y and output windows = windows[olaps$queryHits]
    windows$y = fill_value
    windows[olaps$queryHits]$y = mcols(score_gr[olaps$subjectHits])[[attrib_var]]
    score_dt = data.table::as.data.table(windows)

    return(shift_anchor(score_dt, window_size, anchor))
}


#' Summarizes signal in bins.  The same number of bins per region in qgr is used
#' and widths can vary in qgr, in contrast to
#' \code{\link{viewGRangesWinSample_dt}} where width must be constant across
#' regions.
#'
#' This function is most appropriate where features are expected to vary greatly
#' in size and feature boundaries are important, ie. gene bodies, enhancers or
#' TADs.
#'
#' Columns in output data.table are: standard GRanges columns: seqnames, start,
#' end, width, strand id - matched to names(score_gr). if names(score_gr) is
#' missing, added as 1:length(score_gr) y - value of score from score_gr x -
#' relative bp position
#' @param score_gr GRanges with a "score" metadata column.
#' @param qgr regions to view by window.
#' @param n_tiles numeric >= 1, the number of tiles to use for every region in
#'   qgr.
#' @param attrib_var character name of attribute to pull data from. Default is
#'   "score", compatible with with bigWigs or bam coverage.
#' @param attrib_type one of NULL, qualitative or quantitative.  If NULL will
#'   attempt to guess by casting attrib_var attribute to character or factor.
#'   Default is NULL.
#' @param fill_value numeric or character value to use where queried regions are
#'   empty.  Default is 0 and appropriate for both calculated coverage and
#'   bedgraph/bigwig like files.  Will automatically switch to "MISSING" if data
#'   is guessed to be qualitative.
#' @param anchor character. controls how x value is derived from position for
#'   each region in qgr.  0 may be the left side or center.  If not unstranded,
#'   x coordinates are flipped for (-) strand. One of c("center",
#'   "center_unstranded", "left", "left_unstranded"). Default is "center".
#' @param summary_FUN function. used to aggregate score by tile.  must accept
#'   x=score and w=width numeric vectors as only arguments. default is
#'   weighted.mean.  limma::weighted.median is a good alternative.
#' @return data.table that is GRanges compatible
#' @export
#' @importFrom stats weighted.mean
#' @examples
#' bam_file = system.file("extdata/test.bam",
#'     package = "seqsetvis")
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' # unlike viewGRangesWinSample_dt, width is not fixed
#' # qgr = GenomicRanges::resize(qgr, width = 500, fix = "center")
#' bam_gr = seqsetvis:::fetchBam(bam_file, qgr)
#' bam_dt = viewGRangesWinSummary_dt(bam_gr, qgr, 50)
#'
#' if(Sys.info()['sysname'] != "Windows"){
#'     bw_file = system.file("extdata/MCF10A_CTCF_FE_random100.bw",
#'         package = "seqsetvis")
#'     bw_gr = rtracklayer::import.bw(bw_file, which = qgr)
#'     bw_dt = viewGRangesWinSummary_dt(bw_gr, qgr, 50)
#' }
viewGRangesWinSummary_dt = function (score_gr,
                                     qgr,
                                     n_tiles = 100,
                                     attrib_var = "score",
                                     attrib_type = NULL,
                                     fill_value = 0,
                                     anchor = c("center", "center_unstranded",
                                                "left", "left_unstranded")[1],
                                     summary_FUN = stats::weighted.mean) {
    #reserve bindings for data.table
    x = id = tile_start = tile_end = tile_id =
        tile_width = scored_width = tile_density = score_width = NULL
    stopifnot(is(score_gr, "GRanges"))
    stopifnot(!is.null(mcols(score_gr)[[attrib_var]]))
    stopifnot(is(qgr, "GRanges"))
    stopifnot(is.numeric(n_tiles))
    stopifnot(n_tiles >= 1)
    stopifnot(n_tiles %% 1 == 0)
    stopifnot(anchor %in% c("center", "center_unstranded", "left",
                            "left_unstranded"))
    qgr = .check_qgr(qgr)
    tiles = tile(qgr, n_tiles)
    # lapply(seq_len(tqgr), function(i)as.data.table(tqgr[[i]]))


    # names(tiles) = qgr$id
    # tile's handling of names seems to have changed and now every nest
    # GRanges has  parent's name
    names(tiles) = NULL
    tiles = unlist(tiles)
    tiles$id = names(tiles)
    tiles$tile_id = seq_along(tiles)


    olaps = suppressWarnings(data.table::as.data.table(
        findOverlaps(
            query = tiles,
            subject = score_gr,
            ignore.strand = TRUE
        )
    ))
    #fill in gaps with zeroes
    missing_idx = setdiff(seq_along(tiles), olaps$queryHits)
    if (length(missing_idx) > 0) {
        olaps = rbind(
            olaps,
            data.table::data.table(
                queryHits = missing_idx,
                subjectHits = length(score_gr) + 1
            )
        )[order(queryHits)]
        suppressWarnings({
            patch_gr = GRanges("chrPatchZero",
                               IRanges::IRanges(1, 1))
            mcols(patch_gr)[[attrib_var]] = NA
            score_gr = c(score_gr,
                         patch_gr
            )
        })
    }
    mcols(score_gr) = mcols(score_gr)[attrib_var]
    cov_dt = cbind(as.data.table(score_gr[olaps$subjectHits])[, -c(1, 4:5)],
                   as.data.table(tiles[olaps$queryHits])[, -c(1, 4:5)])
    colnames(cov_dt)[4:5] = c("tile_start", "tile_end")

    cov_dt[start == 1 &
               end == 1, c("start", "end") := list(tile_start, tile_end)]
    #trim score regions to fit in tiles so score weighting is accurate

    cov_dt[start < tile_start, start := tile_start]
    cov_dt[end > tile_end, end := tile_end]
    cov_dt[, score_width := end - start + 1]
    cov_dt[, tile_width := tile_end - tile_start + 1]

    if(all(is.na(cov_dt[[attrib_var]]))){
        cov_dt[[attrib_var]] = fill_value
    }
    set(cov_dt, i = which(is.na(cov_dt[[attrib_var]])), j = attrib_var, value = fill_value)

    if(is.null(attrib_type)){
        if(any(class(cov_dt[[attrib_var]]) %in% c("character", "factor"))){
            attrib_type = "qualitative"
            if(fill_value == 0){
                fill_value = "MISSING"
            }
        }else{
            attrib_type = "quantitative"
        }
    }

    # check for incompletely retrieved regions (zeroes omitted for instance)
    check_dt = cov_dt[, list(scored_width = sum(score_width)), by = list(tile_id, id, tile_width)]
    # check_dt$tile_width = width(tiles)
    #add a dummy interval of score zero to correct width
    repair_dt = check_dt[tile_width != scored_width, list(
        start = -1,
        end = -1,
        score = fill_value,
        tile_start = -1,
        tile_end = -1,
        id = id,
        tile_id = tile_id,
        score_width = tile_width - scored_width,
        tile_width
    )]
    setnames(repair_dt, "score", attrib_var)
    cov_dt = rbind(cov_dt, repair_dt)




    if(attrib_type == "qualitative"){
        density_dt = cov_dt[, list(tile_density = summary_FUN(score_width / tile_width, rep(1, length(width)))),
                            by = c("tile_id", "id", attrib_var)]
    }else if(attrib_type == "quantitative"){
        density_dt = cov_dt[, list(tile_density = summary_FUN(get(attrib_var), score_width)),
                            by = list(tile_id, id)]
    }else{
        stop("attrib_type must be one of qualitative or quantitative")
    }

    tiles = as.data.table(tiles)
    tiles[, x := (seq_len(n_tiles) - .5) / n_tiles, by = id]

    # density_dt[, x := (seq_len(n_tiles) - .5) / n_tiles, by = id]
    if (!all(tiles$id %in% density_dt$id))
        stop("something bad happened when merging density ",
             "data.table back to tiles GRanges.")

    if(attrib_type == "qualitative"){
        score_dt = merge(tiles, density_dt[, list(y = tile_density, id, tile_id, ATTRIB_VAR = get(attrib_var))], by = c("tile_id", "id"))
        setnames(score_dt, "ATTRIB_VAR", attrib_var)
        #score_dt[is.na(quality), quality := fill_value]
    }else{
        score_dt = merge(tiles, density_dt[, list(y = tile_density, id, tile_id)], by = c("tile_id", "id"))
        #score_dt[is.na(y), y := fill_value]
    }
    score_dt$tile_id = NULL
    #slightly different than summary,
    #x is already set and regions are already contiguous.
    #just need to flip x or center as needed.
    switch(
        anchor,
        center = {
            center_val = mean((seq_len(n_tiles) - .5) / n_tiles)
            score_dt[, x := x - center_val, by = id]
            score_dt[strand == "-", x := (-1 * x)]
        },
        center_unstranded = {
            center_val = mean((seq_len(n_tiles) - .5) / n_tiles)
            score_dt[, x := x - center_val, by = id]
        },
        left = {
            score_dt[strand == "-", x := (1 - 1 * x), by = id]
        },
        left_unstranded = {
            #do nothing
        }
    )

    #ensure x is rounded and easily flippable
    score_dt$x = round(score_dt$x, 9)

    score_dt
}

#' orients the relative position of x's zero value and
#' extends ranges to be contiguous
#'
#' @param score_dt data.table, GRanges() sufficient
#' @param window_size numeric, window size used to generate socre_dt
#' @param anchor character, one of c("center", "center_unstranded",
#' "left", "left_unstranded")
#' @return score_dt with x values shifted appropriately and start and end
#' extended to make ranges contiguous
shift_anchor = function(score_dt, window_size, anchor) {
    x = id = NULL
    shift = round(window_size / 2)
    switch(
        anchor,
        center = {
            score_dt[, x := start - min(start) + shift, by = id]
            score_dt[, x := x - round(mean(x)), by = id]
            score_dt[strand == "-", x := -1 * x]
        },
        center_unstranded = {
            score_dt[, x := start - min(start) + shift, by = id]
            score_dt[, x := x - round(mean(x)), by = id]
        },
        left = {
            score_dt[, x := -1]
            score_dt[strand != "-", x := start - min(start) + shift,
                     by = id]
            #flip negative
            score_dt[strand == "-", x := -1 * (end - max(end) - shift),
                     by = id]
        },
        left_unstranded = {
            score_dt[, x := start - min(start) + shift, by = id]
        }
    )

    score_dt[, start := start - shift]
    score_dt[, end := end + window_size - shift - 1]
    return(score_dt)
}

#' prepares GRanges for windowed fetching.
#'
#' output GRanges parallels input with consistent width evenly divisible by
#' win_size.  Has warning if GRanges needed resizing, otherwise no warning
#' and input GRanges is returned unchanged.
#'
#' @param qgr GRanges to prepare
#' @param win_size numeric window size for fetch
#' @param min_quantile numeric [0,1], lowest possible quantile value.  Only
#' relevant if target_size is not specified.
#' @param target_size numeric final width of qgr if known. Default of NULL
#' leads to quantile based determination of target_size.
#' @param skip_centerFix boolean, if FALSE (default) all regions will be resized
#' GenomicRanges::resize(x, w, fix = "center") to a uniform size based on
#' min_quantile to a width divisible by win_size.
#' @return GRanges, either identical to qgr or with suitable consistent width
#' applied.
#' @export
#' @examples
#' qgr = prepare_fetch_GRanges(CTCF_in_10a_overlaps_gr, win_size = 50)
#' #no warning if qgr is already valid for windowed fetching
#' prepare_fetch_GRanges(qgr, win_size = 50)
prepare_fetch_GRanges = function(qgr,
                                 win_size,
                                 min_quantile = .75,
                                 target_size = NULL,
                                 skip_centerFix = FALSE) {
    if (!skip_centerFix &&
        (length(unique(width(qgr))) > 1 || width(qgr)[1] %% win_size != 0)) {
        if (is.null(target_size)) {
            target_size = quantileGRangesWidth(
                qgr = qgr,
                min_quantile = min_quantile,
                win_size = win_size
            )
        }
        if (target_size %% win_size != 0) {
            stop(
                "target_size: ",
                target_size,
                " not evenly divisible by win_size: ",
                win_size
            )
        }

        qgr = centerFixedSizeGRanges(qgr, fixed_size = target_size)
        message(
            "widths of qgr were not ",
            "identical and evenly divisible by win_size.",
            "\nA fixed width of ",
            target_size,
            " was applied based on the data provided."
        )
    }
    if (any(start(qgr) < 1)) {
        warning("Some out of bounds GRanges had to be shifted back to start >= 1.")
        fix_shift = 1 - start(qgr)
        fix_shift = vapply(fix_shift, function(x)
            max(x, 0), FUN.VALUE = 1)
        qgr = IRanges::shift(qgr, fix_shift)
    }
    return(qgr)
}

#' Quantile width determination strategy
#'
#' Returns the lowest multiple of win_size greater than
#' min_quantile quantile of width(qgr)
#'
#' @param qgr GRanges to calculate quantile width for
#' @param min_quantile numeric [0,1] the minimum quantile of width in qgr
#' @param win_size numeric/integer >=1, returned value will be a multiple of
#' this
#' @return numeric that is >= min_quantile and evenly divisible by win_size
#' @export
#' @examples
#' gr = CTCF_in_10a_overlaps_gr
#' quantileGRangesWidth(gr)
#' quantileGRangesWidth(gr, min_quantile = .5, win_size = 100)
quantileGRangesWidth = function(qgr,
                                min_quantile = .75,
                                win_size = 1) {
    stopifnot(is(qgr, "GRanges"))
    stopifnot(is.numeric(min_quantile), is.numeric(win_size))
    stopifnot(min_quantile >= 0 && min_quantile <= 1)
    stopifnot(length(min_quantile) == 1 && length(win_size) == 1)
    stopifnot(win_size %% 1 == 0)
    stopifnot(win_size >= 1)
    qwidth = quantile(width(qgr), min_quantile)
    fwidth = ceiling(qwidth / win_size) * win_size
    return(fwidth)
}
