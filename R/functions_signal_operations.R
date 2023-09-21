#' Transforms set of GRanges to all have the same size.
#'
#' \code{centerFixedSizeGRanges} First calculates the central coordinate of each
#' GRange in \code{grs} and extends in both direction by half of
#' \code{fixed_size}
#' @export
#' @param grs Set of GRanges with incosistent and/or incorrect size
#' @param fixed_size The final width of each GRange returned.
#' @return Set of GRanges after resizing all input GRanges, either shortened
#' or lengthened as required to match \code{fixed_size}
#' @import GenomicRanges
#' @examples
#' library(GenomicRanges)
#' grs = GRanges("chr1", IRanges(1:10+100, 1:10*3+100))
#' centered_grs = centerFixedSizeGRanges(grs, 10)
#' width(centered_grs)
centerFixedSizeGRanges = function(grs, fixed_size = 2000) {
    stopifnot(is(grs, "GRanges"))
    stopifnot(is(fixed_size, "numeric"))
    stopifnot(fixed_size > 0)
    m = floor(start(grs) + width(grs)/2)
    ext = floor(fixed_size/2)
    start(grs) = m - ext
    end(grs) = m + fixed_size - ext - 1
    # resize isn't ideal - repeated applications accumulate rounding shifts
    # grs = GenomicRanges::resize(x = grs, width = fixed_size, fix = "center")
    return(grs)
}


#' applies a spline smoothing to a tidy data.table containing x and y values.
#'
#' \code{applySpline} Is intended for two-dimensional tidy data.tables, as
#' retured by \code{ssvFetchBigwig}
#' @export
#' @param dt a tidy data.table containing two-dimensional data
#' @param n the number of interpolation points to use per input point, see
#' \code{?spline}.  n must be > 1.
#' @param x_ the variable name of the x-values
#' @param y_ the variable name of the y-values
#' @param by_ optionally, any variables that provide grouping to the data.
#' default is none. see details.
#' @param splineFun a function that accepts x, y, and n as arguments and
#' returns a list of length 2 with named elements x and y.
#' \code{stats::spline} by default.
#' see \code{stats::spline} for details.
#'
#' @return a newly derived data.table that is \code{n} times longer than
#' original.
#'
#' @details by_ is quite powerful.  If \code{by_ = c('gene_id', 'sample_id')},
#' splines
#' will be calculated individually for each gene in each sample. alternatively
#' if \code{by_ = c('gene_id')}
#' @seealso \code{\link{ssvFetchBigwig}}
#' @importFrom stats spline
#' @examples
#' #data may be blockier than we'd like
#' ggplot(CTCF_in_10a_profiles_dt[, list(y = mean(y)), by = list(sample, x)]) +
#'     geom_line(aes(x = x, y = y, color = sample))
#'
#' #can be smoothed by applying a spline  (think twice about doing so,
#' #it may look prettier but may also be deceptive or misleading)
#'
#' splined_smooth = applySpline(CTCF_in_10a_profiles_dt, n = 10,
#'     y_ = 'y', by_ = c('id', 'sample'))
#' ggplot(splined_smooth[, list(y = mean(y)), by = list(sample, x)]) +
#'     geom_line(aes(x = x, y = y, color = sample))
applySpline = function(dt, n, x_ = "x", y_ = "y", by_ = c("id", "sample"),
                       splineFun = stats::spline) {
    output_GRanges = FALSE
    if(is(dt, "GRanges")){
        dt = as.data.table(dt)
        output_GRanges = TRUE
    }
    stopifnot(data.table::is.data.table(dt))
    stopifnot(is.character(x_), is.character(y_), is.character(by_))
    stopifnot(is.function(splineFun))
    if (!any(x_ == colnames(dt))) {
        stop("applySpline : x_ (", x_,
             ") not found in colnames of input data.table")
    }
    if (!any(y_ == colnames(dt))) {
        stop("applySpline : y_ (", y_,
             ") not found in colnames of input data.table")
    }
    if (by_[1] != "" | length(by_) > 1)
        if (!all(by_ %in% colnames(dt))) {
            stop("applySpline : by_ (", by_,
                 ") not found in colnames of input data.table")
        }
    dt = dt[order(get(x_))]
    if(by_[1] != ""){
        for(.by_ in by_){
            dt = dt[order(get(.by_))]
        }
    }

    stopifnot(n > 1)
    dupe_x_within_by = suppressWarnings(
        any(dt[, any(duplicated(get(x_))), by = by_]$V1))
    if (dupe_x_within_by)
        warning("applySpline : Duplicate values of x_ (\"", x_,
                "\") exist within groups defined with by_ (\"", by_, "\"). ",
                "This Results in splines through the means of yvalues at",
                " duplicated xs.")
    extra_cols = setdiff(colnames(dt), c(x_, y_, by_))
    # sdt = dt[, list(n = floor(.N * n)), by = by_]
    sdt = dt[, splineFun(x = get(x_), y = get(y_), n = floor(.N * n)), by = by_]
    colnames(sdt)[colnames(sdt) == "x"] = x_
    colnames(sdt)[colnames(sdt) == "y"] = y_

    #repair with columns dropped in by_ apply spline
    #each row will be duplicated n times
    if(length(extra_cols) > 0){
        if(n > 1){
            repair = dt[rep(seq_len(nrow(dt)), each = n),
                        c(extra_cols, by_[by_ != ""]), with = FALSE]
            sdt = cbind(sdt, repair)
        }else{
            # warning("")
            # repair = unique(dt[, c(extra_cols, by_, x_), with = FALSE])
            # repair = dt
            # sdt
            # merge(sdt, repair, by = by_)
            # unique(sdt[, by_, with = FALSE])
            # merge(sdt, repair, by = by_)
        }

    }

    k = colnames(dt) %in% colnames(sdt)
    sdt = sdt[, colnames(dt)[k], with = FALSE]
    if(output_GRanges){
        sdt = GRanges(sdt)
    }
    return(sdt)
}

#' centers profile of x and y.  default is to center by region but across all
#' samples.
#'
#' \code{centerAtMax} locates the coordinate x of the maximum in y and shifts x
#' such that it is zero at max y.
#' @export
#' @param dt data.table
#' @param x_ the variable name of the x-values. default is 'x'
#' @param y_ the variable name of the y-values default is 'y'
#' @param by_ optionally, any variables that provide grouping to the data.
#' default is none.  see details.
#' @param view_size the size in \code{x_} to consider for finding the max
#' of \code{y_}.
#' if length(view_size) == 1, range will be c(-view_size, view_size).
#' if length(view_size) > 1, range will be range(view_size).
#' default value of NULL uses complete range of x.
#' @param replace_x logical, default TRUE.
#' if TRUE x_ will be replaced with position relative to summit.
#' if FALSE x_ will be preserved and x_summitPosition added.
#' @param trim_to_valid valid \code{x_} values are those with a set \code{y_}
#' value in all \code{by_} combinations
#' @param check_by_dupes default assumption is that there should be on set of
#' x_ for a by_ instance.
#' if this is not the case and you want to disable warnings about set this
#' to FALSE.
#' @param x_precision numerical precision of x, default is 3.
#' @details character.  by_ controls at the level of the data centering is
#' applied.  If by_ is "" or NULL, a single max position will be determined
#' for the entire dataset.  If by is "id" (the default) then each region will be
#' centered individually across all samples.
#' @return
#' data.table with x (or xnew if replace_x is FALSE) shifted such that
#' x = 0 matches the maximum y-value define by by_ grouping
#' @examples
#' centerAtMax(CTCF_in_10a_profiles_gr, y_ = 'y', by_ = 'id',
#'   check_by_dupes = FALSE)
#' #it's a bit clearer what's happening with trimming disabled
#' #but results are less useful for heatmaps etc.
#' centerAtMax(CTCF_in_10a_profiles_gr, y_ = 'y', by_ = 'id',
#'   check_by_dupes = FALSE, trim_to_valid = FALSE)
#' #specify view_size to limit range of x values considered, prevents
#' #excessive data trimming.
#' centerAtMax(CTCF_in_10a_profiles_gr, y_ = 'y', view_size = 100, by_ = 'id',
#' check_by_dupes = FALSE)
centerAtMax = function(dt,
                       x_ = "x",
                       y_ = "y",
                       by_ = "id",
                       view_size = NULL,
                       trim_to_valid = TRUE,
                       check_by_dupes = TRUE,
                       x_precision = 3,
                       replace_x = TRUE) {

    ymax = xsummit = xnew = N = NULL  #reserve data.table variables
    output_GRanges = FALSE
    if(is(dt, "GRanges")){
        dt = data.table::as.data.table(dt)
        output_GRanges = TRUE
    }
    if (!data.table::is.data.table(dt)) {
        stop("dt must be of type data.table, was ", class(dt))
    }
    stopifnot(is.character(x_),
              is.character(y_),
              is.character(by_) || is.null(by_))
    stopifnot(is.numeric(view_size) || is.null(view_size))
    stopifnot(is.logical(trim_to_valid),
              is.logical(check_by_dupes),
              is.logical(replace_x))

    if (!any(x_ == colnames(dt))) {
        stop("centerAtMax : x_ (", x_,
             ") not found in colnames of input data.table")
    }
    if (!any(y_ == colnames(dt))) {
        stop("centerAtMax : y_ (", y_,
             ") not found in colnames of input data.table")
    }
    if(length(unique(dt[[x_]])) != length(unique(round(dt[[x_]], x_precision)))){
        stop(
            "centerAtMax : x_precision,",
            x_precision,
            " , is too low, how many decimal places is x defined to?"
        )
    }
    # check_by_dupes = FALSE
    if (is.null(by_)) {
        by_ = ""
        check_by_dupes = TRUE
    }
    if (all(by_ != ""))
        if (!any(by_ %in% colnames(dt))) {
            stop("centerAtMax : by_ (", by_,
                 ") not found in colnames of input data.table")
        }
    if (check_by_dupes) {
        dupe_x_within_by = suppressWarnings(any(dt[, any(duplicated(get(x_))),
                                                   by = by_]$V1))
        if (dupe_x_within_by)
            message("centerAtMax : duplicate values of x_ (", x_,
                    ") exist within groups defined with by_ (", by_,
                    ").\n    If this is the desired functionality, set",
                    "check_by_dupes <- FALSE to hide future messages. ",
                    "If no by_ grouping is intended set by_ <- \"\" as",
                    "well.")
    }
    dt = data.table::copy(dt)
    if (is.null(view_size)) {
        view_size = range(dt[[x_]])
    } else if (length(view_size) == 1) {
        view_size = c(-view_size, view_size)
    }
    view_size = range(view_size)
    closestToZero = function(x) {
        x[order(abs(x))][1]
    }
    dt[, `:=`(ymax, max(get(y_)[get(x_) <= max(view_size) &
                                    get(x_) >= min(view_size)])), by = by_]
    dt[, `:=`(xsummit, closestToZero(get(x_)[get(y_) == ymax])), by = by_]
    dt[, `:=`(xnew, get(x_) - xsummit)]
    dt$xnew = round(dt$xnew, x_precision)
    dt[, `:=`(ymax, NULL)]
    dt[, `:=`(xsummit, NULL)]
    if (trim_to_valid) {
        # valid values of x are those with values in all by_ defined grouping
        xcounts = dt[, .N, xnew]
        xcounts = xcounts[N == max(N)]
        dt = dt[xnew %in% xcounts$xnew]
    }
    if (replace_x) {
        data.table::set(dt, j = x_, value = dt$xnew)
        dt$xnew = NULL
    } else {
        colnames(dt)[colnames(dt) == "xnew"] = paste0(x_, "_summitPosition")
    }
    if(output_GRanges){
        dt = GRanges(dt)
    }

    return(dt)
}

#' Centers query GRanges at maximum signal in prof_dt.
#'
#' @param prof_dt a GRanges or data.table as returned by ssvFetch*.
#' @param qgr the GRanges used to query ssvFetch* as the qgr argument.
#' @param x_ positional variable.  Should almost always be the default, "x".
#' @param y_ the signal value variable.  Likely the default value of "y" but
#'   could be "y_norm" if append_ynorm was applied to data.
#' @param by_ region identifier variable. Should almost always be the default,
#'   "id".
#' @param width Desired width of final regions.  Default is 1.
#'
#' @return a GRanges with same mcols as qgr that has been centered based on
#'   signal in prof_dt and with regions of specified width.
#' @export
#'
#' @examples
#' centerGRangesAtMax(CTCF_in_10a_profiles_dt, CTCF_in_10a_overlaps_gr)
#' centerGRangesAtMax(CTCF_in_10a_profiles_gr, CTCF_in_10a_overlaps_gr)
#'
centerGRangesAtMax = function(prof_dt, qgr, x_ = "x", y_ = "y", by_ = "id", width = 1){
    if(length(by_) > 1) stop("only by_ of length 1 supported.")
    if(is(prof_dt, "GRanges")){
        prof_dt = data.table::as.data.table(prof_dt)
    }
    cent_dt = centerAtMax(prof_dt, y_ = y_, by_ = by_, check_by_dupes = FALSE)
    cent_dt = cent_dt[, .SD[which(get(x_) == min(abs(get(x_))))[1],], c(by_)]
    cent_gr = GenomicRanges::GRanges(cent_dt[, list(seqnames, start = (start + end)/2, end = (start + end)/2)])
    names(cent_gr) = cent_dt[[by_]]

    qgr = prepare_fetch_GRanges_names(qgr)
    GenomicRanges::strand(cent_gr) = GenomicRanges::strand(qgr[names(cent_gr)])
    GenomicRanges::mcols(cent_gr) = GenomicRanges::mcols(qgr[names(cent_gr)])
    cent_gr = cent_gr[names(qgr)]
    GenomicRanges::resize(cent_gr, width, fix = "center")
}

#' findMaxPos
#'
#' @param prof_dt a GRanges or data.table as returned by ssvFetch*.
#' @param qgr the GRanges used to query ssvFetch* as the qgr argument.
#' @param x_ positional variable.  Should almost always be the default, "x".
#' @param y_ the signal value variable.  Likely the default value of "y" but
#'   could be "y_norm" if append_ynorm was applied to data.
#' @param by_ region identifier variable. Should almost always be the default,
#'   "id".
#' @param width Desired width of final regions.  Default is 1.
#'
#' @return data.table of relative x position from center per id
#' @export
#'
#' @examples
#' findMaxPos(CTCF_in_10a_profiles_dt, CTCF_in_10a_overlaps_gr)
#' findMaxPos(CTCF_in_10a_profiles_gr, CTCF_in_10a_overlaps_gr)
findMaxPos = function(prof_dt, qgr, x_ = "x", y_ = "y", by_ = "id", width = 1){
    #binding for data.table
    x = y = id = NULL
    if(length(by_) > 1) stop("only by_ of length 1 supported.")
    if(is(prof_dt, "GRanges")){
        prof_dt = data.table::as.data.table(prof_dt)
    }
    prof_dt[, list(x = x[which.max(y)[1]]), list(id)]
}
