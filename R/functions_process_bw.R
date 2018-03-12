#' Transforms set of GRanges to all have the same size.
#'
#' \code{centerFixedSizeGranges} First calculates the central coordinate of each
#' GRange in \code{grs} and extends in both direction by half of \code{fixed_size}
#' @export
#' @param grs Set of GRanges with incosistent and/or incorrect size
#' @param fixed_size The final width of each GRange returned.
#' @return Set of GRanges after resizing all input GRanges, either shortened
#' or lengthened as required to match \code{fixed_size}
centerFixedSizeGRanges = function(grs, fixed_size = 2000) {
    m = floor(start(grs) + width(grs)/2)
    ext = floor(fixed_size/2)
    start(grs) = m - ext
    end(grs) = m + fixed_size - ext - 1
    return(grs)
}

#' Fetch values from a bigwig appropriate for heatmaps etc.
#'
#' \code{fetchWindowedBigwig} Gets values for each region of the query GRanges (\code{qgr}).
#' Values correspond to the center of each window of size \code{win_size}.  A tidy formatted data.table
#' object is returned suitable for plotting using ggplots.
#' @export
#' @param bw_file The character vector path to bigwig files to read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @return A tidy formatted data.table containing fetched values.
#'
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved from \code{bw_file}
#'
fetchWindowedBigwig = function(bw_file, qgr, win_size = 50) {
    queryHits = id = x = NULL
    if (!all(width(qgr)%%win_size == 0)) {
        stop(paste("all widths of qgr are not evenly divisible by win_size,", win_size))
    }
    if (!all(width(qgr) == width(qgr)[1])) {
        warning(paste("Widths vary between GRanges.  Before aggregating or co-plotting, recommend one of:
                  1) modifying input qgr such that all GRanges have same width
                  2) filtering output so that the same values of x are present for every region"))
    }
    # suppressWarnings({
    bw_gr = rtracklayer::import.bw(bw_file, which = qgr)
    # })
    windows = slidingWindows(qgr, width = win_size, step = win_size)
    if (is.null(qgr$id)) {
        if (!is.null(names(qgr))) {
            qgr$id = names(qgr)
        } else {
            qgr$id = paste0("region_", seq_along(qgr))
        }
    }
    names(windows) = qgr$id
    # print(object.size(windows), units = 'GB')
    windows = unlist(windows)
    windows$id = names(windows)
    mid_gr = function(gr) {
        start(gr) + floor((width(gr) - 1)/2)
    }
    mids = mid_gr(windows)
    start(windows) = mids
    end(windows) = mids
    olaps = suppressWarnings(data.table::as.data.table(findOverlaps(query = windows, subject = bw_gr)))
    # patch up missing/out of bound data with 0
    missing_idx = setdiff(1:length(windows), olaps$queryHits)
    if (length(missing_idx) > 0) {
        olaps = rbind(olaps, data.table::data.table(queryHits = missing_idx, subjectHits = length(bw_gr) + 1))[order(queryHits)]
        bw_gr = c(bw_gr, GRanges(seqnames(bw_gr)[length(bw_gr)], IRanges::IRanges(1, 1), score = 0))
    }
    # set FE and output windows = windows[olaps$queryHits]
    windows$FE = bw_gr[olaps$subjectHits]$score
    bw_dt = data.table::as.data.table(windows)
    bw_dt[, `:=`(x, start - min(start) + win_size/2), by = id]
    bw_dt[, `:=`(x, x - round(mean(x))), by = id]
    shift = round(win_size/2)
    bw_dt[, `:=`(start, start - shift + 1)]
    bw_dt[, `:=`(end, end + win_size - shift)]
    return(bw_dt)
}

#' Iterates a character vector (ideally named) and calls \code{fetchWindowedBigwig}
#' on each.  Appends grouping variable to each resulting data.table and uses rbindlist to
#' efficiently combine results
#'
#' \code{fetchWindowedBigwigList} iteratively calls \code{fetchWindowedBigwig}.
#' See \code{\link{fetchWindowedBigwig}} for more info.
#' @export
#' @param bw_files The character vector paths to bigwig files to read from.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#' interval should be identical and evenly divisible by \code{win_size}.
#' @param bw_names names to use in final data.table to designate source bigwig
#' @param bw_variable_name The column name where bw_names are stored.
#' Default is 'sample'
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @return A tidy formatted data.table containing fetched values.
#'
#' @details if \code{qgr} contains the range chr1:1-100 and \code{win_size} is
#' 10, values from positions chr1 5,15,25...85, and 95 will be retrieved from \code{bw_file}
fetchWindowedBigwigList = function(bw_files, qgr, bw_names = names(bw_files), bw_variable_name = "sample", win_size = 50) {
    if (is.null(bw_names)) {
        bw_names = basename(bw_files)
    }
    if (any(duplicated(bw_names))) {
        stop(paste0("some bw_names are duplicated:\n", paste(collapse = "\n", unique(bw_names[duplicated(bw_names)]))))
    }
    bw_list = lapply(names(bw_files), function(nam) {
        print(paste0("loading ", nam, " ..."))
        f = bw_files[nam]
        dt = fetchWindowedBigwig(bw_file = f, win_size = win_size, qgr = qgr)
        dt[[bw_variable_name]] = nam
        print(paste0("finished loading ", nam, "."))
        dt
    })
    all_bw_dt = data.table::rbindlist(bw_list)
    return(all_bw_dt)
}


#' applies a spline smoothing to a tidy data.table containing x and y values.
#'
#' \code{applySpline} Is intended for two-dimensional tidy data.tables, as
#' retured by \code{fetchWindowedBigwig}
#' @export
#' @param dt a tidy data.table containing two-dimensional data
#' @param n the number of interpolation points to use per input point, see
#' \code{?spline}
#' @param x_ the variable name of the x-values
#' @param y_ the variable name of the y-values
#' @param by_ optionally, any variables that provide grouping to the data.
#' default is none. see details.
#' @param ... additional arguments to pass to \code{\link{spline}}
#' @return a newly derived data.table that is \code{n} times longer than
#' original.
#'
#' @details by_ is quite powerful.  If \code{by_ = c('gene_id', 'sample_id')},
#' splines
#' will be calculated individually for each gene in each sample. alternatively
#' if \code{by_ = c('gene_id')}
#' @seealso \code{\link{fetchWindowedBigwig}}
#' @examples
#' #data may be blockier than we'd like
#' ggplot(CTCF_in_10a_profiles_dt[, .(FE = mean(FE)), by = .(sample, x)]) +
#' geom_line(aes(x = x, y = FE, color = sample))
#' splined_up = applySpline(CTCF_in_10a_profiles_dt, n = 2, y_ = 'FE',
#' by_ = c('id', 'sample'))
#' #can be smoothed by applying a spline  (think twice about doing so,
#' #it may look prettier but may also be deceptive or misleading)
#' ggplot(splined_up[, .(FE = mean(FE)), by = .(sample, x)]) +
#' geom_line(aes(x = x, y = FE, color = sample))
#' splined_smooth = applySpline(CTCF_in_10a_profiles_dt, n = 10,
#' y_ = 'FE', by_ = c('id', 'sample'))
#' ggplot(splined_smooth[, .(FE = mean(FE)), by = .(sample, x)]) +
#' geom_line(aes(x = x, y = FE, color = sample))
#' #another potential use is to down sample
#' splined_down = applySpline(CTCF_in_10a_profiles_dt, n = .5,
#' y_ = 'FE', by_ = c('id', 'sample'))
#' ggplot(splined_down[, .(FE = mean(FE)), by = .(sample, x)]) +
#' geom_line(aes(x = x, y = FE, color = sample))
applySpline = function(dt, n, x_ = "x", y_ = "y", by_ = "", ...) {
    if (!data.table::is.data.table(dt)) {
        stop(paste("dt must be of type data.table, was", class(dt)))
    }
    if (!any(x_ == colnames(dt))) {
        stop(paste("applySpline : x_ (", x_, ") not found in colnames of input data.table"))
    }
    if (!any(y_ == colnames(dt))) {
        stop(paste("applySpline : y_ (", y_, ") not found in colnames of input data.table"))
    }
    if (by_[1] != "" | length(by_) > 1)
        if (!all(by_ %in% colnames(dt))) {
            stop(paste("applySpline : by_ (", by_, ") not found in colnames of input data.table"))
        }

    dupe_x_within_by = suppressWarnings(any(dt[, any(duplicated(get(x_))), by = by_]$V1))
    if (dupe_x_within_by)
        warning(paste0("applySpline : Duplicate values of x_ (", x_, ") exist within groups defined with by_ (", by_, ").\n
                                      This Results in splines through the means of yvalues at duplicated xs."))

    sdt = dt[, spline(get(x_), get(y_), n = floor(.N * n), ...), by = by_]
    colnames(sdt)[colnames(sdt) == "x"] = x_
    colnames(sdt)[colnames(sdt) == "y"] = y_
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
#' @details by_ is quite powerful.  If \code{by_ = c('gene_id', 'sample_id')},
#' splines
#' will be calculated individually for each gene in each sample. alternatively
#' if \code{by_ = c('gene_id')}
#' @examples
#' centerAtMax(CTCF_in_10a_profiles_dt, y_ = 'FE', by_ = 'id',
#'   check_by_dupes = FALSE)
#' #it's a bit clearer what's happening with trimming disabled
#' centerAtMax(CTCF_in_10a_profiles_dt, y_ = 'FE', by_ = 'id',
#'   check_by_dupes = FALSE, trim_to_valid = FALSE)
#' #specify view_size to limit range of x values considered, prevents
#' #excessive data trimming.
#' centerAtMax(CTCF_in_10a_profiles_dt, y_ = 'FE', view_size = 100, by_ = 'id',
#' check_by_dupes = FALSE)
centerAtMax = function(dt, x_ = "x", y_ = "y", by_ = NULL, view_size = NULL, trim_to_valid = TRUE, check_by_dupes = TRUE, replace_x = TRUE) {
    ymax = xsummit = xnew = N = NULL  #reserve data.table variables
    if (!data.table::is.data.table(dt)) {
        stop(paste("dt must be of type data.table, was", class(dt)))
    }
    if (!any(x_ == colnames(dt))) {
        stop(paste("centerAtMax : x_ (", x_, ") not found in colnames of input data.table"))
    }
    if (!any(y_ == colnames(dt))) {
        stop(paste("centerAtMax : y_ (", y_, ") not found in colnames of input data.table"))
    }
    # check_by_dupes = FALSE
    if (is.null(by_)) {
        by_ = ""
        check_by_dupes = TRUE
    }
    if (all(by_ != ""))
        if (!any(by_ %in% colnames(dt))) {
            stop(paste("centerAtMax : by_ (", by_, ") not found in colnames of input data.table"))
        }
    if (check_by_dupes) {
        dupe_x_within_by = suppressWarnings(any(dt[, any(duplicated(get(x_))), by = by_]$V1))
        if (dupe_x_within_by)
            message(paste0("centerAtMax : duplicate values of x_ (", x_, ") exist within groups defined with by_ (", by_, ").\n    If this is the desired functionality, set check_by_dupes <- FALSE to hide future messages. If no by_ grouping is intended set by_ <- \"\" as well."))
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
    # dt[, `:=`(ymax, max(get(y_)[get(x_) <= max(view_size) & get(x_) >= min(view_size)])), by = by_] dt[, `:=`(xsummit,
    # closestToZero(get(x_)[get(y_) == ymax])), by = by_] dt[, `:=`(xnew, get(x_) - xsummit)] dt[, `:=`(ymax, NULL)] dt[, `:=`(xsummit,
    # NULL)] dt[, data.table::`:=`(ymax, max(get(y_)[get(x_) <= max(view_size) & get(x_) >= min(view_size)])), by = by_] dt[,
    # data.table::`:=`(xsummit, closestToZero(get(x_)[get(y_) == ymax])), by = by_] dt[, data.table::`:=`(xnew, get(x_) - xsummit)]
    # dt[, data.table::`:=`(ymax, NULL)] dt[, data.table::`:=`(xsummit, NULL)]
    dt[, `:=`(ymax, max(get(y_)[get(x_) <= max(view_size) & get(x_) >= min(view_size)])), by = by_]
    dt[, `:=`(xsummit, closestToZero(get(x_)[get(y_) == ymax])), by = by_]
    dt[, `:=`(xnew, get(x_) - xsummit)]
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
    return(dt)
}


