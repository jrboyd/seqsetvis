#' easyLoad_FUN takes a character vector of file paths run an arbitrary function defined in load_FUN
#'
#' @param file_paths character vector of paths to narrowPeak files.  If named,
#' those names will be used in output unless overriden by providing file_names.
#' @param load_FUN Arbitrary function that takes at least a file path as argument.  May take other arguments that should be set in call to easyLoad_FUN.
#' @param file_names character vector of names for output list.  If not NULL
#' will override any existing names for file_paths.  Default is NULL.
#' @param n_cores number of cores to use, uses mc.cores option if set or 1.
#' @param ... extra parameters passed to load_FUN
#'
#' @return a named list of results from load_FUN
#' @export
#'
#' @examples
#' bed_f = system.file("extdata/test_loading.bed",
#'     package = "seqsetvis", mustWork = TRUE)
#' easyLoad_bed(bed_f, "my_bed")
easyLoad_FUN = function(file_paths, load_FUN,
                        file_names = NULL,
                        n_cores = getOption("mc.cores", 1),
                        ...){
    if (is.factor(file_paths)) {
        file_paths = as.character(file_paths)
    }
    stopifnot(is.character(file_paths))
    stopifnot(all(file.exists(file_paths) |
                      grepl("://", file_paths)))
    #check names of file-paths
    if (is.null(names(file_paths))) {
        names(file_paths) = basename(file_paths)
    }
    #override names with file_names if not NULL
    if (!is.null(file_names)) {
        if (is.factor(file_names)) {
            file_names = as.character(file_names)
        }
        stopifnot(is.character(file_paths))
        stopifnot(length(file_paths) == length(file_names))
        names(file_paths) = file_names
    }
    if(n_cores == 1){
        grs <- lapply(
            file_paths, load_FUN, ...)
    }else{
        grs = ssv_mclapply(
            file_paths,
            load_FUN,
            mc.cores = n_cores, ...)
    }
    grs
}

#' easyLoad_IDRmerged loads "overlapped-peaks.txt" from IDR.
#'
#' @param file_paths character vector of paths to narrowPeak files.  If named,
#' those names will be used in output unless overriden by providing file_names.
#' @param file_names character vector of names for output list.  If not NULL
#' will override any existing names for file_paths.  Default is NULL.
#' @param n_cores number of cores to use, uses mc.cores option if set or 1.
#' @param max_idr maximum IDR value allowed
#'
#' @return named list of GRanges
#' @export
#'
#' @examples
#' idr_file = system.file("extdata/test_idr.overlapped-peaks.txt",
#'     package = "seqsetvis", mustWork = TRUE)
#' easyLoad_IDRmerged(idr_file)
#' easyLoad_IDRmerged(idr_file, max_idr = .01)
easyLoad_IDRmerged = function(file_paths,
                              file_names = NULL,
                              n_cores = getOption("mc.cores", 1),
                              max_idr = .05){
    #bindings for data.table
    chr1 = start1 = start2 = stop1 = stop2 = sig.value1 = sig.value2 = idr.local = IDR = NULL
    load_idr = function(f, max_idr = .05){
        idr_dt = suppressWarnings({data.table::fread(f)[IDR <= max_idr]})
        GenomicRanges::GRanges(idr_dt[, list(seqnames = chr1, start = pmin(start1, start2), end = pmax(stop1, stop2), sig.value1, sig.value2, idr.local, IDR)])
    }
    easyLoad_FUN(file_paths, load_idr, file_names, n_cores, max_idr = max_idr)
}

#functions that load a list of files as a list of GRanges

#' easyLoad_bed takes a character vector of file paths to bed plus files and
#' returning named list of GRanges.
#'
#' Mainly a utility function for loading MACS2 narrowPeak and broadPeak.
#'
#' @param file_paths character vector of paths to narrowPeak files.  If named,
#' those names will be used in output unless overriden by providing file_names.
#' @param file_names character vector of names for output list.  If not NULL
#' will override any existing names for file_paths.  Default is NULL.
#' @param extraCols named character vector of classes.  passed to
#' rtracklayer::import for format = "BED". default is character().
#' @param n_cores number of cores to use, uses mc.cores option if set or 1.
#' @return a named list of GRanges loaded from file_paths
#' @import rtracklayer
#' @export
#' @examples
#' bed_f = system.file("extdata/test_loading.bed",
#'     package = "seqsetvis", mustWork = TRUE)
#' easyLoad_bed(bed_f, "my_bed")
easyLoad_bed = function(file_paths,
                        file_names = NULL,
                        extraCols = character(),
                        n_cores = getOption("mc.cores", 1)) {
    load_bed = function(f){
        gr = rtracklayer::import(f, format = "BED",
                                 extraCols = extraCols)
        missed = setdiff(c("name", "score", names(extraCols)), colnames(mcols(gr)))
        for(m in missed){
            if(m == "name"){
                mcols(gr)[[m]] = rep("DEFAULT", length(gr))
            }else if(m == "score"){
                mcols(gr)[[m]] = rep(0, length(gr))
            }else{
                def = switch(extraCols[m],
                             numeric = {
                                 0L
                             },
                             character = {
                                 "DEFAULT"
                             },
                             integer = {
                                 0
                             })
                mcols(gr)[[m]] = rep(def, length(gr))
            }

        }
        gr
    }

    easyLoad_FUN(file_paths,
                 load_bed,
                 file_names = file_names,
                 n_cores = n_cores)
}

#' easyLoad_narrowPeak takes a character vector of file paths to narrowPeak
#' files from MACS2 and returns a named list of GRanges.
#' @param file_paths character vector of paths to narrowPeak files.  If named,
#' those names will be used in output unless overriden by providing file_names.
#' @param file_names character vector of names for output list.  If not NULL
#' will override any existing names for file_paths.  Default is NULL.
#' @param n_cores number of cores to use, uses mc.cores option if set or 1.
#' @return a named list of GRanges loaded from file_paths
#' @import rtracklayer
#' @export
#' @examples
#' np_f = system.file("extdata/test_loading.narrowPeak",
#'     package = "seqsetvis", mustWork = TRUE)
#' easyLoad_narrowPeak(np_f, "my_narrowPeak")
easyLoad_narrowPeak = function(file_paths,
                               file_names = NULL,
                               n_cores = getOption("mc.cores", 1)) {
    #from: https://charlesjb.github.io/How_to_import_narrowPeak/
    extraCols_narrowPeak <-
        c(
            signalValue = "numeric",
            pValue = "numeric",
            qValue = "numeric",
            relSummit = "integer"
        )
    easyLoad_bed(file_paths = file_paths,
                 file_names = file_names,
                 extraCols = extraCols_narrowPeak,
                 n_cores = n_cores)
}

#' easyLoad_broadPeak takes a character vector of file paths to narrowPeak
#' files from MACS2 and returns a named list of GRanges.
#' @param file_paths character vector of paths to narrowPeak files.  If named,
#' those names will be used in output unless overriden by providing file_names.
#' @param file_names character vector of names for output list.  If not NULL
#' will override any existing names for file_paths.  Default is NULL.
#' @param n_cores number of cores to use, uses mc.cores option if set or 1.
#' @return a named list of GRanges loaded from file_paths
#' @import rtracklayer
#' @export
#' @examples
#' bp_f = system.file("extdata/test_loading.broadPeak",
#'     package = "seqsetvis", mustWork = TRUE)
#' easyLoad_broadPeak(bp_f, "my_broadPeak")
easyLoad_broadPeak = function(file_paths,
                              file_names = NULL,
                              n_cores = getOption("mc.cores", 1)) {
    #from: https://charlesjb.github.io/How_to_import_narrowPeak/
    extraCols_broadPeak <-
        c(signalValue = "numeric",
          pValue = "numeric",
          qValue = "numeric")
    easyLoad_bed(file_paths = file_paths,
                 file_names = file_names,
                 extraCols = extraCols_broadPeak,
                 n_cores = n_cores)
}

#' easyLoad_seacr takes a character vector of file paths to seacr output bed files and returns a named list of GRanges.
#' @param file_paths character vector of paths to seacr bed files.  If named,
#' those names will be used in output unless overriden by providing file_names.
#' @param file_names character vector of names for output list.  If not NULL
#' will override any existing names for file_paths.  Default is NULL.
#' @param n_cores number of cores to use, uses mc.cores option if set or 1.
#' @return a named list of GRanges loaded from file_paths
#' @import rtracklayer
#' @export
#' @examples
#' bed_f = system.file("extdata/test_loading.seacr.bed",
#'     package = "seqsetvis", mustWork = TRUE)
#' easyLoad_seacr(bed_f, "my_seacr")
easyLoad_seacr = function(file_paths,
                          file_names = NULL,
                          n_cores = getOption("mc.cores", 1)){
    extraCols_seacr = c(
        total_signal = "numeric",
        max_signal = "numeric",
        max_signal_region = "character"
    )
    easyLoad_bed(file_paths = file_paths,
                 file_names = file_names,
                 extraCols = extraCols_seacr,
                 n_cores = n_cores)
}
