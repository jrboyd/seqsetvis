#functions that load a list of files as a list of GRanges

#' easyLoad_narrowPeak takes a character vector of file paths to narrowPeak
#' files from MACS2 and returns a named list of GRanges.
#' @param file_paths character vector of paths to narrowPeak files.  If named,
#' those names will be used in output unless overriden by providing file_names.
#' @param file_names character vector of names for output list.  If not NULL will
#' override any existing names for file_paths.  Default is NULL.
#' @return a named list of GRanges loaded from file_paths
#' @import rtracklayer
#' @export
#' @examples
#' np_f = system.file("extdata/test_loading.narrowPeak",
#'     package = "seqsetvis", mustWork = TRUE)
#' easyLoad_narrowPeak(np_f, "my_narrowPeak")
easyLoad_narrowPeak = function(file_paths, file_names = NULL){
    #from: https://charlesjb.github.io/How_to_import_narrowPeak/
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                              qValue = "numeric", relSummit = "integer")
    easyLoad_bed(file_paths = file_paths,
                     file_names = file_names,
                     extraCols = extraCols_narrowPeak)
}

#' easyLoad_broadPeak takes a character vector of file paths to narrowPeak
#' files from MACS2 and returns a named list of GRanges.
#' @param file_paths character vector of paths to narrowPeak files.  If named,
#' those names will be used in output unless overriden by providing file_names.
#' @param file_names character vector of names for output list.  If not NULL will
#' override any existing names for file_paths.  Default is NULL.
#' @return a named list of GRanges loaded from file_paths
#' @import rtracklayer
#' @export
#' @examples
#' bp_f = system.file("extdata/test_loading.broadPeak",
#'     package = "seqsetvis", mustWork = TRUE)
#' easyLoad_broadPeak(bp_f, "my_broadPeak")
easyLoad_broadPeak = function(file_paths, file_names = NULL){
    #from: https://charlesjb.github.io/How_to_import_narrowPeak/
    extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                             qValue = "numeric")
    easyLoad_bed(file_paths = file_paths,
                     file_names = file_names,
                     extraCols = extraCols_broadPeak)
}

#' easyLoad_bed takes a character vector of file paths to bed plus files and
#' returning named list of GRanges.
#' Mainly a utility function for loading MACS2 narrowPeak and broadPeak.
#' @param file_paths character vector of paths to narrowPeak files.  If named,
#' those names will be used in output unless overriden by providing file_names.
#' @param file_names character vector of names for output list.  If not NULL will
#' override any existing names for file_paths.  Default is NULL.
#' @param extraCols named character vector of classes.  passed to
#' rtracklayer::import for format = "BED". default is character().
#' @return a named list of GRanges loaded from file_paths
#' @import rtracklayer
#' @export
#' @examples
#' bed_f = system.file("extdata/test_loading.bed",
#'     package = "seqsetvis", mustWork = TRUE)
#' easyLoad_bed(bed_f, "my_bed")
easyLoad_bed = function(file_paths, file_names = NULL, extraCols = character()){
    stopifnot(is.character(file_paths))
    stopifnot(all(file.exists(file_paths) | grepl("://", file_paths)))
    #check names of file-paths
    if(is.null(names(file_paths))){
        names(file_paths) = basename(file_paths)
    }
    #override names with file_names if not NULL
    if(!is.null(file_names)){
        stopifnot(is.character(file_paths))
        stopifnot(length(file_paths) == length(file_names))
        names(file_paths) = file_names
    }
    grs <- lapply(file_paths, function(f){
        rtracklayer::import(f, format = "BED",
                            extraCols = extraCols)
    })
    grs
}
