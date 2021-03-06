% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_easyLoad.R
\name{easyLoad_bed}
\alias{easyLoad_bed}
\title{easyLoad_bed takes a character vector of file paths to bed plus files and
returning named list of GRanges.}
\usage{
easyLoad_bed(
  file_paths,
  file_names = NULL,
  extraCols = character(),
  n_cores = getOption("mc.cores", 1)
)
}
\arguments{
\item{file_paths}{character vector of paths to narrowPeak files.  If named,
those names will be used in output unless overriden by providing file_names.}

\item{file_names}{character vector of names for output list.  If not NULL
will override any existing names for file_paths.  Default is NULL.}

\item{extraCols}{named character vector of classes.  passed to
rtracklayer::import for format = "BED". default is character().}

\item{n_cores}{number of cores to use, uses mc.cores option if set or 1.}
}
\value{
a named list of GRanges loaded from file_paths
}
\description{
Mainly a utility function for loading MACS2 narrowPeak and broadPeak.
}
\examples{
bed_f = system.file("extdata/test_loading.bed",
    package = "seqsetvis", mustWork = TRUE)
easyLoad_bed(bed_f, "my_bed")
}
