% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_easyLoad.R
\name{easyLoad_seacr}
\alias{easyLoad_seacr}
\title{easyLoad_seacr takes a character vector of file paths to seacr output bed files and returns a named list of GRanges.}
\usage{
easyLoad_seacr(
  file_paths,
  file_names = NULL,
  n_cores = getOption("mc.cores", 1)
)
}
\arguments{
\item{file_paths}{character vector of paths to seacr bed files.  If named,
those names will be used in output unless overriden by providing file_names.}

\item{file_names}{character vector of names for output list.  If not NULL
will override any existing names for file_paths.  Default is NULL.}

\item{n_cores}{number of cores to use, uses mc.cores option if set or 1.}
}
\value{
a named list of GRanges loaded from file_paths
}
\description{
easyLoad_seacr takes a character vector of file paths to seacr output bed files and returns a named list of GRanges.
}
\examples{
bed_f = system.file("extdata/test_loading.seacr.bed",
    package = "seqsetvis", mustWork = TRUE)
easyLoad_seacr(bed_f, "my_seacr")
}
