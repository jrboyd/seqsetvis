#' generic for methods to convert various objects to a logical matrix indicating
#' membership of items (rows) in sets (columns)
#'
#' @param object the object to convert
#' @return a logical matrix indicating membership of items (rows) in sets (columns)
setGeneric("setPlotMakeMT", function(object){
  standardGeneric("setPlotMakeMT")
})

#' list of character vectors input
#' @param object a list of items that are valid for as.character
#' A list of GRanges are a special case and will be handled as a GRangesList
setMethod("setPlotMakeMT", signature(object = "list"), function(object){
  if (all(sapply(object, class) == "GRanges")) {#GRanges are a special case
    # print("handling list of GRanges like GRangeList.")
    # object = overlapIntervalSets(object)
    return(setPlotMakeMT(GRangesList(object)))
  }
  if (all(sapply(object, class) != "character")) {
    object = lapply(object, as.character)
  }
  if (all(sapply(object, class) == "character")) {
    object = set_list2memb(object)
  } else {
    stop(paste("can't handle list of non-character classes: ", paste(sapply(object, class), collapse = ", ")))
  }
  return(object)
})

#' GRangesList input
#' @param object a list of GRanges, sent to overlapIntervalSets
setMethod("setPlotMakeMT", signature(object = "GRangesList"), function(object){
  setPlotMakeMT(overlapIntervalSets(object))
})

#' mcols from GRange input
#' @param object A single GRanges, will use logical metadata columns
setMethod("setPlotMakeMT", signature(object = "GRanges"), function(object){
  object = mcols(object)
  setPlotMakeMT(object)
})

#' DataFrame input
#' @param object A single DataFrame (typically from GRanges metadata), will use logical columns
setMethod("setPlotMakeMT", signature(object = "DataFrame"), function(object){
  object = as.data.frame(object)
  setPlotMakeMT(object)
})

#' matrix input
#' @param object A matrix, sent to data.frame to handle column naming
setMethod("setPlotMakeMT", signature(object = "matrix"), function(object){
  object = as.data.frame(object)
  setPlotMakeMT(object)
})

#' data.frame input, final output
#' The final method for most inputs, checks column names
#' @param object a data.frame
setMethod("setPlotMakeMT", signature(object = "data.frame"), function(object){
  if (is.null(colnames(object))) {
    colnames(object) = paste0("set_", LETTERS[seq_len(ncol(object))])
  }
  if(all(colnames(object) == paste0("V", 1:ncol(object)))){
    colnames(object) = paste0("set_", LETTERS[seq_len(ncol(object))])
  }
  mat = as.matrix(object)
  rownames(mat) = rownames(object)
  return(mat)
})
#
# ###TODO this is suitable for a S4 generic method
# makeMembershipMatrix = function(object) {
#   if (class(object) == "GRanges") {
#     object = elementMetadata(object)
#   }
#   if (class(object) == "DataFrame") {
#     object = as.data.frame(object)
#   }
#
#   if (class(object) == "list") {
#     if (all(sapply(object, class) == "factor")) {
#       object = lapply(object, as.character)
#     }
#     if (all(sapply(object, class) == "character")) {
#       object = set_list2memb(object)
#     } else {
#       stop(paste("can't handle list of non-character classes: ", paste(sapply(object, class), collapse = ", ")))
#     }
#
#   }
#   if (is.matrix(object)) {
#     object = as.data.frame(object)
#   }
#   if (is.null(colnames(object))) {
#     colnames(object) = paste0("set_", LETTERS[seq_len(ncol(object))])
#   }
#   assertthat::assert_that(class(object) == "data.frame")
#   return(object)
# }
