
setGeneric("setPlotMakeMT", function(object){
  standardGeneric("setPlotMakeMT")
})

setMethod("setPlotMakeMT", signature(object = "list"), function(object){
  if (all(sapply(object, class) == "GRanges")) {#GRanges are a special case
    print("handling list of GRanges like GRangeList.")
    # object = overlapIntervalSets(object)
    return(GRangesList(object))
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

setMethod("setPlotMakeMT", signature(object = "GRangesList"), function(object){
  setPlotMakeMT(overlapIntervalSets(object))
})

setMethod("setPlotMakeMT", signature(object = "GRanges"), function(object){
  object = elementMetadata(object)
  setPlotMakeMT(object)
})

setMethod("setPlotMakeMT", signature(object = "DataFrame"), function(object){
  object = as.data.frame(object)
  setPlotMakeMT(object)
})

setMethod("setPlotMakeMT", signature(object = "matrix"), function(object){
  object = as.data.frame(object)
  setPlotMakeMT(object)
})

setMethod("setPlotMakeMT", signature(object = "data.frame"), function(object){
  if (is.null(colnames(object))) {
    colnames(object) = paste0("set_", LETTERS[seq_len(ncol(object))])
  }
  object
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
