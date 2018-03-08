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
        char_object = object
        object = set_list2memb(char_object)
    } else {
        stop(paste("can't handle list of non-character classes: ", paste(sapply(object, class), collapse = ", ")))
    }
    return(object)
})

#' GRangesList input
#' setOldClass("GRangesList")
#' @param object a list of GRanges, sent to overlapIntervalSets
#' @import GRangesList
setMethod("setPlotMakeMT", signature(object = "GRangesList"), function(object){
    GRlist_object = object
    setPlotMakeMT(overlapIntervalSets(GRlist_object))
})



#' mcols from GRange input
#' setOldClass("GRanges")
#' @param object A single GRanges, will use logical metadata columns
#' @import GRanges
setMethod("setPlotMakeMT", signature(object = "GRanges"), function(object){
    gr_object = object
    mc_object = mcols(gr_object)
    setPlotMakeMT(mc_object)
})



#' DataFrame input
#' setOldClass("DataFrame")
#' @param object A single DataFrame (typically from GRanges metadata), will use logical columns
#' @import DataFrame
setMethod("setPlotMakeMT", signature(object = "DataFrame"), function(object){
    DF_object = object
    df_object = as.data.frame(DF_object)
    setPlotMakeMT(df_object)
})

#' matrix input
#' @param object A matrix, sent to data.frame to handle column naming
#' @examples
#' memb_mat =  matrix(c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE), ncol = 2, byrow = FALSE)
#' setPlotMakeMT(memb_mat)
setMethod("setPlotMakeMT", signature(object = "matrix"), function(object){
    mat_object = object
    df_object = as.data.frame(mat_object)
    setPlotMakeMT(df_object)
})

#' data.frame input, final output
#' The final method for most inputs, checks column names
#' @param object a data.frame
#' @examples
#' memb_df = data.frame(a = c(TRUE, TRUE, FALSE, FALSE), b = c(TRUE, FALSE, TRUE, FALSE))
#' setPlotMakeMT(memb_df)
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
