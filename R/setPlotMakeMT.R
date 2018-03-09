#' generic for methods to convert various objects to a logical matrix indicating
#' membership of items (rows) in sets (columns)
#'
#' @param object the object to convert. Supported types: list (of character or GRanges),
#' GRanges with membership table metadata, GrangesList,
#' data.frame/matrix/DataFrame of membership table
#' @rdname setPlotMakeMT-methods
#' @exportMethod setPlotMakeMT
#' @return a logical matrix indicating membership of items (rows) in sets (columns)
setGeneric("setPlotMakeMT", function(object){
    standardGeneric("setPlotMakeMT")
})


#' list of character vectors input
#' @rdname setPlotMakeMT-methods
#' @aliases setPlotMakeMT,list-method
#' @import GenomicRanges
#' @examples
#' char_list = list(letters[1:3], letters[2:4])
#' setPlotMakeMT(char_list)
#' library(GenomicRanges)
#' gr_list = list(GRanges("chr1", IRanges(1:3*2, 1:3*2)), GRanges("chr1", IRanges(2:4*2, 2:4*2)))
#' setPlotMakeMT(gr_list)
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
#' @rdname setPlotMakeMT-methods
#' @aliases setPlotMakeMT,GRangesList-method
#' @import GenomicRanges
#' @examples
#' library(GenomicRanges)
#' gr_list = list(GRanges("chr1", IRanges(1:3*2, 1:3*2)), GRanges("chr1", IRanges(2:4*2, 2:4*2)))
#' setPlotMakeMT(GRangesList(gr_list))
setMethod("setPlotMakeMT", signature(object = "GRangesList"), function(object){
    GRlist_object = object
    setPlotMakeMT(overlapIntervalSets(GRlist_object))
})



#' mcols from GRange input
#' setOldClass("GRanges")
#' @rdname setPlotMakeMT-methods
#' @aliases setPlotMakeMT,GRanges-method
#' @import GenomicRanges
#' @examples
#' gr = GRanges("chr1", IRanges(1:3*2, 1:3*2))
#' gr$set_a = c(TRUE, TRUE, FALSE)
#' gr$set_b = c(FALSE, TRUE, TRUE)
#' setPlotMakeMT(gr)
setMethod("setPlotMakeMT", signature(object = "GRanges"), function(object){
    gr_object = object
    mc_object = mcols(gr_object)
    setPlotMakeMT(mc_object)
})



#' DataFrame input
#' setOldClass("DataFrame")
#' @rdname setPlotMakeMT-methods
#' @aliases setPlotMakeMT,DataFrame-method
#' @import GenomicRanges
#' @examples
#' gr = GRanges("chr1", IRanges(1:3*2, 1:3*2))
#' gr$set_a = c(TRUE, TRUE, FALSE)
#' gr$set_b = c(FALSE, TRUE, TRUE)
#' setPlotMakeMT(mcols(gr))
setMethod("setPlotMakeMT", signature(object = "DataFrame"), function(object){
    DF_object = object
    df_object = as.data.frame(DF_object)
    setPlotMakeMT(df_object)
})

#' matrix input
#' @rdname setPlotMakeMT-methods
#' @aliases setPlotMakeMT,matrix-method
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
#' @rdname setPlotMakeMT-methods
#' @aliases setPlotMakeMT,data.frame-method
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
