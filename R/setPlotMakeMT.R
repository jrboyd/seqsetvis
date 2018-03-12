#' generic for methods to convert various objects to a logical matrix indicating
#' membership of items (rows) in sets (columns)
#'
#' @param object the object to convert. Supported types: list (of character or GRanges),
#' GRanges with membership table metadata, GrangesList,
#' data.frame/matrix/DataFrame of membership table
#' @rdname ssvMakeMembTable-methods
#' @exportMethod ssvMakeMembTable
#' @return a logical matrix indicating membership of items (rows) in sets (columns)
setGeneric("ssvMakeMembTable", function(object){
    standardGeneric("ssvMakeMembTable")
})


#' list of character vectors input
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,list-method
#' @import GenomicRanges
#' @examples
#' char_list = list(letters[1:3], letters[2:4])
#' ssvMakeMembTable(char_list)
#' library(GenomicRanges)
#' gr_list = list(GRanges("chr1", IRanges(1:3*2, 1:3*2)), GRanges("chr1", IRanges(2:4*2, 2:4*2)))
#' ssvMakeMembTable(gr_list)
setMethod("ssvMakeMembTable", signature(object = "list"), function(object){
    if (all(sapply(object, class) == "GRanges")) {#GRanges are a special case
        # print("handling list of GRanges like GRangeList.")
        # object = overlapIntervalSets(object)
        return(ssvMakeMembTable(GRangesList(object)))
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
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,GRangesList-method
#' @import GenomicRanges
#' @examples
#' library(GenomicRanges)
#' gr_list = list(GRanges("chr1", IRanges(1:3*2, 1:3*2)), GRanges("chr1", IRanges(2:4*2, 2:4*2)))
#' ssvMakeMembTable(GRangesList(gr_list))
setMethod("ssvMakeMembTable", signature(object = "GRangesList"), function(object){
    GRlist_object = object
    ssvMakeMembTable(overlapIntervalSets(GRlist_object))
})



#' mcols from GRange input
#' setOldClass("GRanges")
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,GRanges-method
#' @import GenomicRanges
#' @examples
#' gr = GRanges("chr1", IRanges(1:3*2, 1:3*2))
#' gr$set_a = c(TRUE, TRUE, FALSE)
#' gr$set_b = c(FALSE, TRUE, TRUE)
#' ssvMakeMembTable(gr)
setMethod("ssvMakeMembTable", signature(object = "GRanges"), function(object){
    gr_object = object
    mc_object = mcols(gr_object)
    ssvMakeMembTable(mc_object)
})



#' DataFrame input
#' setOldClass("DataFrame")
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,DataFrame-method
#' @import GenomicRanges
#' @examples
#' gr = GRanges("chr1", IRanges(1:3*2, 1:3*2))
#' gr$set_a = c(TRUE, TRUE, FALSE)
#' gr$set_b = c(FALSE, TRUE, TRUE)
#' ssvMakeMembTable(mcols(gr))
setMethod("ssvMakeMembTable", signature(object = "DataFrame"), function(object){
    DF_object = object
    df_object = as.data.frame(DF_object)
    ssvMakeMembTable(df_object)
})

#' matrix of logicals, membership table
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,matrix-method
#' @examples
#' memb_mat =  matrix(c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE), ncol = 2, byrow = FALSE)
#' ssvMakeMembTable(memb_mat)
setMethod("ssvMakeMembTable", signature(object = "matrix"), function(object){
    mat_object = object
    df_object = as.data.frame(mat_object)
    ssvMakeMembTable(df_object)
})

#' data.frame input, final output
#' The final method for most inputs, checks column names
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,data.frame-method
#' @examples
#' memb_df = data.frame(a = c(TRUE, TRUE, FALSE, FALSE), b = c(TRUE, FALSE, TRUE, FALSE))
#' ssvMakeMembTable(memb_df)
setMethod("ssvMakeMembTable", signature(object = "data.frame"), function(object){
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
