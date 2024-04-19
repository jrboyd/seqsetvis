#' generic for methods to convert various objects to a logical matrix indicating
#' membership of items (rows) in sets (columns)
#' @export
#' @param object the object to convert. Supported types: list (of character
#' or GRanges), GRanges with membership table metadata, GRangesList,
#' data.frame/matrix/DataFrame of membership table
#' @rdname ssvMakeMembTable-methods
#' @exportMethod ssvMakeMembTable
#' @import methods
#' @return a logical matrix indicating membership of items (rows) in sets
#' (columns)
setGeneric("ssvMakeMembTable", function(object){
    standardGeneric("ssvMakeMembTable")
})


#' list of character vectors input
#' @export
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,list-method
#' @import GenomicRanges
#' @import methods
#' @examples
#' char_list = list(letters[1:3], letters[2:4])
#' ssvMakeMembTable(char_list)
#' library(GenomicRanges)
#' gr_list = list(GRanges("chr1", IRanges(1:3*2, 1:3*2)),
#'     GRanges("chr1", IRanges(2:4*2, 2:4*2)))
#' ssvMakeMembTable(gr_list)
setMethod("ssvMakeMembTable", signature(object = "list"), function(object){
    #GRanges are a special case
    if (all(vapply(object, class, "character") == "GRanges")) {
        return(ssvMakeMembTable(GRangesList(object)))
    }
    if (any(vapply(object, class, "character") != "character")) {
        warning("Converting all non-character items to characters.")
        object = lapply(object, as.character)
    }
    if (all(vapply(object, class, "character") == "character")) {
        char_object = object
        object = set_list2memb(char_object)
    } else {
        #This may not be reachable anymore but I'm leaving as a backstop
        stop("Can't handle list of non-character classes as object: ",
             paste(vapply(object, class, "character"), collapse = ", "))
    }
    return(object)
})

#' GRangesList input
#' @export
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,GRangesList-method
#' @import GenomicRanges
#' @examples
#' library(GenomicRanges)
#' gr_list = list(GRanges("chr1", IRanges(1:3*2, 1:3*2)),
#'     GRanges("chr1", IRanges(2:4*2, 2:4*2)))
#' ssvMakeMembTable(GRangesList(gr_list))
setMethod("ssvMakeMembTable",
          signature(object = "GRangesList"), function(object){
              GRlist_object = object
              ssvMakeMembTable(ssvOverlapIntervalSets(GRlist_object))
          })



#' GRanges with mcols input
#' @export
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,GRanges-method
#' @import GenomicRanges
#' @examples
#' gr = GRanges("chr1", IRanges(1:3*2, 1:3*2))
#' gr$set_a = c(TRUE, TRUE, FALSE)
#' gr$set_b = c(FALSE, TRUE, TRUE)
#' ssvMakeMembTable(gr)
setMethod("ssvMakeMembTable", signature(object = "GRanges"), function(object){
    mc_object = GenomicRanges::mcols(object)
    rownames(mc_object) = names(object)
    ssvMakeMembTable(mc_object)
})


#' DataFrame input
#' @export
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
#' @export
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,matrix-method
#' @examples
#' memb_mat =  matrix(c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE),
#'     ncol = 2, byrow = FALSE)
#' ssvMakeMembTable(memb_mat)
setMethod("ssvMakeMembTable", signature(object = "matrix"), function(object){
    mat_object = object
    df_object = as.data.frame(mat_object)
    ssvMakeMembTable(df_object)
})

#' data.frame input, final output
#' The final method for all inputs, checks column names and returns logical
#' matrix
#' @export
#' @rdname ssvMakeMembTable-methods
#' @aliases ssvMakeMembTable,data.frame-method
#' @examples
#' memb_df = data.frame(a = c(TRUE, TRUE, FALSE, FALSE),
#'     b = c(TRUE, FALSE, TRUE, FALSE))
#' ssvMakeMembTable(memb_df)
setMethod("ssvMakeMembTable",
          signature(object = "data.frame"), function(object){
              col_is_logical = vapply(X = seq(ncol(object)),
                                   FUN = function(i)is.logical(object[,i]),
                                   FUN.VALUE = TRUE)
              if(!all(col_is_logical)){
                  message("Non-logical class columns have been dropped.")
                  object = object[, col_is_logical, drop = FALSE]
              }

              if (is.null(colnames(object))) {
                  colnames(object) = paste0("set_", LETTERS[seq_len(ncol(object))])
              }
              if(all(colnames(object) == paste0("V", seq_len(ncol(object))))){
                  colnames(object) = paste0("set_", LETTERS[seq_len(ncol(object))])
              }
              mat = as.matrix(object)
              rownames(mat) = rownames(object)
              return(mat)
          })

#' Convert any object accepted by ssvMakeMembTable to a factor
#' To avoid ambiguity,
#'
#' see \code{\link{ssvMakeMembTable}}
#'
#' @param object a valid object for conversion to a membership table and then
#' factor
#' @return a 2 column ("id" and "group") data.frame.
#' "id" is factor of item names if any or simply order of items.
#' "group" is a factor of set combinations
#' @export
#' @examples
#' ssvFactorizeMembTable(CTCF_in_10a_overlaps_gr)
#' ssvFactorizeMembTable(list(1:4, 2:3, 4:6))
ssvFactorizeMembTable = function(object){
    memb = ssvMakeMembTable(object)
    id_lev = names(rev(sort(rowSums(memb))))
    group_lev = names(rev(sort(colSums(memb))))
    grps = apply(memb, 1, function(x){
        paste(names(x)[x], collapse = " & ")
    })
    grp_counts = table(grps)
    grp_counts = rev(sort(grp_counts))
    grp_df = data.frame(id = names(grps), group = grps)
    rownames(grp_df) = NULL
    grp_df$id = factor(grp_df$id, levels = id_lev)
    grp_df$group = factor(grp_df$group, levels = names(grp_counts))
    return(grp_df)
}
