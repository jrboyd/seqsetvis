
.ssvAnnotateSubjectGRanges = function(annotation_source, subject_gr, annotation_name = NULL){
    if(!is(annotation_source, "GRanges")){
        stop("GRanges expected as annotation_source.")
    }
    if(!is(subject_gr, "GRanges")){
        stop("GRanges expected as subject_gr.")
    }
    if(!is.null(annotation_name)){
        if(!is.character(annotation_name)){
            stop("character or NULL expected as annotation_name.")
        }
    }
    anno_olap = findOverlaps(query = annotation_source, subject = subject_gr)
    subject_hit_ids = names(subject_gr)[subjectHits(anno_olap)]
    anno_hit_mcols = mcols(annotation_source[queryHits(anno_olap)])

    df1 = cbind(DataFrame(id__ = subject_hit_ids), anno_hit_mcols)
    #append annotation_name if specified
    if(!is.null(annotation_name)){
        not_id = colnames(df1) != "id__"
        colnames(df1)[not_id] = paste0(colnames(df1)[not_id], ".", annotation_name)
    }

    df2 = mcols(subject_gr)
    df2$id__ = rownames(df2)
    #merge subject and annotation, restore rownames
    shared_cn = setdiff(intersect(colnames(df1), colnames(df2)), "id__")
    if(length(shared_cn) > 0){
        stop("Annotation shares colnames with subject_gr. Try specifying annotation_name to avoid this.\n  Offending colnames:\n",
             paste(shared_cn, collapse = "\n"))
    }
    dfm = merge(df2, df1, by = "id__", all.x = TRUE)
    rownames(dfm) = dfm$id__
    dfm$id__ = NULL
    mcols(subject_gr) =  dfm[names(subject_gr),]
    subject_gr
}

#' ssvAnnotateSubjectGRanges
#'
#' @param annotation_source A single GRanges, a list of GRanges, or a GRangesList
#' @param subject_gr The base GRanges to add annotation mcols to.
#' @param annotation_name Optional name for single GRanges. Required for list inputs if list does not have names.
#'
#' @return GRanges with the same regions as `subject_gr` but with addtional mcols added from annotation_source.
#' @export
#' @rdname ssvAnnotateSubjectGRanges
#'
#' @examples
#' library(GenomicRanges)
#' np_grs = CTCF_in_10a_narrowPeak_grs
#' olap_gr = ssvOverlapIntervalSets(np_grs)
#' # annotating with a signle GRanges is OK
#' ssvAnnotateSubjectGRanges(np_grs$MCF10A_CTCF, olap_gr)
#' # provide a name if that's useful
#' ssvAnnotateSubjectGRanges(np_grs$MCF10A_CTCF, olap_gr, annotation_name = "MCF10A")
#' # a named list adds each annotation
#' ssvAnnotateSubjectGRanges(np_grs, olap_gr)
#' # overriding list names is an option
#' ssvAnnotateSubjectGRanges(np_grs, olap_gr, LETTERS[1:3])
#' # GRangeList are handled like a standard list
#' ssvAnnotateSubjectGRanges(GRangesList(np_grs), olap_gr, LETTERS[1:3])
setGeneric("ssvAnnotateSubjectGRanges", function(annotation_source, subject_gr, annotation_name = NULL){ standardGeneric("ssvAnnotateSubjectGRanges")})

#' @rdname ssvAnnotateSubjectGRanges
#' @export
setMethod("ssvAnnotateSubjectGRanges", signature(annotation_source = "GRanges"), .ssvAnnotateSubjectGRanges)

#' @rdname ssvAnnotateSubjectGRanges
#' @export
setMethod("ssvAnnotateSubjectGRanges", signature(annotation_source = "list"), function(annotation_source, subject_gr, annotation_name = NULL){
    if(!all(sapply(annotation_source, is, class2 = "GRanges"))){
        stop("All items in annotation_source list must be GRanges.")
    }
    if(!is.null(annotation_name)){
        if(!is.character(annotation_name)){
            stop("character or NULL expected as annotation_name.")
        }
        if(!length(annotation_name) == length(annotation_source)){
            stop("annotation_source and annotation_name must be same length.")
        }
        names(annotation_source) = annotation_name
    }
    if(is.null(names(annotation_source))){
        stop("annotation_source must have names set.")
    }
    if(any(duplicated(names(annotation_source)))){
        stop("names of annotation_source must be unique.")
    }
    for(name in names(annotation_source)){
        subject_gr = ssvAnnotateSubjectGRanges(
            annotation_source = annotation_source[[name]],
            subject_gr = subject_gr,
            annotation_name = name
        )
    }
    subject_gr
})

#' @rdname ssvAnnotateSubjectGRanges
#' @export
setMethod("ssvAnnotateSubjectGRanges", signature(annotation_source = "GRangesList"),
          function(annotation_source, subject_gr, annotation_name = NULL){
              ssvAnnotateSubjectGRanges(
                  as.list(annotation_source),
                  subject_gr = subject_gr,
                  annotation_name = annotation_name)
          })
