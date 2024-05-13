
.default_multi_resolver_FUN = function(x, variable.name = NULL){
    if(is.character(x)){
        x = paste0(sort(x)[1], "-from_", length(x))
    }else if(is.numeric(x)){
        x = max(x)
    }else if(is.factor(x)){
        x = sort(x)[1]
    }else if(is.logical(x)){
        x = any(x)
    }
    x
}

.ssvAnnotateSubjectGRanges = function(annotation_source,
                                      subject_gr,
                                      annotation_name = NULL,
                                      multi_resolver_FUN = "default"){
    if(!is(annotation_source, "GRanges")){
        stop("GRanges expected as annotation_source.")
    }
    if(!is(subject_gr, "GRanges")){
        stop("GRanges expected as subject_gr.")
    }
    if(is.null(names(subject_gr))){
        names(subject_gr) = seq_along(subject_gr)
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
    if(any(duplicated(df1$id__))){
        if(is.function(multi_resolver_FUN)){
            message("Resolving multiple overlapping annotatation_source items with provided function...")
        }else if(multi_resolver_FUN == "default"){
            message("Resolving multiple overlapping annotatation_source items with default function (override with multi_resolver_FUN)...")
            multi_resolver_FUN = .default_multi_resolver_FUN
        }else{
            stop('multi_resolver_FUN must be a user provided function or "default".')
        }
        resolve_grp = function(grp){
            resolve_var_name = function(var_name){
                multi_resolver_FUN(x = grp[[var_name]], variable.name = var_name)
            }
            cn_todo = colnames(grp[-1])
            names(cn_todo) = cn_todo
            lapply(cn_todo, resolve_var_name)
        }
        res = lapply(
            split(df1, df1$id__),
            resolve_grp
        )
        test_len = vapply(res, function(x){max(lengths(x))}, FUN.VALUE = 1)
        if(any(test_len > 1)){
            print(res[which(test_len > 1)[1]])
            stop("multi_resolver_FUN did not resolve to single values. The first problematic entry is printed above.")

        }
        res = lapply(res, as.data.frame)
        df1 = do.call(rbind, res)
        # df1 = DataFrame(df1)
        df1 = cbind(DataFrame(id__ = rownames(df1)), df1)

    }
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
#' @param multi_resolver_FUN Optional function to resolve multiple overlapping
#'   annotation source regions per subject region. This function must accept 2
#'   arguments. `x` is the values in a single mcol attribute and `variable.name`
#'   is the name of variable. A single value must be returned or an error will
#'   be generated. The default of "default" can handle numeric, logical,
#'   character, and factor types.
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
setGeneric("ssvAnnotateSubjectGRanges",
           function(annotation_source,
                    subject_gr,
                    annotation_name = NULL,
                    multi_resolver_FUN = "default"){
               standardGeneric("ssvAnnotateSubjectGRanges")
           })

#' @rdname ssvAnnotateSubjectGRanges
#' @export
setMethod("ssvAnnotateSubjectGRanges", signature(annotation_source = "GRanges"), .ssvAnnotateSubjectGRanges)

#' @rdname ssvAnnotateSubjectGRanges
#' @export
setMethod("ssvAnnotateSubjectGRanges", signature(annotation_source = "list"),
          function(annotation_source, subject_gr,
                   annotation_name = NULL, multi_resolver_FUN = "default"){
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
                      annotation_name = name,
                      multi_resolver_FUN = multi_resolver_FUN
                  )
              }
              subject_gr
          })

#' @rdname ssvAnnotateSubjectGRanges
#' @export
setMethod("ssvAnnotateSubjectGRanges", signature(annotation_source = "GRangesList"),
          function(annotation_source, subject_gr, annotation_name = NULL, multi_resolver_FUN = "default"){
              ssvAnnotateSubjectGRanges(
                  as.list(annotation_source),
                  subject_gr = subject_gr,
                  annotation_name = annotation_name,
                  multi_resolver_FUN = multi_resolver_FUN)
          })
