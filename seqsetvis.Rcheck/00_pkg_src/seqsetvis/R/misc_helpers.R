#' converts a valid r color name ("black", "red", "white", etc.) to a hex value
#' @param color_name character. one or more r color names.
col2hex = function(color_name) {
  rgb(t(col2rgb(color_name))/255)
}

#' convert a list of sets, each list item should be a character vector
#' denoting items in sets
#' @param set_list a list of character vectors.  default names will be added if missing
set_list2memb = function(set_list) {
  if (is.null(names(set_list))) {
    names(set_list) = paste0("set_", LETTERS[seq_along(set_list)])
  }
  rn = unique(unlist(set_list))
  cn = names(set_list)
  memb = matrix(F, nrow = length(rn), ncol = length(cn))
  rownames(memb) = rn
  colnames(memb) = cn
  for (column in cn) {
    memb[set_list[[column]], column] = T
  }
  return(memb)
}
