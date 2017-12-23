# library(peakvisr)
library(testthat)
a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1))
b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1))
c = GRanges("chr1", IRanges(8:10*10+5, 8:10*10+6))

#test all plotting functions
method_names = c("setPlotVenn", "setPlotEuler", "setPlotBars", "setPlotPie", "setPlotHeatmap")
for(met_name in method_names){
  # test_that(paste(met_name , "exists"), {
  #   expect_true(exists(met_name))
  #   expect_failure(expect_error(get(met_name)))
  # })
  if(!exists(met_name)){
    warning(paste("function", met_name, "couldn't be found! not tested."))
    next
  }
  met = get(met_name)
  if(!is.function(met)){
    warning(paste("function", met_name, "wasn't a function! not tested."))
    next
  }


  test_that(paste(met_name, "accepts overlapIntervalSets output."), {

    olap = overlapIntervalSets(list("a" = a, "b" = b, "c" = c))
    p = met(olap)
    p
    expect_s3_class(p, class = "ggplot")
  })

  # test_that(paste(met_name, "various other paramters don't throw error"), {
  #   olap = overlapIntervalSets(list("a" = a, "b" = b, "c" = c))
  #   p = met(olap, circle.col = c("red", "blue", "green"), fill_alpha = .1,
  #              counts_txt_size = 10, show_outside_count = T,
  #              counts_as_labels = T)
  #   expect_s3_class(p, class = "ggplot")
  #   p = met(olap, circle.col = c("red", "blue", "green"), fill_circles = F,
  #              counts_txt_size = 10, show_outside_count = T,
  #              counts_as_labels = T)
  #   expect_s3_class(p, class = "ggplot")
  # })

  test_that(paste(met_name, "accepts various membership table types."), {
    olap = overlapIntervalSets(list("a" = a, "b" = b, "c" = c))
    DF = mcols(olap) #DataFrame
    p = met(DF)
    p
    expect_s3_class(p, class = "ggplot")
    df = as.data.frame(DF) #data.frame
    p = met(df)
    p
    expect_s3_class(p, class = "ggplot")
    mat = as.matrix(df) #matrix
    p = met(mat)
    p
    expect_s3_class(p, class = "ggplot")
  })

  test_that(paste(met_name, "accepts various membership table types. missing names."), {
    olap = overlapIntervalSets(list("a" = a, "b" = b, "c" = c))
    DF = mcols(olap) #DataFrame

    df = as.data.frame(DF) #data.frame
    colnames(df) = NULL
    p = met(df)
    p
    expect_s3_class(p, class = "ggplot")
    mat = as.matrix(df) #matrix
    p = met(mat)
    p
    expect_s3_class(p, class = "ggplot")
  })

  test_that(paste(met_name, "accepts set lists. names present and missing"), {
    olap = overlapIntervalSets(list("a" = a, "b" = b, "c" = c))
    DF = mcols(olap) #DataFrame
    df = as.data.frame(DF) #data.frame
    rownames(df) = 1:nrow(df)
    sets = apply(df, 2, function(x)rownames(df)[x])
    p = met(sets)
    p
    expect_s3_class(p, class = "ggplot")
    names(sets) = NULL #missing names
    p = met(sets)
    p
    expect_s3_class(p, class = "ggplot")
  })

}


