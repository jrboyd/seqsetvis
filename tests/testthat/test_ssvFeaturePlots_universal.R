testthat::context("FeaturePlots")
library(seqsetvis)
library(testthat)
library(GenomicRanges)
a = GRanges("chr1", IRanges(1:7*10, 1:7*10+1))
b = GRanges("chr1", IRanges(5:10*10, 5:10*10+1))
c = GRanges("chr1", IRanges(8:10*10+5, 8:10*10+6))

#test all plotting functions
method_names = c("ssvFeatureVenn", "ssvFeatureEuler", "ssvFeatureBars", "ssvFeaturePie", "ssvFeatureBinaryHeatmap", "ssvFeatureUpset")
for(met_name in method_names){
  if(!exists(met_name)){
    warning("function", met_name, "couldn't be found! not tested.")
    next
  }
  met = get(met_name)
  if(!is.function(met)){
    warning("function", met_name, "wasn't a function! not tested.")
    next
  }


  test_that(paste(met_name, "accepts ssvOverlapIntervalSets output."), {

    olap = ssvOverlapIntervalSets(list("a" = a, "b" = b, "c" = c))
    p = met(olap)
    p
    expect_s3_class(p, class = "ggplot")
  })

  test_that(paste(met_name, "accepts various membership table types."), {
    olap = ssvOverlapIntervalSets(list("a" = a, "b" = b, "c" = c))
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
    olap = ssvOverlapIntervalSets(list("a" = a, "b" = b, "c" = c))
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
    olap = ssvOverlapIntervalSets(list("a" = a, "b" = b, "c" = c))
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


