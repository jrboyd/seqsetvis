#set aside first
#iterate others and assign set name to hits in first
#remaining set name is no_hit
#return annotatated first gr
intersect_serial = function(grs){
  first_gr = grs[[1]]
  grs = grs[-1]
  first_gr$annotation = "no_hit"
  for(i in 1:length(grs)){
    i_gr = grs[[i]]
    olaps = findOverlaps(first_gr, i_gr)
    first_gr[queryHits(olaps)]$annotation = names(grs)  [i]
  }
  # pie(table(first_gr$annotation))
  # as.data.frame(first_gr)
  return(first_gr)
}

#merge and reduce all
#create membership table, all F, ncol = length(grs), nrow = length(reduced)
#iterate each i_gr in grs, memb[hits_reduced ,i] = T
#cbind membership and reduced then return
intersect_flat = function(grs){
  rdc_gr = reduce(unlist(GRangesList(grs)))
  for(i in 1:length(grs)){
    nam = names(grs)[i]
    elementMetadata(rdc_gr)[[nam]] = F
    olaps = findOverlaps(rdc_gr, grs[[i]])
    elementMetadata(rdc_gr)[[nam]][queryHits(olaps)] = T
  }
  # limma::vennDiagram(elementMetadata(rdc_gr)[,1:3], cex = .8)
  # limma::vennDiagram(elementMetadata(rdc_gr)[,c(2,4)], cex = .8)
  # mat = ifelse(as.matrix(elementMetadata(rdc_gr)), 1, 0)
  # mat = mat + runif(length(mat), min = -.02, .02)
  # gplots::heatmap.2(mat[sample(nrow(mat), 5000),]+.2, col = c("white", "black"), trace = "n")
  return(rdc_gr)
}

# load('peak_gr.save')
# grs = peak_gr

intersectR = function(grs, ext = 0, use_first = F){
  save(grs, ext, use_first, file = "last_intersectR.save")
  if(use_first){
    base_gr = grs[[1]]
    elementMetadata(base_gr) = NULL
    grs = grs[-1]
  }else{
    require(magrittr)
    base_gr = lapply(grs, function(x){elementMetadata(x) = NULL; x}) %>% GRangesList %>% unlist %>% reduce
  }
  start(base_gr) = start(base_gr) - ext
  end(base_gr) = end(base_gr) + ext
  base_gr = reduce(base_gr)
  start(base_gr) = start(base_gr) + ext
  end(base_gr) = end(base_gr) - ext
  suppressWarnings({
    for(i in 1:length(grs)){
      nam = names(grs)[i]
      elementMetadata(base_gr)[[nam]] = F
      olaps = findOverlaps(base_gr, grs[[i]])
      elementMetadata(base_gr)[[nam]][queryHits(olaps)] = T
    }
    if(use_first){
      base_gr$no_hit = apply(elementMetadata(base_gr), 1, function(x)all(!x))
    }
    base_gr$group = "no_hit"
    for(i in rev(1:length(grs))){
      i_gr = grs[[i]]
      olaps = findOverlaps(base_gr, i_gr)
      base_gr[queryHits(olaps)]$group = names(grs)  [i]
    }
  })
  return(base_gr)
}


