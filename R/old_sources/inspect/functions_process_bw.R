library(rtracklayer)
library(pbapply)
library(data.table)
fetch_windowed_bw = function(bw_file, win_size = 50, qgr = NULL){
  suppressWarnings({
  if(is.null(qgr)){
    bw_gr = import.bw(bw_file)  
  }else{
    bw_gr = import.bw(bw_file, which = qgr)
  }
  })
  windows = slidingWindows(qgr, width = win_size, step = win_size)
  print(object.size(windows), units = "GB")
  print(object.size(bw_gr), units = "GB")
  windows = unlist(windows)
  mid_gr = function(gr){
    start(gr) + floor((width(gr) - 1)/2)
  }
  mids = mid_gr(windows)
  start(windows) = mids
  end(windows) = mids
  olaps = suppressWarnings(as.data.table(findOverlaps(windows, bw_gr)))
  #patch up missing/out of bound data with 0
  missing_idx = setdiff(1:length(windows), olaps$queryHits)
  if(length(missing_idx) > 0){
    olaps = rbind(olaps, data.table(queryHits = missing_idx, subjectHits = length(bw_gr) + 1))[order(queryHits)]
    bw_gr = c(bw_gr, GRanges("chr1", IRanges(1, 1), score = 0))
  }
  #set FE and output
  windows$FE = bw_gr[olaps$subjectHits]$score
  return(windows)
}

#tolerates overlapping regions in qgr
fetch_windowed_bw_as_dt = function(bw_file, win_size = 50, qgr = NULL){
  suppressWarnings({
    if(is.null(qgr)){
      bw_gr = import.bw(bw_file)  
    }else{
      bw_gr = import.bw(bw_file, which = qgr)
    }
  })
  windows = slidingWindows(qgr, width = win_size, step = win_size)
  names(windows) = qgr$id
  print(object.size(windows), units = "GB")
  print(object.size(bw_gr), units = "GB")
  windows = unlist(windows)
  windows$id = names(windows)
  mid_gr = function(gr){
    start(gr) + floor((width(gr) - 1)/2)
  }
  mids = mid_gr(windows)
  start(windows) = mids
  end(windows) = mids
  olaps = suppressWarnings(as.data.table(findOverlaps(windows, bw_gr)))
  #patch up missing/out of bound data with 0
  missing_idx = setdiff(1:length(windows), olaps$queryHits)
  if(length(missing_idx) > 0){
    olaps = rbind(olaps, data.table(queryHits = missing_idx, subjectHits = length(bw_gr) + 1))[order(queryHits)]
    bw_gr = c(bw_gr, GRanges("chr1", IRanges(1, 1), score = 0))
  }
  #set FE and output
  # windows = windows[olaps$queryHits]
  windows$FE = bw_gr[olaps$subjectHits]$score
  bw_dt = as.data.table(windows)
  bw_dt[, x := start - min(start), by = id]
  invisible(bw_dt)
}

#extracts data.table, one row per qgr from bw_gr
bw_gr2dt = function(bw_gr, qgr, win_size = 50){
  bw_dt = as.data.table(bw_gr)
  ol = findOverlaps(bw_gr, qgr)
  bw_dt$hit = "tmp"
  bw_dt = bw_dt[queryHits(ol)]
  bw_dt$hit = as.character(qgr$id[subjectHits(ol)])
  bw_dt[, x := start - min(start), by = hit]
  bw_dt = bw_dt[x %% win_size == 0]
  return(bw_dt)
}

gg_bw_banded_quantiles = function(bw_dt, hsv_min = 0, hsv_max = .7, 
                                  n_quantile = 18, quantile_min = .05, 
                                  quantile_max = .95, is_centered = T, win_size = 50){
  #hsv_min = 0; hsv_max = .7; n_quantile = 18; quantile_min = .05; quantile_max = .95; is_centered = T; win_size = 50
  require(data.table)
  q2do = 0:n_quantile/n_quantile
  q2do = round(quantile_min + q2do * (quantile_max - quantile_min), digits = 3)
  dt = bw_dt[, .(qs = .(.(quantile(FE, q2do)))), by = x]
  dt = cbind(dt, as.data.table(t(sapply(dt$qs, function(x)x[[1]]))))
  dt$qs = NULL
  dtm = melt(dt, id.vars = "x")
  setkey(dtm, variable)
  q2do_str = paste0(q2do * 100, "%")
  dt_low = dtm[q2do_str[-length(q2do_str)], .(x, low_q = variable, low = value)]
  dt_high = dtm[q2do_str[-1], .(x, high_q = variable, high = value)]
  dt_c = cbind(dt_low, dt_high[, -1])
  dt_c = dt_c[, .(x, low, high, q_range = paste0(sub("%", "", low_q), "-", high_q))]
  q_o = unique(dt_c$q_range)
  dt_c$q_range = factor(dt_c$q_range, levels = q_o)
  if(!is.null(win_size)){
    dt_c[, xe := x + as.integer(win_size)]
    dt_c = melt(dt_c, id.vars = c("low", "high", "q_range"), value.name = "x")
    dt_c = dt_c[order(variable, decreasing = T)][order(x)][order(q_range)]
  }
  if(is.numeric(is_centered)){
    dt_c[, x := x - is_centered]
  }else if(is.logical(is_centered)){
    if(is_centered){
      dt_c[, x := x - median(x)]
    }
  }
  cols = rainbow(length(q_o), start = hsv_min, end = hsv_max)
  names(cols) = q_o
  dt_c$rn = 1:nrow(dt_c)
  dt_c[, c("q_low", "q_high") := tstrsplit(sub("%", "", q_range), split = "-", keep = 1:2), by = rn]
  dt_c[, q_num := (as.numeric(q_low) + as.numeric(q_high)) / 2]
  dt_c[, c("rn", "q_low", "q_high") := NULL]
  ggplot(dt_c) + 
    geom_ribbon(aes(x = x, ymin = low, ymax = high, fill = q_range)) +
    labs(fill = "quantile band", 
         y = "FE", x = "bp", 
         title = "Enrichment Profiles", 
         subtitle = "aggregated by quantile range") +
    scale_fill_manual(values = cols) +
    scale_color_manual(breaks = q2do, palette = cols) 
    
}
