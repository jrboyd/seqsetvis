# if(exists("all_res")) remove(all_res)
# 
# if(exists("all_res")){
#   start_i = length(unique(all_res$gene_id)) + 1
# }else{
#   start_i = 1
# } 
# res_list = pblapply(start_i:length(to_prof_gr), function(i){
#   my_gr = to_prof_gr[i]
#   chr = as.character(seqnames(my_gr))
#   s = start(my_gr) - 2000
#   e = end(my_gr) + 2000
#   res = add_bigwig_plot.tiles(bigwig_files = runx_bws, 
#                               chr = chr, start = s, end = e, 
#                               just_data = T, 
#                               n_bins = 200, 
#                               bin_method = "mean")
#   
#   res$index_o = rep(1:(nrow(res)/2), each = 2)
#   res$index_se = rep(1:2, nrow(res)/2)
#   # prof = res[, .(x = mean(x), y = unique(y)) , by = c("group", "id")]
#   
#   # res$id = NULL
#   res$gene_id = my_gr$gene_id
#   res$class = my_gr$gene_type
#   if(!exists("all_res")){
#     all_res <<- res
#   }else{
#     all_res <<- rbind(all_res, res)
#   }
#   # peak_dt = as.data.table(subsetByOverlaps(mda_grs, my_gr))
#   # peak_dt$y = -1
#   # peak_dt[group == "RUNX2"]$y = -2
#   # p = ggplot() + geom_line(data = res, aes(x = x, y = y, col = group)) + 
#   #   geom_segment(data = peak_dt, aes(x = start, xend = end, y = y, yend = y, col = group)) + 
#   #   # geom_line(data = prof, aes(x = x, y = y, col = group)) + 
#   #   coord_cartesian(ylim = c(-2,25))
#   # print(p)
# })
# # dev.off()
# # }
# save(all_res, file = "JG_profs.2kb.save")