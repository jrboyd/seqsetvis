library(data.table)


# gtf_f = "~/R/GRIDseq_test/ref/Homo_sapiens.GRCh38.83.gtf"
paste_gtf.dt = function(gtf_f, core_attribs = c("gene_id", "gene_name"), additional_attribs = c()){
  # system.time({
  ref_dt = fread(gtf_f, sep = "\t")
  colnames(ref_dt) = c("seqnames", "source", "feature_type", "start", "end", "dot", "strand", "dot2", "attribs")
  ref_dt = ref_dt[feature_type == "gene"]
  attribs_dt = ref_dt[, tstrsplit(attribs, split = ";")]
  attribs_dt$index = 1:nrow(attribs_dt)
  attribs_dt = melt(attribs_dt, id.vars = "index", na.rm = T)
  attribs_dt$value = gsub("\"", "", attribs_dt$value)
  attribs_dt$value = gsub("^ ", "", attribs_dt$value)
  attribs_dt[, c("attrib_name", "attrib_value") := tstrsplit(value, split = " ", keep = 1:2)]
  attribs_dt[, c("variable", "value") := NULL]
  setkey(attribs_dt, index, attrib_name)
  parsed_dt = data.table(attribs_dt[.(1:nrow(ref_dt), core_attribs[1]), attrib_value])
  for(attrib in c(core_attribs[-1], additional_attribs)){
    parsed_dt = cbind(parsed_dt, data.table(attribs_dt[.(1:nrow(ref_dt), attrib), attrib_value]))
  }
  colnames(parsed_dt) = c(core_attribs, additional_attribs)
  out_dt = cbind(ref_dt[, c(1,4,5,7)], parsed_dt)
  return(out_dt)
  # })
}

.compare_parse = function(){
  gtf_f = "~/R/GRIDseq_test/ref/Homo_sapiens.GRCh38.83.gene.gtf"
  print("new and fast")
  print(system.time({dt1 = paste_gtf.dt(gtf_f)}))
  source("~/R/jrb_R_scripts/parse_gtf.R")
  print("old and slow")
  print(system.time({
    dt2 = as.data.table(parse_gtf(gtf_f, additional_attrib = "gene_biotype"))
    dt2[, gene_biotype := sub(";", "", gene_biotype)]
  }))
}

# cbind(table(dt1$gene_biotype), table(dt1$gene_biotype))
