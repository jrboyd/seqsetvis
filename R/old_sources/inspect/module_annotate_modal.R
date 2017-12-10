ref_path = "~/ShinyApps/shiny_peak_data/references"

annotateModal <- function(failed = FALSE) {
  modalDialog(
    h3(span("Unprocessed", style = "color: red;"), 'gene reference'),
    selectInput(inputId = "SelectAnnotateReference", 
                label = "Preloaded References", choices = dir(ref_path), selected = "hg38.gencode.v27",
                selectize = T),
    tags$hr(),
    DT::dataTableOutput("DTRawReferenceSet", width = "auto"),
    tags$hr(),
    h3(span("Processed", style = "color: red;"), 'gene reference'),
    fixedRow(
      column(width = 4,
             radioButtons("RadioAnnotateFeature", label = "Feature Type", choices = c("promoter", "gene body"), selected = "promoter")
      ),
      column(width = 4,
             conditionalPanel("input.RadioAnnotateFeature == 'promoter'", 
                              numericInput("NumUpstreamExtension", "Upstream TSS Extension", value = 2*10^3, min = 0, max = 10^5, step = 500),
                              numericInput("NumDownstreamExtension", "Downstream TSS Extension", value = 2*10^3, min = 0, max = 10^5, step = 500)
             )
      ),
      column(width = 4
             
      )
    ),
    tags$hr(),
    DT::dataTableOutput("DTProcessedReferenceSet", width = "auto"),
    tags$hr(),
    h3("Annotated features", style = "color: red;"),
    numericInput("NumAnnotateMaxDistance", "Max Annotation Distance", value = 10^5, min = 0, max = 10^9, step = 10000),
    checkboxInput("CheckDiscard", label = "Discard Un-annotated?", value = F),
    tags$hr(),
    DT::dataTableOutput("DTAnnotatedSet", width = "auto"),
    if (failed)
      div(tags$b("Hey man, this didn't work! Not sure why.", style = "color: red;")),
    
    footer = tagList(#fluidRow(
      # column(width = 6,
      #        
      #        
      # ),
      # column(width = 3,
      #        actionButton("BtnAnnotate", label = "Annotate")
      # ),
      # column(width = 3,
      actionButton("BtnCancelAnnotate", "Cancel"),
      actionButton("BtnConfirmAnnotate", "Confirm")
      # )
      #)
    ),
    size = "l",
    title = "Annotateing"
  )
}

server_annotateModal = function(input, output, session, get_annotateing_DF, set_annotateing_DF, roots_reference, load_file = load_peak_wValidation){
  RawRefDT = reactiveVal(NULL)
  ProcRefDT = reactiveVal(NULL)
  AnnotatedDT = reactiveVal(NULL)
  annotateModal()
  observeEvent(input$BtnAnnotateSet, {
    if(is.null(get_annotateing_DF())){
      showNotification("No data selected.", type = "error")
      return()
    }
    showModal(annotateModal())
  })
  observeEvent(input$SelectAnnotateReference, {
    ref_file = paste(ref_path, input$SelectAnnotateReference, sep = "/")
    if(!file.exists(ref_file)) return()
    load(paste(ref_path, input$SelectAnnotateReference, sep = "/"))
    df = ref_dt
    RawRefDT(df)
  })
  observeEvent({
    RawRefDT()
    input$DTRawReferenceSet_rows_all
    input$RadioAnnotateFeature
    input$NumUpstreamExtension
    input$NumDownstreamExtension
  }, {
    if(is.null(RawRefDT())){
      ProcRefDT(NULL)
      return()
    } 
    dt = RawRefDT()[input$DTRawReferenceSet_rows_all]
    # dt = copy(ref_dt)
    if(nrow(dt) == 0) return(data.frame())
    if(input$RadioAnnotateFeature == "promoter"){
      up_ex = as.integer(input$NumUpstreamExtension)
      down_ex = as.integer(input$NumDownstreamExtension)
      dt[strand == "+", end := start]
      dt[strand == "+", start := start - up_ex]
      dt[strand == "+", end := end + down_ex]
      dt[strand == "-", start := end]
      dt[strand == "-", start := start - down_ex]
      dt[strand == "-", end := end + up_ex]
    }else if(input$RadioAnnotateFeature == "gene body"){
      
    }
    ProcRefDT(dt)
  })
  
  observeEvent({
    ProcRefDT()
    input$NumAnnotateMaxDistance
    input$CheckDiscard
  }, {
    if(is.null(ProcRefDT())){
      AnnotatedDT(NULL)
      return()
    }
    if(is.na(input$NumAnnotateMaxDistance)){
       return()
    }
    ref_dt = ProcRefDT()
    fgr = get_annotateing_DF()
    discard_misses = input$CheckDiscard
    max_dist = input$NumAnnotateMaxDistance
    
    save(ref_dt, fgr, discard_misses, max_dist, file = "last_annotation.save")
    load('last_annotation.save')
    qgr = GRanges(fgr)
    ref_gr = GRanges(ref_dt)
    dists = distanceToNearest(qgr, ref_gr)
    fgr$gene_name = ref_dt[subjectHits(dists)]$gene_name
    fgr$gene_type = ref_dt[subjectHits(dists)]$gene_type
    fgr$distance = elementMetadata(dists)$distance
    over_max = fgr$distance > max_dist
    if(discard_misses){
      fgr = fgr[!over_max,]
    }else{
      fgr[over_max,]$gene_name = "no_hit"
      fgr[over_max,]$gene_type = "no_hit"
    }
    # distance2str = function(dist){
    #   as.
    # }
    dt = data.table(distance = fgr$distance, str = "")
    dt[distance < 10^3, str := paste(distance, "bp")]
    dt[distance >= 10^3 & distance < 10^4, str := paste(round(distance / 10^3,2), "kbp")]
    dt[distance >= 10^4 & distance < 10^5, str := paste(round(distance / 10^3,1), "kbp")]
    dt[distance >= 10^5 & distance < 10^6, str := paste(round(distance / 10^3,0), "kbp")]
    dt[distance >= 10^6 & distance < 10^7, str := paste(round(distance / 10^6,2), "Mbp")]
    dt[distance >= 10^7 & distance < 10^8, str := paste(round(distance / 10^6,1), "Mbp")]
    dt[distance >= 10^8 & distance < 10^9, str := paste(round(distance / 10^6,0), "Mbp")]
    dt[distance >= 10^9, str := paste(round(distance / 10^9,2), "Mbp")]
    fgr$dist_str = dt$str
    fdt = as.data.table(fgr)
    fdt[, i := 1:nrow(fdt)]
    fdt[, uniq_str := paste0(1:.N, "/", .N), by = gene_name]
    # fdt[, .(i, str = paste0(1:.N, "/", .N)), by = gene_name]
    # as.datafdt[order(i),]
    fdt[, id := paste0(gene_name, " ", uniq_str, " ", dist_str)]
    AnnotatedDT(fdt[order(i),])
  })
  
  observeEvent(input$BtnCancelAnnotate, {
    RawRefDT(NULL)
    ProcRefDT(NULL)
    AnnotatedDT(NULL)
    removeModal()
  })
  output$DTRawReferenceSet = DT::renderDataTable({
    # df = 
    if(is.null(df)) return(NULL)
    DT::datatable(RawRefDT(), 
                  filter = list(position = "top", clear = TRUE, plain = F),
                  options = list(
                    scrollX = T,
                    pageLength = 5), rownames = F)
  })
  
  output$DTProcessedReferenceSet = DT::renderDataTable({
    DT::datatable(ProcRefDT(), 
                  filter = list(position = "top", clear = TRUE, plain = F),
                  options = list(
                    scrollX = T,
                    pageLength = 5), rownames = F)
  })
  
  output$DTAnnotatedSet = DT::renderDataTable({
    # df = RawRefDT()
    # if(is.null(df)) return(NULL)
    DT::datatable(AnnotatedDT(), 
                  filter = list(position = "top", clear = TRUE, plain = F),
                  options = list(
                    scrollX = T,
                    pageLength = 5), rownames = F)
  })
  
  
  observeEvent({
    input$RadioAnnotateFeature
    input$NumUpstreamExtension
    input$NumDownsreamExtension
  }, {
    showNotification("updated processed annotation")
  })
  
  shinyFileChoose(input, 'FilesAnnotationReference', 
                  roots= roots_reference, 
                  filetypes=c("bed", "gff", "gtf"))
  
  # observeEvent(RawRefDT(), {
  #   df = RawRefDT()
  #   # showNotification("check update DTRawReferenceSetElements")
  #   if(is.null(df)) return(NULL)
  #   output$DTRawReferenceSetElements = renderUI({
  #     tagList(
  #       # shinyFilesButton("FilesAnnotationReference", label = "Find Reference", title = "Find bed/gtf/gff with reference info.", multiple = F),
  #       
  #     )
  #   })
  # })
  # 
  observeEvent(input$BtnAnnotate, {
    showNotification("BtnAnnotate", type = "message")
    df = RawRefDT()
    df = df[sample(input$DTRawReferenceSet_rows_all, size = input$NumAnnotateNumberOfRegions), ]
    RawRefDT(df)
  })
  observeEvent(input$BtnConfirmAnnotate, {
    set_annotateing_DF(AnnotatedDT()[input$DTAnnotatedSet_rows_all,])
    RawRefDT(NULL)
    # ProcRefDT(NULL)
    # AnnotatedDT(NULL)
    removeModal()
    # showNotification("Confirm annotate.", type = "message")
  })
}


paste_gtf.dt = function(gtf_f, filtered_type = "gene", core_attribs = c("gene_id", "gene_name"), additional_attribs = c()){
  require(data.table)
  ref_dt = fread(gtf_f, sep = "\t")
  colnames(ref_dt) = c("seqnames", "source", "feature_type", "start", "end", "dot", "strand", "dot2", "attribs")
  ref_dt = ref_dt[feature_type == filtered_type]
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
}

#assumed colnames for MACS2 peak files
peak_cn = c("seqnames", "start", "end", "id", "score", "strand", "FE", "p-value", "q-value", "summit_pos")
get_col_classes.df = function(df){
  sapply(1:ncol(df), function(i)class(df[[i]]))
}
get_col_classes = function(file, skipFirst = F){
  df = read.table(file, stringsAsFactors = F, header = F, nrows = 10, skip = ifelse(skipFirst, 1, 0))
  get_col_classes.df(df)
}
test_for_header = function(peak_file){
  withFirst = get_col_classes(peak_file, skipFirst = F)
  skipFirst = get_col_classes(peak_file, skipFirst = T)
  #definitely no header
  if(all(withFirst == skipFirst)){
    return(F)
  }
  if(all(withFirst == "character")){
    return(T)
  }
  warning(paste("can't determine if", peak_file, "has a header, guess not."))
  return(F)
}
load_peak_wValidation = function(peak_file, with_notes = F){
  has_header = test_for_header(peak_file)
  df = read.table(peak_file, stringsAsFactors = F, header = has_header)
  col_classes = get_col_classes.df(df)
  if(!has_header){
    if(ncol(df) == length(peak_cn)){
      if(with_notes){
        showNotification("assuming file is narrowPeak.", type = "warning")
      }else{
        print("assuming file is narrowPeak.")
      }
      colnames(df) = peak_cn  
    }else{
      if(with_notes){
        showNotification("file not narrowPeak, loading as minimal bed file.", type = "warning")
      }else{
        print("file not narrowPeak, loading as minimal bed file.")
      }
      bed_cn = peak_cn[1:4]
      nc = min(ncol(df), length(bed_cn))
      colnames(df)[1:nc] = peak_cn[1:nc]
    }
  }else{#try to make colnames GRanges compatible
    if(all(col_classes[1:5] == c("character", "integer", "integer", "integer", "character"))){
      if(with_notes){
        showNotification("file looks like saved GRanges.", type = "warning")
      }else{
        print("file looks like saved GRanges.")
      }
      colnames(df)[1:5] = c("seqnames", "start", "end", "width", "strand")
    }else if(all(col_classes[1:ncol(df)] == c("character", "integer", "integer", 
                                              "character", "integer", "character", 
                                              "numeric", "numeric", "numeric", "integer")[1:ncol(df)])){
      if(with_notes){
        showNotification("file looks like bed or encode peak", type = "warning")
      }else{
        print("file looks like bed or encode peak")
      }
    }else{
      colnames(df)[1:3] = c("seqnames", "start", "end")
      if(with_notes){
        showNotification("forced to assume first 3 columns are minimal bed, might break.", type = "warning")
      }else{
        print("forced to assume first 3 columns are minimal bed, might break.")
      }
    }
  }
  print(head(df))
  print(paste0(nrow(df), " total rows..."))
  invisible(df)
}