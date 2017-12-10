
filterModal <- function(failed = FALSE) {
  modalDialog(
    span('(Please filter the selected sample)'),
    DT::dataTableOutput("DTPeaksFilter", width = "auto"),
    if (failed)
      div(tags$b("Hey man, this didn't work! Not sure why.", style = "color: red;")),
    
    footer = tagList(fluidRow(
      column(width = 8,
             uiOutput("DTPeaksFilterElements")),
      column(width = 4,
             actionButton("BtnCancelFilter", "Cancel"),
             actionButton("BtnConfirmFilter", "Confirm")
      )
    )
    ),
    size = "l",
    title = "Filtering"
  )
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

server_filterModal = function(input, output, session, get_filtering_DF, set_filtering_DF, get_file_path, load_file = load_peak_wValidation){
  FilteringDF = reactiveVal(NULL)
  
  observeEvent(input$BtnFilterSet, {
    if(is.null(get_filtering_DF())){
      showNotification("No data selected.", type = "error")
      return()
    }
    #dataTable instance used to filter set data after adding
    FilteringDF(get_filtering_DF())
    showModal(filterModal())
  })
  observeEvent(input$BtnCancelFilter, {
    FilteringDF(NULL)
    removeModal()
  })
  output$DTPeaksFilter = DT::renderDataTable({
    df = FilteringDF()
    if(is.null(df)) return(NULL)
    DT::datatable(df, 
                  filter = list(position = "top", clear = TRUE, plain = F),
                  options = list(
                    scrollX = T,
                    pageLength = 10), rownames = F)
  })
  
  observeEvent(FilteringDF(), {
    df = FilteringDF()
    # showNotification("check update DTPeaksFilterElements")
    if(is.null(df)) return(NULL)
    output$DTPeaksFilterElements = renderUI({
      tagList(
        numericInput("NumFilterNumberOfRegions", "Desired Number of Regions", value = nrow(df), min = 0, max = nrow(df), step = 100),
        actionButton("BtnFilterByTruncate", label = "Truncate"),
        actionButton("BtnFilterByRandom", label = "Random Sample"),
        actionButton("BtnFilterReloadFile", label = "Reload File")
      )
    })
  })
  
  observeEvent(input$BtnFilterByRandom, {
    df = FilteringDF()
    df = df[sample(input$DTPeaksFilter_rows_all, size = input$NumFilterNumberOfRegions), ]
    FilteringDF(df)
  })
  observeEvent(input$BtnFilterByTruncate, {
    df = FilteringDF()
    df = df[input$DTPeaksFilter_rows_all[1:input$NumFilterNumberOfRegions],]
    FilteringDF(df)
  })
  observeEvent(input$BtnFilterReloadFile, {
    filepath = get_file_path()
    df = load_file(filepath)
    FilteringDF(df)
  })
  observeEvent(input$BtnConfirmFilter, {
    set_filtering_DF(FilteringDF()[input$DTPeaksFilter_rows_all,])
    FilteringDF(NULL)
    removeModal()
    # showNotification("Confirm filter.", type = "message")
  })
}
