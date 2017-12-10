
filterModal <- function(failed = FALSE) {
  modalDialog(
    span('(Please filter the selected sample)'),
    DT::dataTableOutput("DTPeaksFilter", width = "auto"),
    if (failed)
      div(tags$b("Hey man, this didn't work! Not sure why.", style = "color: red;")),
    
    footer = tagList(
      uiOutput("DTPeaksFilterElements"),
      # actionButton("BtnFilterByNumber", label = "Truncate"),
      actionButton("BtnCancelFilter", "Cancel"),
      actionButton("BtnConfirmFilter", "Confirm")
    ),
    size = "l",
    title = "Filtering"
  )
}

server_filtering_modal = function(input, output, session, df_reactiveVal){
  FilteringDF = reactiveVal(NULL)
  
  observeEvent(input$BtnFilterSet, {
    # if(length(input$ChooseIntervalSets$selected) != 1){
    #   showNotification("No data selected in Step 2 right panel.", type = "error")
    #   return()
    # }
    #dataTable instance used to filter set data after adding
    df = df_reactiveVal()
    FilteringDF(as.data.frame(df))
    showModal(filterModal())
  })
  observeEvent(input$BtnCancelFilter, {
    FilteringDF(NULL)
    removeModal()
  })
  output$DTPeaksFilter = DT::renderDataTable({
    df = FilteringDF()
    if(is.null(df)) return(NULL)
    # sdf = rbind(head(df), rep(".", ncol(df)), tail(df))
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
        numericInput("NumFilterNumberOfRegions", "Truncate Number of Regions", value = nrow(df), min = 0, max = nrow(df), step = 100),
        actionButton("BtnFilterByNumber", label = "Truncate"),
        actionButton("BtnFilterReloadFile", label = "Reload File")
      )
    })
  })
  
  observeEvent(input$BtnFilterByNumber, {
    df = FilteringDF()
    df = df[input$DTPeaksFilter_rows_all[1:input$NumFilterNumberOfRegions],]
    FilteringDF(df)
    # output$DTPeaksFilter = DT::renderDataTable({
    #   # sdf = rbind(head(df), rep(".", ncol(df)), tail(df))
    #   DT::datatable(df, 
    #                 filter = list(position = "top", clear = TRUE, plain = F),
    #                 options = list(
    #                   scrollX = T,
    #                   pageLength = 10), rownames = F)
    # })
    
  })
  observeEvent(input$BtnFilterReloadFile, {
    filepath = SetsLoaded_FilePaths()[[input$ChooseIntervalSets$selected]]
    df = load_peak_wValidation(filepath)
    FilteringDF(df)
    #not necessary if setting df triggers update
    # output$DTPeaksFilter = DT::renderDataTable({
    #   # sdf = rbind(head(df), rep(".", ncol(df)), tail(df))
    #   DT::datatable(df, 
    #                 filter = list(position = "top", clear = TRUE, plain = F),
    #                 options = list(
    #                   scrollX = T,
    #                   pageLength = 10), rownames = F)
    # })
  })
  observeEvent(input$BtnConfirmFilter, {
    # tmp = SetsLoaded_DataFrames()
    # df = tmp[[input$ChooseIntervalSets$selected]]
    # tmp[[input$ChooseIntervalSets$selected]] = 
    df_reactiveVal(FilteringDF())
    # SetsLoaded_DataFrames(tmp)
    FilteringDF(NULL)
    removeModal()
    # showNotification("Confirm filter.", type = "message")
  })
}

