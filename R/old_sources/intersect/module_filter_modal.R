
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

server_filterModal = function(input, output, session, get_filtering_DF, set_filtering_DF, get_file_path){
  FilteringDF = reactiveVal(NULL)
  
  observeEvent(input$BtnFilterSet, {
    if(length(input$ChooseIntervalSets$selected) != 1){
      showNotification("No data selected in Step 2 right panel.", type = "error")
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
    df = load_peak_wValidation(filepath)
    FilteringDF(df)
  })
  observeEvent(input$BtnConfirmFilter, {
    set_filtering_DF(FilteringDF()[input$DTPeaksFilter_rows_all,])
    FilteringDF(NULL)
    removeModal()
    # showNotification("Confirm filter.", type = "message")
  })
}
