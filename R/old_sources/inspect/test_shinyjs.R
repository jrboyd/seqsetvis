shinyApp(
  ui <- fluidPage(
    useShinyjs(), 
    actionButton("BtnActivateTmp", label = "Activate"),
    actionButton("BtnDisableTmp", label = "Disable"),
    actionButton("BtnToggleTmp", label = "Toggle"),
    actionButton("BtnTmp", label = "Tmp")
  ),
  server = function(input, output){
    observeEvent(input$BtnActivateTmp, {
      showNotification("enable")
      shinyjs::enable("BtnTmp")
    })
    
    observeEvent(input$BtnDisableTmp, {
      showNotification("disable")
      shinyjs::disable("BtnTmp")
    })
    
    observeEvent(input$BtnToggleTmp, {
      showNotification('toggle')
      shinyjs::toggleState("BtnTmp")
    })
    
    observeEvent(input$BtnTmp, {
      showNotification("TMP PRESSED OMG", type = "message")
    })
    
  }
)
