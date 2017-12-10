
renameModal <- function(starting_name, failed = FALSE, title = "Rename", size = "s") {
  modalDialog(
    # h3(span("Unprocessed", style = "color: red;"), 'gene reference'),
    if (failed)
      div(tags$b("Hey man, this didn't work! Not sure why.", style = "color: red;")),
    textInput("TxtRenameBigwig", label = "bigwig name", value = starting_name),
    footer = tagList(#fluidRow(
      actionButton("BtnCancelRename", "Cancel"),
      actionButton("BtnConfirmRename", "Confirm")
    ),
    size = size,
    title = title
  )
}

server_renameModal = function(input, output, session, set_name){
  observeEvent(input$BtnCancelRename, {
    removeModal()
  })
  observeEvent(input$BtnConfirmRename, {
    set_name(input$TxtRenameBigwig)
    removeModal()
  })
}