ui_tab_process = function(){
  tagList(
    h3("Bed Description"),
    DT::dataTableOutput("BedSummary"),
    textOutput("BedLength"),
    h3("Bigwigs Description"),
    DT::dataTableOutput("BigWigSummary"),
    actionButton("BtnFinishProcess", label = "Process Profiles")
    
  )
}

