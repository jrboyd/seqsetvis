ui_tab_setup = function(){
  tagList( 
    h3("Pick a bed file (ideally after annotating with peak_annotatR)"),
    tags$hr(),
    shinyFilesButton(id = "FilesLoadSet", label = "Find Files on Server", title = "Find Peaks to Annotate", multiple = F),
    fileInput(inputId = "UploadLoadSet", label = "Browse Local Files"),
    tags$hr(),
    fixedRow(
      column(width = 1,
             br(),
             br(),
             actionButton("BtnFilterSet", "Filter"),
             br(),
             br(),
             actionButton("BtnAnnotateSet", "Annotate"),
             br(),
             br(),
             shinySaveButton(id = "FilesSaveSet", label = "Save", title = "Save annotated/filtered bed")
      ),
      column(width = 10,
             DT::dataTableOutput("SetPreview")
      )
    ),
    tags$hr(),
    h3("Pick bigwigs to visualize at bed intervals"),
    tags$hr(),
    shinyFilesButton(id = "FilesLoadBigwig", label = "Find bigwig on Server", title = "Find Peaks to Annotate", multiple = F),
    textInput("TxtAddBigWig", label = "Bigwig name"),
    actionButton(inputId = "BtnAddBigwig", label = "Add Bigwig"),
    tags$hr(),
    fixedRow(
      column(width = 1,
             br(),
             br(),
             actionButton(inputId = "BtnRemoveBigWig", "Remove"),
             br(),
             br(),
             actionButton(inputId = "BtnRenameBigWig", "Rename")
      ),
      column(width = 10,
             DT::dataTableOutput("AddedBigWigs")
      )
    ),
    
    
    tags$hr(),
    actionButton("BtnFinishSetup", label = "Finish Setup")
  )
}
