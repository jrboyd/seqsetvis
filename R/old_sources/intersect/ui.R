# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(shinythemes)
source("chooser.R")

radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL){
  
  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}


shinyUI(fluidPage(
  theme = shinytheme("spacelab"),
  # Application title
  titlePanel("peak intersectR"),
  # radioButtons(inputId = "StrategyRadio", label = "Strategy Type", choices = c("serial (pie chart friendly)" = "serial", "flat (venn diagram friendly)" = "flat")),
  br("step1 - file selection"),
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(width = 4,
                 fixedRow(
                   column(width = 12,
                          fileInput(inputId = "BtnUploadPeakfile", label = "Browse Local Files"),
                          h5("or"),
                          shinyFilesButton(id = "FilesLoadData", label = "Find Files on Server", title = "Find Peaks to Annotate", multiple = F)
                   )
                 ),
                 br(),
                 fixedRow(
                   column(width = 6,
                          textInput(inputId = "TxtFileName", value = "", placeholder = "file description", label = "Descriptive Name")
                   )
                 ),
                 fixedRow(
                   column(width = 3, 
                          actionButton(inputId = "BtnAddFile", label = "Add")
                   ),
                   column(width = 3,
                          actionButton(inputId = "BtnCancelFile", label = "Cancel")
                   )
                 )#,
    ),
    mainPanel(
      tags$h3("Preview of Loaded File : ", tags$span("this can be filtered and sorted in Step 2.", style="color:red")),
      DT::dataTableOutput("DTPeaksPreview")
    )
  ),
  br("step2 - organize"),
  fixedRow(
    column(width = 4,
           uiOutput("SetChooser")
    ),
    column(width = 4,
           actionButton("BtnDeleteSet", label = "Delete"),
           actionButton("BtnRenameSet", label = "Rename"),
           actionButton("BtnFilterSet", label = "Filter")#,
           # bsModal(id = "ModalFilter", 
           #         title = "Filter", 
           #         trigger = "BtnFilterSet", 
           #         DT::dataTableOutput(outputId = "DTPeakSFilter"))
           
    )
  ),
  br("step3 - analyze"),
  sidebarLayout(
    sidebarPanel(width = 4,
                 fixedRow(
                   column(width = 4,
                          radioButtons(inputId = "StrategyRadio", label = "Strategy Type", 
                                       choices = c("flat (venn)" = "flat", "serial (pie)" = "serial"))
                   ),
                   column(width = 4,
                          radioButtons(inputId = "RadioPlotType", label = "Plot Type", choices = c("bars", "pie", "venn", "euler", "heatmap"))
                   ),
                   column(width = 4,
                          sliderInput(inputId = "SliderMergeExtension", label = "Extension", min = 0, max = 5*10^3, value = 0, round = 2, step = 100),
                          uiOutput(outputId = "NumericMergeExtensionOut")
                          
                   )
                 ),
                 radioTooltip(id = "StrategyRadio", choice = "serial", title = "pie chart appropriate. annotate each interval in first file to one and only one category (1 category per additional file).", placement = "right", trigger = "hover"),
                 radioTooltip(id = "StrategyRadio", choice = "flat", title = "venn diagram appropriate. combine and flatten all intervals then assign each interval to as many categories as apply (1 category per file)", placement = "right", trigger = "hover"),
                 actionButton("BtnQuickFlat", label = "Load Example"),
                 conditionalPanel(condition = "input.StrategyRadio == 'serial'",
                                  tags$div(
                                    HTML(paste("", tags$span(style="color:red", "Strategy is serial: First file will be annotated serially by all others.  Earlier files take precedence over later."), sep = ""))
                                  )),
                 # actionButton("BtnAnalyze", label = "Analyze!"),
                 downloadButton(outputId = "DownloadResults", label = "Download Results"),
                 downloadButton(outputId = "DownloadPlot", label = "Download Plot"),
                 # hidden(actionButton("BtnSaveResults", label = "Download Results")),
                 shinySaveButton(id = "FilesSaveResults", label = "Save Results (for inspectR)", title = "Save Results")
    ),
    mainPanel( 
      plotOutput("AnalysisPlot", width = "400", height = "400")
    )
    
  ),
  verbatimTextOutput("DebugTxt", placeholder = T)
  # todoListInput("ListOrganize", leftLabel = "My cool stuff", rightLabel = "b", leftChoices = 1:3, rightChoices = 0)
)
)

