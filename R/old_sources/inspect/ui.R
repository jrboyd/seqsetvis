source("ui_setup.R")

# ui.R definition
ui <- fluidPage(
  #Activate shinyjs
  useShinyjs(), 
  #Exgend shinjs to support tabsetsPanel operations
  extendShinyjs(text = tab_jscode),
  
  inlineCSS(tab_css),
  # extendShinyjs(text = selected_jscode),
  # HTML(selected_html ),
  # Set theme
  theme = shinytheme("spacelab"),
  # radioButtons("ToggleHelp", "Need Help?", choices = c("no", "yes", "YES!!!"), inline = T),
  fluidRow(
    h3("Example Datasets:"),
    actionButton(inputId = "ExampleMCF7_bza", label = "MCF7_bza"),
    actionButton(inputId = "ExampleKasumi", label = "Kasumi"),
    actionButton(inputId = "ExampleBivalency", label = "Bivalency"),
    actionButton(inputId = "ExampleRunxAndCTCF", label = "Runx & CTCF"),
    actionButton(inputId = "ExampleMCF10A", label = "MCF10A with Runx1"),
    actionButton(inputId = "ExampleU937", label = "U937: RUNX1 and IKZF1")
  ),
  br(),
  tabsetPanel(id = "navbar",
              tabPanel(title = "1) setup", 
                       value = "tab1",
                       # value = "1", 
                       ui_tab_setup()),
              tabPanel(title = "2) process", 
                       value = "tab2",
                       # value = "2", 
                       ui_tab_process()),
              tabPanel(title = "3) analyze", 
                       value = "tab3",
                       # value = "3", 
                       ui_tab_analyze())
  )#,
  # actionButton("stop", "Stop!")
  
)