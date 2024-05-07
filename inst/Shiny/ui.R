library(shinydashboard)
library(shiny)
library(zoo)
library(knitr)
library(ggplot2) 
library(shinythemes) 
library(OpenImageR)
library(dplyr)
library("shinyWidgets")
library(DT)
library(openxlsx)
library(patchwork)
library(shinyjs)
library(shinybusy)
ui <- dashboardPage(
  
  #theme = shinytheme("paper"),
  dashboardHeader(title = "ORCA",
                  tags$li(a(onclick = "onclick =window.open('https://github.com/qBioTurin/ORCA')",
                            href = NULL,
                            icon("github"),
                            title = "GitHub",
                            style = "cursor: pointer;"),
                          class = "dropdown"),
                  ## DOWNLOAD
                  tags$li(class = "dropdown",
                          tags$style("#mydropdown{cursor: pointer; background-color: #852fb4 ;
                              color: white;
                              border: 1px solid #852fb4; height: 50px;left:0;}"),
                          dropdownButton(inputId = "mydropdown",label = "",
                                         right = T,
                                         circle = FALSE,
                                         icon = icon("download"),
                                         tags$li(shinyjs::useShinyjs(),
                                                 downloadButton(outputId = 'downloadReport',
                                                                href = "report.html",
                                                                download = "report.html",
                                                                label = "Report",title = "Report generation",
                                                                style = "cursor: pointer; width: 98%;
                                      text-align: center; vertical-align: middle;
                                      border: 1px solid #9809AF;",
                                                                class="dlButton")
                                         ),
                                         tags$li(shinyjs::useShinyjs(),
                                                 downloadButton(outputId = 'downloadRDSwholeAnalysis',
                                                                label = "RDs",title = "RDS storing the whole analysis",
                                                                style = "cursor: pointer; width: 98%;
                                      text-align: center; vertical-align: middle;
                                      border: 1px solid #9809AF;",
                                                                class="dlButton")
                                         )
                          ),
                          block = TRUE)
  ),
  dashboardSidebar(
    sidebarMenu(id = "SideTabs",
                menuItem('Home',
                         tabName = 'Home',
                         icon = icon('home')
                ),
                menuItem("Data Analysis",
                         tabName = 'DataAnalysis',
                         icon = icon('chart-line'),
                         menuItem('Border Identification',
                                  tabName = 'border',
                                  menuSubItem("Upload Video", tabName = "uploadVideo")
                         )
                ),
                menuItem('Statistical analysis',
                         tabName = 'StatAnalysis_tab',
                         icon = icon('magnifying-glass-chart')
                ),
                menuItem('Dataverse',
                         tabName = 'Dataverse_tab',
                         icon = icon('eye')
                ),
                menuItem('Load analysis',
                         tabName = 'LoadAnalysis',
                         icon = icon('upload')
                )
    )
  ),
  dashboardBody(
    tabItems(
      ## HOME ####
      tabItem(
        tabName = "Home",
        h1("ORCA: Omni Reproducible Cell Analysis"),
        h1("  "),
        fluidRow(
          column(width = 6,
                 p(img(src = "ORCAlogo.png",
                       #p(img(src = system.file("Shiny/www/images","ORCAlogo.png", package = "ORCA"),
                       height="30%", width="30%"), align = "center"),
          ),
          column(width = 6,
                 h2(em("A cellular biologist’s toolbox for data analysis.")),
                 h4("ORCA  provides an exhaustive platform where scientists can analyze raw:"),
                 tags$ol(
                   tags$li(
                     h4(strong("Western Blot")," (WB), ")
                   ),
                   tags$li(
                     h4(strong("Reverse Transcription-quantitative PCR ")," (RT-qPCR),")
                   ),
                   tags$li(
                     h4(strong("Enzyme-Linked ImmunoSorbent Assay ")," (ELISA),")
                   ),
                   tags$li(
                     h4(strong("Endocytosis")," and,")
                   ),
                   tags$li(
                     h4(strong("Cytotoxicity experiments"),".")
                   )
                 )
          )
        ),
        p(img(src = "Logo_QBio.png",
              height="15%", width="15%",style = "margin:100px 0px"), align = "center")
        # h4("ORCA consists of two modules:"),
        # tags$ol(
        #   tags$li(
        #     h4("The ", strong("Data Analysis module"), "includes tools specifically developed or adapted for the elaboration
        #    of raw Western Blot (WB), Reverse Transcription-quantitative PCR (RT-qPCR) and Enzyme-Linked ImmunoSorbent Assay (ELISA) experiments.")
        #   ),
        #   tags$li(
        #     h4("The ",strong("Model Integration module"),"supports scientists in the process of integration of lab data resulting from any type of experiment into a computational model.")
        #   )
        # ),
        # h4(em("Check", a("here", href="https://www.google.com/")," for a brief video presentation of the ORCA framework, or ",
        #       a("here", href="https://www.google.com/"),"to download the user guide.")),
      ),
      ###### BEGIN LOAD ANALYSIS ####
      tabItem(tabName = "LoadAnalysis",
              h2("Load analysis"),
              box(
                width = 12,
                fluidRow(
                  column(
                    10,
                    fileInput(
                      inputId = "loadAnalysis_file",
                      label = "",
                      placeholder = "Select the RDs files storing ORCA analyses",
                      width = "80%", 
                      multiple = TRUE)
                  ),
                  column(
                    1,
                    actionButton( label = "Load",style = "margin-top: 20px;",
                                  icon = shiny::icon("upload"),
                                  inputId = "loadAnalysis_Button" )
                  )
                ),
                fluidRow(
                  column(
                    width = 10,
                    offset = 1,
                    verbatimTextOutput("loadAnalysis_Error")
                  )
                )
              )
      ),
      
      ###### BEGIN DATA ANALYSIS ####
      ## BEGIN data integration: Video  #######
      # First tab content
      tabItem(tabName = "uploadVideo",
              shinyjs::useShinyjs(),
              h2("Load video"),
              fluidRow(
                column(10,
                       fileInput(inputId = "uploaded_video",
                                 label = "",
                                 placeholder = "Choose .lif file",accept = ".lif",
                                 width = "80%"),
                      
                ),
                column(1,
                       style = "margin-top: 20px;",
                       use_busy_spinner(spin = "fading-circle"),
                       actionButton("uploadVideo_button", "Upload Video", icon = shiny::icon("upload")),
                )
              ),
              
              fluidRow(
                uiOutput("imageGallery")
                
              ),
              
              fluidRow(
                column(
                  width = 10,offset = 1,
                  verbatimTextOutput("LoadingError_Video")
                )
              )
      ),
      # Second tab content
      tabItem(tabName = "cellsVisual",
              h2("Cells Identification"),
              fluidRow(
              )
      )
      ## END data analysis: Video 
      
      # Here ends the ui body
    )
    
  )
  
)
