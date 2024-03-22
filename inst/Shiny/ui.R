library(shinydashboard)
library(shiny)
library(shinyalert)
library(zoo)
library(knitr)
library(ggplot2)
library(shinythemes)
library(OpenImageR)
library(dplyr)
library(shinyWidgets)
library(DT)
library(openxlsx)
library(patchwork)

ui <- dashboardPage(
  dashboardHeader(title = "ORCA",
                  tags$li(a(onclick = "window.open('https://github.com/qBioTurin/ORCA')",
                            href = NULL,
                            icon("github"),
                            title = "GitHub",
                            style = "cursor: pointer;"),
                          class = "dropdown"),
                  tags$li(class = "dropdown",
                          tags$style("#mydropdown {cursor: pointer; background-color: #852fb4; color: white; border: 1px solid #852fb4; height: 50px; left: 0;}"),
                          dropdownButton(inputId = "mydropdown", label = "", right = TRUE, circle = FALSE, icon = icon("download"),
                                         tags$li(shinyjs::useShinyjs(),
                                                 downloadButton(outputId = 'downloadReport', label = "Report", title = "Report generation", style = "cursor: pointer; width: 98%; text-align: center; vertical-align: middle; border: 1px solid #9809AF;", class = "dlButton")),
                                         tags$li(shinyjs::useShinyjs(),
                                                 downloadButton(outputId = 'downloadRDSwholeAnalysis', label = "RDs", title = "RDS storing the whole analysis", style = "cursor: pointer; width: 98%; text-align: center; vertical-align: middle; border: 1px solid #9809AF;", class = "dlButton"))
                          ),
                          block = TRUE)
  ),
  dashboardSidebar(
    sidebarMenu(id = "SideTabs",
                menuItem('Home', tabName = 'Home', icon = icon('home')),
                menuItem("Data Analysis", tabName = 'DataAnaslysis', icon = icon('chart-line'),
                         menuItem('Western Blot analysis', tabName = 'wb',
                                  menuSubItem("Upload Image", tabName = "uploadIm"),
                                  menuSubItem("Protein Bands", tabName = "plane"),
                                  menuSubItem("Profile Plots", tabName = "grey"),
                                  menuSubItem("Quantification", tabName = "quantification")),
                         menuItem('RT-qPCR analysis', tabName = 'pcr',
                                  menuSubItem("Upload data", tabName = "uploadPCR"),
                                  menuSubItem("Quantification", tabName = "tablesPCR")),
                         menuItem('ELISA analysis', tabName = 'elisa',
                                  menuSubItem("Upload data", tabName = "uploadELISA"),
                                  menuSubItem("Quantification", tabName = "tablesELISA")),
                         menuItem('Endocytosis assay', tabName = 'endoc',
                                  menuSubItem("Upload data", tabName = "uploadENDOC"),
                                  menuSubItem("Quantification", tabName = "tablesENDOC")),
                         menuItem('Cytotoxicity assay', tabName = 'cytotox',
                                  menuSubItem("Upload data", tabName = "uploadCYTOTOX"),
                                  menuSubItem("Quantification", tabName = "tablesCYTOTOX"))
                ),
                
                menuItem('Statistical analysis', tabName = 'StatAnalysis_tab', icon = icon('magnifying-glass-chart')),
                menuItem('Dataverse', tabName = 'Dataverse_tab', icon = icon('eye')),
                menuItem('Load analysis', tabName = 'LoadAnalysis', icon = icon('upload'))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "Home",
              h1("ORCA: Omni Reproducible Cell Analysis"),
              h1(" "),
              fluidRow(
                column(width = 6,
                       p(img(src = "ORCAlogo.png", height = "30%", width = "30%"), align = "center")
                ),
                column(width = 6,
                       h2(em("A cellular biologist’s toolbox for data analysis.")),
                       h4("ORCA provides an exhaustive platform where scientists can analyze raw:"),
                       tags$ol(
                         tags$li(h4(strong("Western Blot"), " (WB),")),
                         tags$li(h4(strong("Reverse Transcription-quantitative PCR"), " (RT-qPCR),")),
                         tags$li(h4(strong("Enzyme-Linked ImmunoSorbent Assay"), " (ELISA),")),
                         tags$li(h4(strong("Endocytosis"), " and,")),
                         tags$li(h4(strong("Cytotoxicity experiments"), "."))
                       )
                )
              ),
              p(img(src = "Logo_QBio.png", height = "15%", width = "15%", style = "margin:100px 0px"), align = "center")
      ),
      tabItem(tabName = "LoadAnalysis",
              h2("Load analysis"),
              box(
                width = 12,
                fluidRow(
                  column(
                    10,
                    fileInput(inputId = "loadAnalysis_file", label = "", placeholder = "Select the RDs files storing ORCA analyses", width = "80%", multiple = TRUE)
                  ),
                  column(
                    1,
                    actionButton(label = "Load", style = "margin-top: 20px;", icon = shiny::icon("upload"), inputId = "loadAnalysis_Button")
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
      tabItem(tabName = "uploadIm",
              h2("Upload Image"),
              fluidRow(
                column(9,
                       fileInput(
                         inputId = "imImport",
                         label = "Select a tif file",
                         placeholder = "Select a tif file",
                         width = "90%"
                       )
                ),
                column(2,
                       actionButton(label = "Load",style = "margin-top: 20px;",
                                    icon = shiny::icon("upload"),
                                    inputId = "LoadingTif"
                       )
                ),
                tags$style(type='text/css', "#LoadingTif { width:100%; margin-top: 20px;}")
              ),
              fluidRow(
                column(9, offset = 1,
                       tags$h5("In section «Upload Image», you have to upload the original file.tif and then you are directly redirected to «Protein Band» section.",
                               style = "margin-top: 20px; text: center") # Stile aggiunto per un po' di margine
                )
              )
      ),
      tabItem(tabName = "plane",
              h2("Select Protein Bands"),
              fluidRow(
                uiOutput("TiffBox")
              ),
              fluidRow(
                box( width = 6,
                     title = tagList(shiny::icon("gear", verify_fa = FALSE), 
                                     "Protein Band Selection Coordinates"),
                     h5(em("Click on the SampleName column to associate the lane with a sample name. Let us note that equal sample names are not allowed.")),
                     DTOutput("PlanesStructureTable")
                ),
                box( width = 6,
                     actionButton(inputId = "panelSelect_button", label = "Select Protein Bands"),
                     actionButton(inputId = "ResetPan", label = 'Reset Protein Bands'),
                     actionButton(inputId = "GenLanes", label = 'Generate Plots'),
                     verbatimTextOutput("rectCoordOutput")
                )
              )
      ),
      tabItem(tabName = "grey",
              h2("Signal Profiles"),
              fluidRow(
                column(width = 6, offset = 0.5,
                       selectInput(inputId = "LaneChoice",
                                   label = "Choose a lane:",
                                   choices = c("")
                       ) 
                )
              ),
              box(width = 12,
                  plotOutput("DataPlot")
              ),
              fluidRow(
                box(width = 6,
                    tabsetPanel(id = "tabs",
                                tabPanel("Vertical cut", value= "V",textOutput("V"),
                                         sliderInput(inputId = "truncV", label = h4("Vertical truncation:"),
                                                     min = 0, max = 0, value = c(0,0),step = 1
                                         ),
                                         actionButton( 
                                           label = "Cut", inputId = "actionButton_TruncV",
                                           icon = icon("cut")
                                         )
                                ),
                                tabPanel("Horizontal cut", value= "H",textOutput("H"),
                                         sliderInput(inputId = "truncH", label = h4("Horizontal truncation:"),
                                                     min = 0, max = 0, value = 0),
                                         actionButton( label = "Cut", inputId = "actionButton_TruncH",
                                                       icon = icon("cut") 
                                         )
                                )
                    )
                ),
                box(width=6,
                    DTOutput('AUC'),
                    fluidRow(
                      column(width = 3,
                             actionButton(label = "Reset all", 
                                          inputId = "actionButton_ResetPlanes")
                      ),
                      column(width = 4,
                             downloadButton(label = "Download Analysis & Excel", 
                                            outputId = "downloadWBAnalysis",
                                            icon = icon("download"))
                      ),
                      column(width = 3,
                             actionButton(inputId = "NextWBQuantif",
                                          label = 'Proceed to Quantification',
                                          align = "right",
                                          icon = shiny::icon("forward"))
                      )
                    )
                )
              )
      ),
      tabItem(tabName = "quantification",
              h2("WB quantification"),
              fluidRow(
                box( width = 6,
                     title = tagList(shiny::icon("gear", verify_fa = FALSE), "Set the WB analysis as normalizer"),
                     h5("To keep the necessary samples (rows) just clicking on it."),
                     column(9,
                            fileInput(
                              inputId = "NormWBImport",
                              label = "",
                              placeholder = "Select an WB RDs file generated through the Profile Plots step",
                              width = "90%"
                            )
                     ),
                     column(3,
                            actionButton(
                              label = "Load",
                              icon = shiny::icon("upload"),
                              width = "100%",
                              inputId = "actionB_loadingNormWB"
                            )
                     ),
                     tags$style(type='text/css',
                                "#actionB_loadingNormWB { width:100%; margin-top: 20px;}"
                     ),
                     fluidRow(
                       column(
                         width = 10,
                         offset = 1,
                         verbatimTextOutput("LoadingErrorNormWB")
                       )
                     ),
                     fluidRow(
                       column(
                         width = 10,
                         offset = 1,
                         DTOutput('AUC_WBnorm')
                       )
                     )
                ),
                box( width = 6,
                     title = tagList(shiny::icon("gear", verify_fa = FALSE),
                                     "Set the WB analysis to normalize"),
                     h5("To keep the necessary samples (rows) just clicking on it."),
                     column(9,
                            fileInput(inputId = "WBImport",
                                      label = "",
                                      placeholder = "Select an WB RDs file generated through the Profile Plots step",
                                      width = "90%"
                            )
                     ),
                     column(3,
                            actionButton( label = "Load",
                                          icon = shiny::icon("upload"),
                                          width = "100%",
                                          inputId = "actionB_loadingWB"
                            )
                     ),
                     tags$style(type='text/css',
                                "#actionB_loadingWB { width:100%; margin-top: 20px;}"
                     ),
                     fluidRow(
                       column(
                         width = 10,
                         offset = 1,
                         verbatimTextOutput("LoadingErrorWB")
                       )
                     ),
                     fluidRow(
                       column(
                         width = 10,
                         offset = 1,
                         DTOutput('AUC_WB')
                       )
                     )
                )
              ),
              box( width = 12,
                   title = tagList(shiny::icon("gear", verify_fa = FALSE), "WB quantification"),
                   h3("Relative Density"),
                   fluidRow(
                     column(
                       width = 8,
                       selectInput("IdLaneNorm_RelDens",
                                   label = "Select the Lane to use for the relative normalization",
                                   choices = "Nothing selected",selected = "Nothing selected")
                     ),
                     column(
                       width = 12,
                       DTOutput('AUC_RelDens')
                     )
                   ),
                   h3("Adjusted Relative Density"),
                   fluidRow(
                     column(
                       width = 12,
                       DTOutput('AUC_AdjRelDens')
                     )
                   ),
                   plotOutput("plot_AdjRelDens"),
                   fluidRow(
                     column(width = 1, offset = 7,
                            downloadButton( label = "Download the analysis", 
                                            outputId = "downloadButton_WBquant",
                                            icon = icon("download") 
                            )
                     ),
                     column(width = 1,offset = 9,
                            downloadButton( label = "Download xlsx", 
                                            outputId = "downloadButtonExcel_WBquant",
                                            icon = icon("download") )
                     )
                   )
              )
      )
    )
  )
)

