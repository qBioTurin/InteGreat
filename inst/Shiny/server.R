#shiny.maxRequestSize=1000*1024^2
#shiny.launch.browser = .rs.invokeShinyWindowExternal

Sys.setenv("DATAVERSE_SERVER" = "dataverse.harvard.edu")
APIkey_path = system.file("Data",".APIkey", package = "ORCA")

source(system.file("Shiny","AuxFunctions.R", package = "ORCA"))
# source("./inst/Shiny/AuxFunctions.R")

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  DataAnalysisModule <- reactiveValues(wbResult = NULL,
                                       wbquantResult = NULL,
                                       endocResult = NULL,
                                       elisaResult = NULL,
                                       pcrResult = NULL,
                                       cytotoxResult = NULL,
                                       videoResult = NULL)
  
  DataIntegrationModule <- reactiveValues(dataLoaded = NULL,
                                          data = NULL,
                                          wbTabs = NULL, 
                                          pcrTabs = NULL,
                                          cytotoxTabs= NULL,
                                          endocTabs=NULL,
                                          otherTabs = NULL,
                                          otherTabsMean = NULL)
  
  MapAnalysisNames =c("WB", "WB comparison", "Endocytosis", "ELISA", "RT-qPCR", "Cytotoxicity") 
  names(MapAnalysisNames) =c("wbResult", "wbquantResult", "endocResult", "elisaResult", "pcrResult", "cytotoxResult") 
  
  
  #### Video analysis ####
  
  videoResult = reactiveValues(Initdata = NULL,
                               data = NULL)
  videoResult0 = list(Initdata = NULL,
                      data = NULL)
  
  # save everytime there is a change in the results
  VIDEOresultListen <- reactive({
    reactiveValuesToList(videoResult)
  })
  
  observeEvent(input$LoadVIDEO_Button,{
    
    if( !is.null(videoResult$Initdata) )
    { ### alert!!! if it is already present! 
      showModal(modalDialog(
        title = "Important message",
        "Do you want to update the Video already present?",
        easyClose = TRUE,
        footer= tagList(actionButton("confirmUploadVideo", "Update"),
                        modalButton("Cancel")
        )
      ))
    }
    
    output$LoadingError_Video <- renderText({
      validate(
        need(!is.null(input$PCRImport) && file.exists(input$PCRImport$datapath) ,
             "Please select a Video file!!" )
      )
      
      mess = readfile(
        filename = input$VideoImport$datapath,
        type = "Video"
      )
      
      validate(
        need(!setequal(names(mess),c("message","call")) ,
             mess[["message"]] )
      )
      
      videoResult$Initdata = mess
      
      "The Video has been uploaded  with success"
    })
  })

  observe({
    DataAnalysisModule$videoResult = reactiveValuesToList(videoResult)
  })
  
  #### END Video analysis ####
  
  ### Loading files ####
  UploadDataAnalysisModuleAllFalse  = reactiveValues(FlagALL = F,
                                                     FlagUpdate = F,
                                                     FlagWB = F,
                                                     FlagPRCC = F,
                                                     FlagELISA = F,
                                                     FlagCYTOTOX = F,
                                                     FlagENDOC = F)
  UploadDataAnalysisModule = reactiveValues(FlagALL = F,
                                            FlagUpdate = F,
                                            FlagWB = F,
                                            FlagPRCC = F,
                                            FlagELISA = F,
                                            FlagCYTOTOX = F,
                                            FlagENDOC = F)
  
  
  # upload in the statistic module
  observeEvent(input$loadStatAnalysis_file_Button,{
    output$loadStatAnalysis_Error <- renderText({
      validate(
        need(!is.null(input$loadStatAnalysis_file) && all(file.exists(input$loadStatAnalysis_file$datapath)) ,
             "Please select one RDs file generated throught the Data Analysis module." )
      )
      
      datapaths = input$loadStatAnalysis_file$datapath
      for(dpath in 1:length(datapaths)){
        mess = readRDS(datapaths[dpath])
        
        validate(
          need(all(names(mess) %in% names(DataAnalysisModule)) ||
                 all(names(mess) %in% names(elisaResult)) ||
                 all(names(mess) %in% names(wbquantResult)) || 
                 all(names(mess) %in% names(pcrResult)) ||
                 all(names(mess) %in% names(cytotoxResult)) ||
                 all(names(mess) %in% names(endocResult)) ,
               paste(mess[["message"]],"\n The file must be RDs saved throught the Data Analysis module." ))
        )
        
        DataStatisticModule$Flag = T
        
        if( all(names(mess) %in% names(wbquantResult)) || all(names(mess) %in% names(DataAnalysisModule)) ){
          DataStatisticModule$WB[[dpath]] <- mess$AdjRelDensitiy %>% mutate(DataSet = dpath)
        }else if( all(names(mess) %in% names(pcrResult)) || all(names(mess) %in% names(DataAnalysisModule))){
          DataAnalysisModule$PRCC[[dpath]]  <- mess
        }else if(all(names(mess) %in% names(endocResult)) || all(names(mess) %in% names(DataAnalysisModule))){
          DataAnalysisModule$ENDOC[[dpath]]  <- mess
        }else if(all(names(mess) %in% names(elisaResult)) || all(names(mess) %in% names(DataAnalysisModule))){
          DataAnalysisModule$ELISA[[dpath]]  <- mess
        }else if(all(names(mess) %in% names(cytotoxResult)) || all(names(mess) %in% names(DataAnalysisModule))){
          DataAnalysisModule$CYTOTOX[[dpath]]  <- mess
        }
        
      }
      
      "The RDs files have been uploaded  with success."
      
    })
  })
  
  # general upload in the app
  observeEvent(input$loadAnalysis_Button,{
    # output$loadAnalysis_Error <- renderText({
    #   validate(
    #     need(!is.null(input$loadAnalysis_file) && all(file.exists(input$loadAnalysis_file$datapath)) ,
    #          "Please select one RDs file generated throught the Data Analysis module." )
    #   )
    #   
    #   mess = readRDS(input$loadAnalysis_file$datapath)
    #   
    #   messNames = names(mess)
    #   if("Flags"%in% messNames) messNames = messNames[ messNames != "Flags"]
    #   
    #   validate(
    #     need(all(messNames %in% names(DataAnalysisModule)) ||
    #            all(messNames %in% names(elisaResult)) ||
    #            all(messNames %in% names(wbResult)) || 
    #            all(messNames %in% names(pcrResult)) ||
    #            all(messNames %in% names(cytotoxResult)) ||
    #            all(messNames %in% names(endocResult)) ,
    #          paste(mess[["message"]],"\n The file must be RDs saved throught the Data Analysis module." ))
    #   )
    #   
    #   if(all(messNames %in% names(DataAnalysisModule)) ){
    #     DataAnalysisModule <- mess
    #     UploadDataAnalysisModule$FlagALL = T
    #   }else if( all(messNames %in% names(wbResult)) ){
    #     DataAnalysisModule$wbResult <- mess
    #     UploadDataAnalysisModule$FlagWB = T
    #   }else if( all(messNames %in% names(pcrResult)) ){
    #     DataAnalysisModule$pcrResult <- mess
    #     UploadDataAnalysisModule$FlagPRCC = T
    #   }else if(all(messNames %in% names(endocResult)) ){
    #     DataAnalysisModule$endocResult <- mess
    #     UploadDataAnalysisModule$FlagENDOC = T
    #   }else if(all(messNames %in% names(elisaResult)) ){
    #     DataAnalysisModule$elisaResult <- mess
    #     UploadDataAnalysisModule$FlagELISA = T
    #   }else if(all(messNames %in% names(cytotoxResult)) ){
    #     DataAnalysisModule$cytotoxResult <- mess
    #     UploadDataAnalysisModule$FlagCYTOTOX = T
    #   }
    #   
    #   UploadDataAnalysisModule$FlagUpdate = T
    #   
    #   "The RDs file has been uploaded  with success."
    #   
    # })
  })
  
}
