#shiny.maxRequestSize=1000*1024^2
#shiny.launch.browser = .rs.invokeShinyWindowExternal

Sys.setenv("DATAVERSE_SERVER" = "dataverse.harvard.edu")
APIkey_path = system.file("Data",".APIkey", package = "ORCA")

#source(system.file("Shiny","AuxFunctions.R", package = "ORCA"))
source("./inst/Shiny/AuxFunctions.R")

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  alert <- reactiveValues(alertContext = "")
  
  DataAnalysisModule <- reactiveValues(wbResult = NULL,
                                       wbquantResult = NULL,
                                       endocResult = NULL,
                                       elisaResult = NULL,
                                       pcrResult = NULL,
                                       cytotoxResult = NULL,
                                       facsResult = NULL)
  
  DataIntegrationModule <- reactiveValues(dataLoaded = NULL,
                                          data = NULL,
                                          wbTabs = NULL, 
                                          pcrTabs = NULL,
                                          cytotoxTabs= NULL,
                                          endocTabs=NULL,
                                          otherTabs = NULL,
                                          otherTabsMean = NULL)
  
  MapAnalysisNames =c("WB", "WB comparison", "Endocytosis", "ELISA", "RT-qPCR", "Cytotoxicity", "FACS") 
  names(MapAnalysisNames) =c("wbResult", "wbquantResult", "endocResult", "elisaResult", "pcrResult", "cytotoxResult", "facsResult") 
  
  ### WB analysis ####
  
  PanelData = data.frame(SampleName = character(),
                         xmin = numeric(), ymin = numeric(), 
                         xmax = numeric(), ymax = numeric())
  
  
  wbResult <- reactiveValues(
    Normalizer = NULL,
    Im = NULL,
    Planes = NULL,
    TruncatedPanelsValue = NULL,
    PanelsValue = NULL,
    Plots = NULL,
    TruncatedPlots = NULL,
    pl = NULL,
    AUCdf=data.frame(SampleName = "-", Truncation = "-", AUC = "-"  ))
  
  wbResult0 <- list(   Normalizer = NULL,
                       Im = NULL,
                       Planes = NULL,
                       TruncatedPanelsValue = NULL,
                       PanelsValue = NULL,
                       Plots = NULL,
                       TruncatedPlots = NULL,
                       pl = NULL,
                       AUCdf=data.frame(SampleName = "-", Truncation = "-", AUC = "-" )
  )
  
  
  WBresultList <- reactive({
    reactiveValuesToList(wbResult)
  })
  observeEvent(WBresultList(), {
    DataAnalysisModule$wbResult = reactiveValuesToList(wbResult)
    DataAnalysisModule$wbResult$Flags = reactiveValuesToList(Flags)
  })
  
  Flags <- reactiveValues( ShowTif = F, 
                           LanesCut = F,
                           CutTab="V",
                           IDlane = 0)
  prev_vals <- NULL
  PanelStructures <- reactiveValues(data = PanelData )
  NumberOfPlanes <- reactiveValues(N = 0)
  PlaneSelected <- reactiveValues(First = NULL)
  

  observeEvent(input$LoadingTif, {
    alert$alertContext <- "WB-reset"
    if(!is.null(wbResult$Im) ) { 
      shinyalert(
        title = "Important message",
        text = "Do you want to update the WB data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else loadTifFile()
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "WB-reset") {  
      resetPanel("WB", 
                 Flags, 
                 PanelStructures, 
                 NumberOfPlanes, 
                 PlaneSelected, 
                 wbResult, 
                 output, PanelData)
      updateSelectInput(session, "LaneChoice", choices = "")
    
      loadTifFile()
    }
  })
  
  loadTifFile <- function() {
    alert$alertContext <- ""
    file <- !is.null(input$imImport) && file.exists(input$imImport$datapath)
    mess = readfile(filename = input$imImport$datapath, type = "tif", file)
    
    
    if(setequal(names(mess), c("message", "call"))) {
      showAlert("Error", mess[["message"]], "error", 5000)
    } else {
      Flags$ShowTif <- TRUE
      wbResult$Im = mess
      updateTabsetPanel(session, "SideTabs", selected = "plane")
      showAlert("Success", "The image was uploaded successfully!", "success", 1000)
    }
  }
  
  observe({
    if(Flags$ShowTif) {
      wbResult$Im -> ListIm 
      im = ListIm$RGB
      
      output$TiffBox <- renderUI({
        column(12,align="center",
               box(plotOutput("TifPlot2",
                              hover = "plot_hover",
                              brush = "plot_brush"),
                   width = 12,
                   height = dim(im)[1]+0.1*dim(im)[1])
        )
      })
      
      output$TifPlot2 <- renderPlot({ 
        plot(c(1,dim(im)[2]),c(1,dim(im)[1]), type='n',ann=FALSE)
        rasterImage(im,1,1,dim(im)[2],dim(im)[1])
        if (nrow(PanelStructures$data) > 0) {
          r <- PanelStructures$data
          rect(r$xmin, r$ymin, r$xmax, r$ymax, border = "red")
        }
      }, width  = dim(im)[2],height = dim(im)[1] )
      
      Flags$ShowTif = FALSE
    }
  })
  
  observeEvent(input$panelSelect_button,{
    e <- input$plot_brush
    if (!is.null(e)) {
      vals <- data.frame(xmin = round(e$xmin, 1),
                         ymin = round(e$ymin, 1),
                         xmax = round(e$xmax, 1),
                         ymax = round(e$ymax, 1))
      
      if (!identical(vals,prev_vals))  
      {
        NumberOfPlanes$N = NumberOfPlanes$N + 1
        if(NumberOfPlanes$N > 1){
          vals$ymax = prev_vals$ymax
          vals$ymin = prev_vals$ymin
          
          newH = vals$ymax - vals$ymin
          newW = vals$xmax - vals$xmin
          prevH = prev_vals$ymax - prev_vals$ymin
          prevW = prev_vals$xmax - prev_vals$xmin
          
          if( abs((newH - prevH) + (newW - prevW)) > 1 )
          {
            NumberOfPlanes$N = 1
            PanelStructures$data <- data.frame(SampleName = "1",vals)
          }else{
            PanelStructures$data <- rbind(PanelStructures$data,
                                          cbind(data.frame(
                                            SampleName = paste(nrow(PanelStructures$data)+1) ),
                                            vals))
          }
        }else{
          PanelStructures$data <- rbind(PanelStructures$data,
                                        cbind(data.frame(
                                          SampleName = paste(nrow(PanelStructures$data)+1) ),
                                          vals))
        }
        prev_vals <<- vals
      }
    } else showAlert("Error", "please select before a protein bands", "error", 5000)
  })
  
  observeEvent(input$ResetPan,{
    NumberOfPlanes$N = 1
    PanelStructures$data <- PanelData
  })
  
  output$rectCoordOutput <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
    }
    xy_range_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
             " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
    }
    
    paste0(
      "Mouse hover: ", xy_str(input$plot_hover),
      "New band coordinates: ", xy_range_str(input$plot_brush)
    )
  })
  
  output$PlanesStructureTable <- renderDT(
    {
      datatable(PanelStructures$data,
                editable = list(target = "cell", 
                                disable = list(columns = 1:4)),
                options = list(lengthChange = FALSE, autoWidth = TRUE),
                rownames= FALSE
      )
    }
  )
  
  observeEvent(input$PlanesStructureTable_cell_edit, {
    cells = input$PlanesStructureTable_cell_edit
    
    data = PanelStructures$data %>% filter(SampleName == PanelStructures$data[cells$row,"SampleName"])
    
    if(cells$value %in% data$SampleName){
      k = table(data$SampleName)[cells$value]
      cells$value = paste0(cells$value, " (",k,")")
    }
    cells$col = 1
    
    PanelStructures$data <- editData( PanelStructures$data , cells, 'PlanesStructureTable')
    wbResult$Planes = PanelStructures$data
  })
  
  observeEvent(wbResult$AUCdf,{
    output$AUC <-  renderDT({
      wbResult$AUCdf %>%
        dplyr::select(SampleName,Truncation, AUC)
    },
    selection = 'none',
    rownames= FALSE)
  })
  
  observeEvent(input$GenLanes,{
    if(NumberOfPlanes$N >1){
      Planes = PanelStructures$data
      Planes[,-1] = round(Planes[,-1])
      wbResult$Planes = Planes
      Flags$LanesCut= T
    }else{
      Flags$LanesCut= F
      showAlert("Error", "please select before a protein bands", "error", 5000)
    }
    
    if(Flags$LanesCut)
    {
      updateTabsetPanel(session, "SideTabs",
                        selected = "grey")
      
      im = wbResult$Im$WB
      PanelData = wbResult$Planes
      
      PanelsValue = do.call("rbind",
                            lapply(1:dim(PanelData)[1],
                                   function(i,im,PanelData){
                                     p = PanelData[i,]
                                     Nrow = dim(im)[1]
                                     Ncol= dim(im)[2]
                                     plane = abs(im[(Nrow-p$ymax):(Nrow-p$ymin),p$xmin:p$xmax]-1)
                                     
                                     GreyPlane = apply(plane,1,"mean")
                                     data.frame(Values = GreyPlane,
                                                ID = paste0(i,". ",p$SampleName),
                                                Y = 1:length(GreyPlane) )
                                   },
                                   im,PanelData)
      )
      
      
      pl <- ggplot(PanelsValue, aes(x =Y,y=Values)) +
        geom_line() + 
        theme_bw() +
        facet_wrap(~ID) + 
        lims(y=c(0,max(PanelsValue$Values)))
      
      wbResult$PanelsValue <- PanelsValue
      wbResult$Plots <- pl
      
      if (!is.null(PanelsValue$ID) && length(PanelsValue$ID) > 0) {
        updateSelectInput(session, "LaneChoice",
                          choices = unique(PanelsValue$ID),
                          selected = unique(PanelsValue$ID)[1])
      } else {
        updateSelectInput(session, "LaneChoice",
                          choices = c("No lanes available"),
                          selected = "No lanes available")
      }
      
      output$DataPlot <- renderPlot({pl})
      
      aucList = lapply(unique(PanelsValue$ID), function(IDlane) AUCfunction(wbResult0$AUCdf,PanelsValue,SName = IDlane) )
      wbResult$AUCdf <- do.call(rbind,aucList)
    }
  })
  
  output$downloadWBAnalysis <- downloadHandler(
    filename = function() {
      paste('WBanalysis-', Sys.Date(), '.zip', sep='')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      
      nomeRDS <- paste0("WBquant_analysis-", Sys.Date(), ".rds")
      nomeXLSX <- paste0("WBquant_analysis-", Sys.Date(), ".xlsx")
      
      tempRdsPath <- file.path(tempDir, nomeRDS)
      tempXlsxPath <- file.path(tempDir, nomeXLSX)
      results <- DataAnalysisModule$wbResult
      saveRDS(results, file = tempRdsPath)
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "WB", PanelStructures)
      
      zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
    },
  )
  
  observeEvent(input$actionButton_ResetPlanes,{
    
    wbResult$AUCdf =  wbResult$AUCdf %>% filter(Truncation == "-") 
    wbResult$TruncatedPanelsValue = wbResult0$TruncatedPanelsValue
    wbResult$TruncatedPlots = wbResult0$TruncatedPlots 
    
    output$AUC <-  renderDT({wbResult$AUCdf %>%
        dplyr::select(SampleName,Truncation, AUC)},
        selection = 'none', 
        rownames= FALSE,
    )
    
    output$DataPlot <- renderPlot({wbResult$Plots})
  })
  
  
  observeEvent(c(input$actionButton_TruncV,input$actionButton_TruncH),{
    
    if( !is.null(wbResult$PanelsValue))
    {
      Flags$IDlane -> IDlane
      if(!is.null(wbResult$TruncatedPanelsValue ))
      {
        pl <- wbResult$TruncatedPlots    
        wbResult$TruncatedPanelsValue -> PanelsValue
      }
      else{
        wbResult$PanelsValue -> PanelsValue
        pl<-wbResult$Plots
      }
      
      maxPanelsValue=max(wbResult$PanelsValue$Values)
      wbResult$AUCdf -> AUCdf
      AUCdf.new <- AUCdf[length(AUCdf$Truncation),]
      
      lastTrunc = AUCdf %>% 
        group_by(SampleName) %>%
        filter(SampleName == IDlane, row_number()==n() ) %>%
        ungroup() %>%
        dplyr::select(Truncation) 
      
      if(length(lastTrunc$Truncation) > 0 & lastTrunc$Truncation != "-")
        AUCdf.new$Truncation <- lastTrunc$Truncation
      
      AUCdf.new$SampleName <- IDlane
      
      if(Flags$CutTab=="V")
      {
        MinTrunc<-input$truncV[1]
        MaxTrunc<-input$truncV[2]
        AUCdf.new$Truncation <- paste(AUCdf.new$Truncation ,";\n X = [", MinTrunc," ; ", MaxTrunc ,"]",collapse = "")
        PanelsValue<- PanelsValue[!((PanelsValue$Y < MinTrunc | PanelsValue$Y > MaxTrunc) & PanelsValue$ID == IDlane),]
        PanelsValue$Values[PanelsValue$ID == IDlane] <- PanelsValue$Values[PanelsValue$ID == IDlane] -min(PanelsValue$Values[PanelsValue$ID == IDlane]) 
        
        updateSliderInput(session,"truncV",
                          min = min(PanelsValue$Y[PanelsValue$ID == IDlane]),
                          max = max(PanelsValue$Y[PanelsValue$ID == IDlane]),
                          value = c(min(PanelsValue$Y[PanelsValue$ID == IDlane]),
                                    max(PanelsValue$Y[PanelsValue$ID == IDlane]) ) )
      }
      else if(Flags$CutTab=="H")
      {
        TruncY<-input$truncH[1]
        PanelsValue <- PanelsValue[!(PanelsValue$Values<TruncY & PanelsValue$ID == IDlane),]
        PanelsValue$Values[PanelsValue$ID == IDlane] <- PanelsValue$Values[PanelsValue$ID == IDlane] - TruncY
        AUCdf.new$Truncation <- paste(AUCdf.new$Truncation ,";\n Y = ", TruncY)
        
        updateSliderInput(session,"truncH",
                          min = min(PanelsValue$Values[PanelsValue$ID == IDlane]),
                          max = max(PanelsValue$Values[PanelsValue$ID == IDlane]),
                          value = c(min(PanelsValue$Values[PanelsValue$ID == IDlane]),
                                    max(PanelsValue$Values[PanelsValue$ID == IDlane]) ) )
      }
      
      pl <- ggplot(PanelsValue, aes(x =Y,y=Values)) +
        geom_line() + 
        theme_bw() +
        facet_wrap(~ID)+ 
        lims(y=c(0,maxPanelsValue))
      
      wbResult$TruncatedPanelsValue <- PanelsValue
      wbResult$TruncatedPlots <- pl
      output$DataPlot <- renderPlot({pl})
      AUCdf<-AUCfunction(AUCdf.new=AUCdf.new,wbResult$AUCdf,PanelsValue,SName = IDlane)
      
      output$AUC <- renderDT({
        AUCdf  %>% 
          dplyr::select(SampleName,Truncation, AUC) 
      },
      selection = 'none', 
      rownames= FALSE
      )
      wbResult$AUCdf <- AUCdf
    }
  })
  
  observeEvent(list(input$LaneChoice),{
    if(Flags$LanesCut & !is.null(wbResult$PanelsValue))
    {
      if(!is.null(wbResult$TruncatedPanelsValue ))
      {
        pl <- wbResult$TruncatedPlots    
        wbResult$TruncatedPanelsValue -> PanelsValue
      }
      else{
        wbResult$PanelsValue -> PanelsValue
        pl<-wbResult$Plots
      }
      
      Plots.Lane <- PanelsValue[which(PanelsValue$ID == input$LaneChoice),]
      colnames(Plots.Lane) = c("Y","ID","X")
      
      cat(input$LaneChoice,"\n")
      updateSliderInput(session,"truncV",
                        min = min(Plots.Lane$X),
                        max = max(Plots.Lane$X),
                        value = c(min(Plots.Lane$X), max(Plots.Lane$X) ) ) 
      updateSliderInput(session,"truncH",
                        min = round(min(Plots.Lane$Y),digits = 3),
                        max = round(max(Plots.Lane$Y),digits = 3),
                        value = round(min(Plots.Lane$Y),digits = 3) )
      
    }
    
  } )  
  
  observe({  Flags$CutTab <- input$tabs })
  
  observeEvent(c(input$truncV,input$truncH,input$LaneChoice), {
    if(!is.null(wbResult$PanelsValue))
    {
      if(!is.null(wbResult$TruncatedPanelsValue ))
      {
        pl <- wbResult$TruncatedPlots    
        wbResult$TruncatedPanelsValue -> PanelsValue
      }
      else{
        wbResult$PanelsValue -> PanelsValue
        pl<-wbResult$Plots
      }
      
      IDlane = input$LaneChoice
      Flags$IDlane <- IDlane
      
      if(Flags$CutTab=="V")
      {
        MinTrunc<-input$truncV[1]
        MaxTrunc<-input$truncV[2]
        
        vline.dat <- data.frame( ID = as.factor(rep(PanelsValue$ID,2)),
                                 vl = 0 )
        
        vline.dat  <-  vline.dat[ vline.dat$ID == IDlane, ]
        vline.dat$vl <- c(MinTrunc,MaxTrunc)
        
        pl <- pl + geom_vline(data=vline.dat,aes(xintercept=vl),linetype="dashed")
      }else if(Flags$CutTab=="H")
      {
        TruncY<-input$truncH[1]
        hline.dat <- data.frame(ID=as.factor(PanelsValue$ID), hl =min(PanelsValue$Y))
        hline.dat  <-  hline.dat[ hline.dat$ID == IDlane, ]
        hline.dat$hl <- TruncY
        
        pl <- pl + geom_hline(data=hline.dat,aes(yintercept = hl),linetype="dashed")
      }
      
      output$DataPlot <- renderPlot({pl})
    }  
  })
  
  observeEvent(c(input$actionButton_TruncV,input$actionButton_TruncH),{
    
    if( !is.null(wbResult$PanelsValue))
    {
      Flags$IDlane -> IDlane
      if(!is.null(wbResult$TruncatedPanelsValue ))
      {
        pl <- wbResult$TruncatedPlots    
        wbResult$TruncatedPanelsValue -> PanelsValue
      }
      else{
        wbResult$PanelsValue -> PanelsValue
        pl<-wbResult$Plots
      }
      
      maxPanelsValue=max(wbResult$PanelsValue$Values)
      wbResult$AUCdf -> AUCdf
      AUCdf.new <- AUCdf[length(AUCdf$Truncation),]
      #AUCdf.new$ExpName = "-"
      lastTrunc = AUCdf %>% 
        group_by(SampleName) %>%
        filter(SampleName == IDlane, row_number()==n() ) %>%
        ungroup() %>%
        dplyr::select(Truncation) 
      
      if(length(lastTrunc$Truncation) > 0 & lastTrunc$Truncation != "-")
        AUCdf.new$Truncation <- lastTrunc$Truncation
      
      AUCdf.new$SampleName <- IDlane
      
      if(Flags$CutTab=="V")
      {
        MinTrunc<-input$truncV[1]
        MaxTrunc<-input$truncV[2]
        AUCdf.new$Truncation <- paste(AUCdf.new$Truncation ,";\n X = [", MinTrunc," ; ", MaxTrunc ,"]",collapse = "")
        PanelsValue<- PanelsValue[!((PanelsValue$Y < MinTrunc | PanelsValue$Y > MaxTrunc) & PanelsValue$ID == IDlane),]
        PanelsValue$Values[PanelsValue$ID == IDlane] <- PanelsValue$Values[PanelsValue$ID == IDlane] -min(PanelsValue$Values[PanelsValue$ID == IDlane]) 
        # pl <- ggplot(PanelsValue, aes(x =Y,y=Values)) +
        #   geom_line() + theme_bw() +
        #   facet_wrap(~ID)
        updateSliderInput(session,"truncV",
                          min = min(PanelsValue$Y[PanelsValue$ID == IDlane]),
                          max = max(PanelsValue$Y[PanelsValue$ID == IDlane]),
                          value = c(min(PanelsValue$Y[PanelsValue$ID == IDlane]),
                                    max(PanelsValue$Y[PanelsValue$ID == IDlane]) ) )
      }
      else if(Flags$CutTab=="H")
      {
        TruncY<-input$truncH[1]
        PanelsValue <- PanelsValue[!(PanelsValue$Values<TruncY & PanelsValue$ID == IDlane),]
        PanelsValue$Values[PanelsValue$ID == IDlane] <- PanelsValue$Values[PanelsValue$ID == IDlane] - TruncY
        AUCdf.new$Truncation <- paste(AUCdf.new$Truncation ,";\n Y = ", TruncY)
        # pl <- ggplot(PanelsValue, aes(x =Y,y=Values)) +
        #   geom_line() + 
        #   theme_bw() +
        #   facet_wrap(~ID)+ 
        #   lims(y=c(minPanelsValue,maxPanelsValue))
        
        updateSliderInput(session,"truncH",
                          min = min(PanelsValue$Values[PanelsValue$ID == IDlane]),
                          max = max(PanelsValue$Values[PanelsValue$ID == IDlane]),
                          value = c(min(PanelsValue$Values[PanelsValue$ID == IDlane]),
                                    max(PanelsValue$Values[PanelsValue$ID == IDlane]) ) )
      }
      
      pl <- ggplot(PanelsValue, aes(x =Y,y=Values)) +
        geom_line() + 
        theme_bw() +
        facet_wrap(~ID)+ 
        lims(y=c(0,maxPanelsValue))
      
      wbResult$TruncatedPanelsValue <- PanelsValue
      wbResult$TruncatedPlots <- pl
      output$DataPlot <- renderPlot({pl})
      AUCdf<-AUCfunction(AUCdf.new=AUCdf.new,wbResult$AUCdf,PanelsValue,SName = IDlane)
      
      output$AUC <- renderDT({
        AUCdf  %>% 
          dplyr::select(SampleName,Truncation, AUC) 
      },
      selection = 'none', 
      rownames= FALSE
      # editable = list(target = "cell", 
      #                 disable = list(columns = 1:4))
      )
      wbResult$AUCdf <- AUCdf
    }
  })
  
  ## next buttons
  observeEvent(input$NextWBQuantif,{
    if(!is.null(wbResult$AUCdf))
      wbquantResult$WBanalysis = reactiveValuesToList(wbResult)
    
    updateTabsetPanel(session, "SideTabs",
                      selected = "quantification")
  })
  
  ## quantification WB
  wbquantResult = reactiveValues(NormWBanalysis = NULL,
                                 NormWBanalysis_filtered = NULL,
                                 WBanalysis = NULL,
                                 WBanalysis_filtered = NULL,
                                 RelDensitiy = NULL,
                                 AdjRelDensitiy = NULL
  )
  wbquantResult0 = list(NormWBanalysis = NULL,
                        NormWBanalysis_filtered = NULL,
                        WBanalysis = NULL,
                        WBanalysis_filtered = NULL,
                        RelDensitiy = NULL,
                        AdjRelDensitiy = NULL
  )
  FlagsWBquant = reactiveValues(BothUploaded = F)
  
  
  observeEvent(input$DataverseUpload_Button,{
    
    if(input$selectAnalysis_DV != ""){
      
      if(input$Title_DV == "" && input$Author_DV == "" && input$Description_DV == "" &&
         input$AuthorAff_DV== "" && input$ContactN_DV == "" && input$ContactEmail_DV == "")
        output$LoadingError_DATAVERSE = showAlert("Error", "Error: missing information", "error", 5000)
      else{
        # creation of a temporary folder
        tempfolder = paste0(tempdir(check = T),"/ORCA")
        system(paste0("mkdir ",tempfolder))
        system(paste0("mkdir ",tempfolder,"/dataset"))
        # create the metadata
        result <- rjson::fromJSON(file = system.file("docker","metadata.json", package = "ORCA") )
        result$dataset_title = input$Title_DV
        result$dataset_description = paste0(input$Description_DV,"\n This dataset has been obtained using the ORCA application, specifically the module: ", input$selectAnalysis_DV)
        result$author_name = input$Author_DV
        result$author_affiliation= input$AuthorAff_DV
        result$dataset_contact_name = input$ContactN_DV
        result$dataset_contact_email = input$ContactEmail_DV
        # result$subject = as.vector(result$subject)
        write(rjson::toJSON(result), file=paste0(tempfolder,"/metadata.json") )
        
        # move the file in the temporary folder
        
        saveExcel(filename = paste0(tempfolder,"/dataset/",
                                    gsub(pattern = " ", 
                                         x = input$selectAnalysis_DV,
                                         replacement = ""),".xlsx"),
                  ResultList = DataAnalysisModule[[ names(MapAnalysisNames[MapAnalysisNames == input$selectAnalysis_DV])]] ,
                  analysis = input$selectAnalysis_DV )
        
        saveRDS(DataAnalysisModule[[ names(MapAnalysisNames[MapAnalysisNames == input$selectAnalysis_DV])]] ,
                file = paste0(tempfolder,"/dataset/ORCA_",
                              gsub(pattern = " ", 
                                   x = input$selectAnalysis_DV,
                                   replacement = ""),".RDs"))
        # docker run
        ORCA::docker.run(params = paste0("--volume ", tempfolder,
                                         ":/Results/ -d qbioturin/orca-upload-dataverse python3 main.py /Results/metadata.json /Results/dataset") 
        )
        
      }
      
    }
  })
  
  observeEvent(input$NextWBQuantif,{
    if(!is.null(wbResult$AUCdf))
      wbquantResult$WBanalysis = reactiveValuesToList(wbResult)
    
    updateTabsetPanel(session, "SideTabs",
                      selected = "quantification")
  })
  
  wbquantResult = reactiveValues(NormWBanalysis = NULL,
                                 NormWBanalysis_filtered = NULL,
                                 WBanalysis = NULL,
                                 WBanalysis_filtered = NULL,
                                 RelDensitiy = NULL,
                                 AdjRelDensitiy = NULL
  )
  
  observeEvent(input$actionB_loadingNormWB, {
    mess = readfile(
      filename = input$NormWBImport$datapath,
      type = "RDs",
      namesAll = namesAll
    )
    
    if(is.list(mess) && !is.null(mess$message)) {
      showAlert("Error", mess[["message"]], "error", 5000)
      return() 
    }
    
    validate(
      need(!setequal(names(mess),c("message","call")) ,
           mess[["message"]])
    )
    
    choices = paste0(mess$AUCdf$SampleName, " with ", mess$AUCdf$Truncation)
    wbquantResult$NormWBanalysis = mess
    showAlert("Success", "The RDS has been uploaded with success", "success", 2000)
  })
  
  observeEvent(input$actionB_loadingWB,{
    mess = readfile(
        filename = input$WBImport$datapath,
        type = "RDs",
        namesAll = namesAll
      )
      
    if(is.list(mess) && !is.null(mess$message)) {
      showAlert("Error", mess[["message"]], "error", 5000)
      return() 
    }
    
    validate(
      need(!setequal(names(mess),c("message","call")) ,
           mess[["message"]])
    )
    
    wbquantResult$WBanalysis = mess
    wbquantResult$WBanalysis_filtered = NULL
    showAlert("Success", "The RDS has been uploaded with success", "success", 2000)
  })
  
  observe({
    if(!is.null(wbquantResult$WBanalysis) & !is.null(wbquantResult$NormWBanalysis))
      FlagsWBquant$BothUploaded = T
  })
  
  observe({
    if(is.null(wbquantResult$NormWBanalysis)){
      table = wbResult0$AUCdf
    }else{
      table = wbquantResult$NormWBanalysis$AUCdf
    }
    output$AUC_WBnorm <- renderDT(
      table , 
      #filter = 'top', server = FALSE, 
      selection = "multiple", 
      # editable = list(target = "cell", 
      #                 disable = list(columns = 1:3)),
      options = list(lengthChange = FALSE, autoWidth = TRUE),
      rownames= FALSE
    )
  })
  
  observe({
    if(is.null(wbquantResult$WBanalysis)){
      table = wbResult0$AUCdf
    }else{
      table = wbquantResult$WBanalysis$AUCdf
    }
    output$AUC_WB <- renderDT(
      table,
      #filter = 'top', server = FALSE, 
      selection = "multiple", 
      # editable = list(target = "cell", 
      #                 disable = list(columns = 1:3)),
      options = list(lengthChange = FALSE, autoWidth = TRUE),
      rownames= FALSE
    )
  })
  
  observeEvent(input$AUC_WB_rows_selected,{
    if(!is.null(wbquantResult$WBanalysis) ){
      indexesWB = input$AUC_WB_rows_selected
      AUCdf = wbquantResult$WBanalysis$AUCdf
      
      if(length(indexesWB) > 0){
        wbquantResult$WBanalysis_filtered = AUCdf[indexesWB,]
      }else{
        wbquantResult$WBanalysis_filtered = AUCdf
      }
    }
  })
  
  observeEvent(input$AUC_WBnorm_rows_selected,{
    if(!is.null(wbquantResult$NormWBanalysis)){
      indexesWB = input$AUC_WBnorm_rows_selected
      AUCdf = wbquantResult$NormWBanalysis$AUCdf
      
      if(length(indexesWB) > 0){
        wbquantResult$NormWBanalysis_filtered = AUCdf[indexesWB,]
        
        choices = paste0( AUCdf[indexesWB,]$SampleName, "; truncated ", AUCdf[indexesWB,]$Truncation)
        selected = input$IdLaneNorm_RelDens
        updateSelectInput("IdLaneNorm_RelDens",
                          session = session,
                          choices = choices,
                          selected = selected)
      }else{
        wbquantResult$NormWBanalysis_filtered = AUCdf
      }
    }
  })
  
  observeEvent(list(FlagsWBquant$BothUploaded, input$IdLaneNorm_RelDens,input$AUC_WB_rows_selected,input$AUC_WBnorm_rows_selected),{
    table =  wbResult0$AUCdf 
    
    if(!is.null(wbquantResult$WBanalysis_filtered) & !is.null(wbquantResult$NormWBanalysis_filtered)){
      IdLaneNorm_RelDens = input$IdLaneNorm_RelDens
      IdLaneNorm_RelDens = strsplit(IdLaneNorm_RelDens,
                                    split = "; truncated ")[[1]]
      
      tbWBnorm = wbquantResult$NormWBanalysis_filtered %>%
        filter(SampleName ==IdLaneNorm_RelDens[1],
               Truncation == IdLaneNorm_RelDens[2]) %>%
        rename(AUC_Norm = AUC,
               Truncation_Norm = Truncation,
               SampleName_Norm = SampleName)
      
      tbWB = wbquantResult$WBanalysis_filtered
      
      if(!is.null(tbWBnorm) & !is.null(tbWB) & dim(tbWBnorm)[1]==1 ){
        if(!all(table(tbWB$SampleName)==1) ){
          output$LoadingErrorWB <- renderText({"No rows with equal sample name are allowed"})
        }
        else{ 
          table = tbWB
          table$AUC_Norm = tbWBnorm$AUC_Norm
          table$RelDens = table$AUC/table$AUC_Norm
          table = table %>%
            dplyr::select(SampleName, Truncation, AUC, AUC_Norm, RelDens) 
        }
      }
    }
    
    wbquantResult$RelDensitiy = table
    
    output$AUC_RelDens <- renderDT(
      table,
      filter = 'top',
      server = FALSE,
      options = list(lengthChange = FALSE, autoWidth = TRUE),
      rownames= FALSE
    )
    
  })
  
  observeEvent(list(FlagsWBquant$BothUploaded, input$AUC_WB_rows_selected,input$AUC_WBnorm_rows_selected),{
    table = data.frame(SampleName = "-",
                       Truncation = "-", 
                       Truncation_Norm = "-",
                       AUC = "-", 
                       AUC_Norm = "-",
                       AdjRelDens = "-")
    
    if(!is.null(wbquantResult$WBanalysis_filtered) & !is.null(wbquantResult$NormWBanalysis_filtered)){
      
      tbWB = wbquantResult$WBanalysis_filtered
      tbWBnorm = wbquantResult$NormWBanalysis_filtered
      
      if(!all(table(tbWBnorm$SampleName)==1) ){
        output$LoadingErrorNormWB <- renderText({"No rows with equal sample name are allowed"})
      }else if(!all(table(tbWB$SampleName)==1) ){
        output$LoadingErrorWB <- renderText({"No rows with equal sample name are allowed"})
      }
      else{ # we admit only one SampleName
        
        tbWBnorm = tbWBnorm  %>%
          rename(AUC_Norm = AUC,
                 Truncation_Norm = Truncation)
        
        table = merge(tbWBnorm,tbWB, by = "SampleName" ,all = T )
        
        table$AdjRelDens = table$AUC/table$AUC_Norm
        table = table %>% 
          dplyr::select( SampleName, Truncation, Truncation_Norm, AUC, AUC_Norm, AdjRelDens) 
        
        wbquantResult$AdjRelDensitiy = table
        output$AUC_AdjRelDens <- renderDT(
          table ,
          server = FALSE,
          options = list(lengthChange = FALSE, autoWidth = TRUE),
          rownames= FALSE
        )
        
        if(dim(table)[1] > 1 ){
          barPlotAdjRelDens = table %>% 
            mutate(Normalizer = paste0("Sample: ",SampleName ),
                   WB = paste0("Sample: ",SampleName))  %>%
            ggplot() +
            geom_bar(aes(x = SampleName,
                         y = AdjRelDens,
                         fill = Normalizer ),
                     stat = "identity" ) +
            #facet_grid(~WB)+
            theme_bw()
        }else{
          barPlotAdjRelDens = ggplot()
        }
        output$plot_AdjRelDens <- renderPlot({
          barPlotAdjRelDens
        })
      }
    }
  })
  
  output$downloadWBquantAnalysis <- downloadHandler(
    filename = function() {
      paste0('WBquantanalysis-', Sys.Date(), '.zip')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      tempRDSPath <- file.path(tempDir, paste0("Analysis-", Sys.Date(), ".rds"))
      tempExcelPath <- file.path(tempDir, paste0("Analysis-", Sys.Date(), ".xlsx"))
      
      resultsRDS <- DataAnalysisModule$wbquantResult
      saveRDS(resultsRDS, tempRDSPath)
      
      resultsExcel <- DataAnalysisModule$wbquantResult
      saveExcel(filename = tempExcelPath, ResultList=resultsExcel, analysis = "WB comparison")
      
      zip(file, files = c(tempRDSPath, tempExcelPath), flags = "-j")
      manageSpinner(FALSE)
    }
  )
  
  toListenWBquant <- reactive({
    reactiveValuesToList(wbquantResult)
  })
  observeEvent(toListenWBquant(),{
    DataAnalysisModule$wbquantResult = reactiveValuesToList(wbquantResult)
  })
  
  ### End WB analysis ####
  
  #### PCR analysis ####
  
  pcrResult = reactiveValues(
    Initdata = NULL,
    selectPCRcolumns = NULL,
    data = NULL,
    PCRnorm = NULL,
    BaselineExp = NULL,
    plotPRC = NULL,
    NewPCR = NULL)
  pcrResult0 = list(
    Initdata = NULL,
    selectPCRcolumns = NULL,
    data = NULL,
    PCRnorm = NULL,
    BaselineExp = NULL,
    plotPRC = NULL,
    NewPCR = NULL)
  
  # save everytime there is a change in the results
  PCRresultListen <- reactive({
    reactiveValuesToList(pcrResult)
  })
  observeEvent(PCRresultListen(), {
    DataAnalysisModule$pcrResult = reactiveValuesToList(pcrResult)
    DataAnalysisModule$pcrResult$Flags = reactiveValuesToList(FlagsPCR)
  })
  
  ## next buttons
  observeEvent(input$NextQuantif,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "tablesPCR")
  })
  observeEvent(input$NextpcrPlots,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "plotsPCR")
  })
  
  FlagsPCR <- reactiveValues(norm=F, 
                             baseline = F)
  
  
  observeEvent(input$LoadPCR_Button,{
    alert$alertContext <- "PCR-reset"
    if( !is.null(pcrResult$Initdata) ) {
      shinyalert(
        title = "Important message",
        text = "Do you want to update the WB data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else loadExcelFilePCR()
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "PCR-reset") {  
      resetPanel("PCR", flags = FlagsPCR, result = pcrResult)
      loadExcelFilePCR()
    }
  })
  
  loadExcelFilePCR <- function() {
    alert$alertContext <- ""
    
    mess = readfile(
      filename = input$PCRImport$datapath,
      type = "Excel", 
      isFileUploaded = !is.null(input$PCRImport) && file.exists(input$PCRImport$datapath),
    )
    
    
    if (!is.null(mess$message) && mess$call == "") {
      showAlert("Error", mess$message, "error", 5000)
      return()
    }
    
    pcrResult$Initdata = mess
    
    updateSelectInput(session,
                      "PCR_gene",
                      choices = c("", colnames(pcrResult$Initdata)),
                      selected = "")
    updateSelectInput(session,
                      "PCR_sample",
                      choices = c("", colnames(pcrResult$Initdata)),
                      selected = "")
    updateSelectInput(session,
                      "PCR_value",
                      choices = c("", colnames(pcrResult$Initdata)),
                      selected = "")
    updateSelectInput(session,
                      "PCR_time",
                      choices = c("", colnames(pcrResult$Initdata)),
                      selected = "")
    showAlert("Success", "The RT-qPCR excel has been uploaded with success", "success", 2000)
    
  }
  
  observeEvent(list(input$PCR_value,input$PCR_gene,input$PCR_sample,input$PCR_time),{
    if( !is.null(pcrResult$Initdata) ){
      selectPCRcolumns = c(input$PCR_gene,input$PCR_sample,input$PCR_value,input$PCR_time)
      selectPCRcolumns = selectPCRcolumns[selectPCRcolumns!= ""]
      
      PCR = pcrResult$Initdata
      colNames = colnames(PCR)
      output$PCRpreview = renderTable({
        if(length(selectPCRcolumns)!=0 ){
          tmp = PCR[,selectPCRcolumns]
          #colnames(tmp) = c("Gene", "Sample", "Value")[1:length(colnames(tmp))]
          head(tmp) 
        }
        else
          NULL
      })
      
      if(length(selectPCRcolumns)==3 ){
        tmp = PCR[,selectPCRcolumns]
        colnames(tmp) = c("Gene", "Sample", "Value")
        tmp$Time = ""
        pcrResult$data = tmp
        pcrResult$selectPCRcolumns = selectPCRcolumns
      }else if(length(selectPCRcolumns)==4 ){
        tmp = PCR[,selectPCRcolumns]
        colnames(tmp) = c("Gene", "Sample", "Value","Time")
        pcrResult$data = tmp
        pcrResult$selectPCRcolumns = selectPCRcolumns
      }else{
        pcrResult$data = NULL
      }
    }
    
  })
  
  observe({
    if(!is.null(pcrResult$data)){
      
      PCR = pcrResult$data
      
      AllGenes = unique(PCR$Gene)
      Exp = unique(PCR$Sample)
      
      updateSelectInput(session, "PCRbaseline",
                        choices = Exp )
      updateCheckboxGroupInput(session,"PCRnorm",
                               choices = AllGenes )
      
    }else{
      updateSelectInput(session, "PCRbaseline",
                        choices = "" )
      updateCheckboxGroupInput(session,"PCRnorm",
                               choices = "" )
    }
  })
  
  observeEvent(input$PCRnorm,{
    pcrResult$PCRnorm = input$PCRnorm
    FlagsPCR$norm = T
  })
  observeEvent(input$PCRbaseline,{
    pcrResult$BaselineExp = input$PCRbaseline
    FlagsPCR$baseline = T
  })
  
  observe({
    if(FlagsPCR$baseline & FlagsPCR$norm & !is.null(pcrResult$data)){
      
      pcrResult$BaselineExp -> BaselineExp
      pcrResult$PCRnorm -> PCRnorm
      pcrResult$data -> PCR
      
      NewPCR = PCR %>% 
        na.omit()%>%
        group_by(Sample,Gene,Time) %>%
        dplyr::summarise(Mean = mean(Value),
                         Sd = sd(Value)) %>%
        ungroup()
      
      HousekGenePCR = NewPCR %>%
        filter(Gene %in% PCRnorm)%>%
        rename(HousekGene = Gene,HousekGeneMean=Mean, HousekGeneSd=Sd) 
      
      PCRstep2 = merge(HousekGenePCR,NewPCR %>% filter(!Gene %in% PCRnorm),all.y = T,by=c("Sample","Time") )
      
      #PCRstep3 = merge(BaselinePCR,PCRstep2,all.y = T,by=c("Gene","Time") )
      
      
      PCRstep3 = PCRstep2 %>%
        group_by(Sample,Gene,Time) %>%
        dplyr::mutate(dCt = Mean - HousekGeneMean)%>%
        ungroup()
      
      BaselinePCR = PCRstep3 %>% 
        filter(Sample == BaselineExp) %>%
        rename(BaselineMean=Mean, BaselineSd=Sd,BaselinedCt = dCt) %>%
        dplyr::select(-Sample, -HousekGene, -HousekGeneMean, -HousekGeneSd)
      
      PCRstep4 = merge(BaselinePCR,PCRstep3,all.y = T,by=c("Gene","Time") )
      
      PCRstep5 = PCRstep4 %>%
        group_by(Sample,Gene,Time) %>%
        dplyr::summarize(
          ddCt = dCt - BaselinedCt,
          Q = 2^{-ddCt},
          Sd = Sd,
          Mean = Mean)%>%
        ungroup()
      
      
      NormPCR = PCRstep5 %>%
        filter(Gene %in% PCRnorm ) %>%
        rename(Norm = Gene,
               NormQ = Q,
               NormSd = Sd,
               NormMean = Mean)
      
      # CompPRC = merge(OnePCR,NormPCR)
      # 
      # CompPRC = CompPRC %>% group_by(Sample,Gene,Norm) %>%
      #   dplyr::summarise(Qnorm = Q/NormQ,
      #                    SDddct = sqrt(Sd^2+NormSd^2),
      #                    SDrq = log(2)*Qnorm*SDddct) %>%
      #   ungroup()
      
      AllGenes = unique(PCR$Gene)
      
      #pcrResult$CompPRC = CompPRC
      pcrResult$NewPCR = PCRstep5
      
      output$PCRtables <- renderUI({
        plot_output_list <- lapply(AllGenes, function(i) {
          tablename <- paste("tablename", i, sep="")
          tableOutput(tablename)
        })
        do.call(tagList, plot_output_list)
      })
      
      # output$PCRtablesComp <- renderUI({
      #   plot_output_list <- lapply(AllGenes[-which(AllGenes %in% PCRnorm)], function(i) {
      #     tablename <- paste("CompTablename", i, sep="")
      #     tableOutput(tablename)
      #   })
      #   do.call(tagList, plot_output_list)
      # })
      
      plot1 = ggplot(data = PCRstep5,
                     aes(x= as.factor(Time), y = ddCt, col = Sample)) + 
        facet_wrap(~Gene, ncol = 1) +
        geom_jitter(width = 0.1, height = 0,size = 2)+
        theme_bw()+
        labs(x = "Time", y = "DDCT")
      
      plot2 = ggplot(data = PCRstep5,
                     aes(x= as.factor(Time), y = Q, col = Sample)) + 
        facet_wrap(~Gene, ncol = 1) +
        geom_jitter(width = 0.1, height = 0,size = 2)+
        theme_bw()+
        labs(x = "Time", y = "2^(-DDCT)")
      
      pcrResult$plotPRC =  plot1/plot2
      
      output$PCRplot <- renderPlot({
        pcrResult$plotPRC
      })
      
      for (i in AllGenes[!AllGenes %in% PCRnorm]){
        local({
          my_i <- i
          tablename <- paste("tablename", my_i, sep="")
          output[[tablename]] <- renderTable({
            PCRstep5 %>% filter(Gene == my_i) %>% rename(DDCT = ddCt, `2^(-DDCT)` = Q)
          })
          
          # ComparisonPCR = list()
          # if(my_i %in% AllGenes[-which(AllGenes %in% PCRnorm)]){
          #   tablename <- paste("CompTablename", my_i, sep="")
          #   output[[tablename]] <- renderTable({
          #     CompPRC %>% 
          #       filter(Gene == my_i) %>%
          #       arrange(Norm,Sample)
          #   })
          # }
        })    
      }
      
    }
  })
  
  observe({
    DataAnalysisModule$pcrResult = reactiveValuesToList(pcrResult)
  })
  
  output$downloadRTPCRAnalysis <- downloadHandler(
    filename = function() {
      paste('RTqPCRanalysis-', Sys.Date(), '.zip', sep='')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      
      nomeRDS <- paste0("RTqPCRanalysis-", Sys.Date(), ".rds")
      nomeXLSX <- paste0("RTqPCRanalysis-", Sys.Date(), ".xlsx")
      
      tempRdsPath <- file.path(tempDir, nomeRDS)
      tempXlsxPath <- file.path(tempDir, nomeXLSX)
      
      results <- DataAnalysisModule$pcrResult
      saveRDS(results, file = tempRdsPath)
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "RT-qPCR")
      
      zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
    },
  )
  
  #### END PCR analysis ####
  
  #### ELISA analysis ####
  observeEvent(input$NextElisaQuantif,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "tablesELISA")
  })
  #
  
  elisaResult = reactiveValues(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    ELISAcell_EXP = "",
    ELISAcell_SN = NULL,
    ELISAcell_COLOR = NULL,
    MapBaseline = NULL,
    MapBlank = NULL,
    Tablestandcurve = NULL,
    Regression = NULL)
  
  elisaResult0 = list(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    ELISAcell_EXP = "",
    ELISAcell_SN = NULL,
    ELISAcell_COLOR = NULL,
    MapBaseline = NULL,
    MapBlank = NULL,
    Tablestandcurve = NULL,
    Regression = NULL)
  
  left_data_elisa <- reactiveVal()
  right_data_elisa <- reactiveVal()
  
  # save everytime there is a change in the results
  ELISAresultListen <- reactive({
    reactiveValuesToList(elisaResult)
  })
  observeEvent(ELISAresultListen(), {
    DataAnalysisModule$elisaResult = reactiveValuesToList(elisaResult)
    DataAnalysisModule$elisaResult$Flags = reactiveValuesToList(FlagsELISA)
  })
  
  ##
  FlagsELISA <- reactiveValues(cellCoo = NULL,
                               AllExp = "",
                               BASEselected = "",
                               STDCselected = "",
                               BLANCHEselected = "",
                               EXPselected = "",
                               EXPcol = NULL)
  
  observeEvent(input$LoadELISA_Button,{
    alert$alertContext <- "ELISA-reset"
    if( !is.null(elisaResult$Initdata) ) {
      shinyalert(
        title = "Important message",
        text = "Do you want to update the ELISA data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else loadExcelFileELISA()
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "ELISA-reset") {  
      resetPanel("ELISA", flags = FlagsELISA, result = elisaResult)
      
      loadExcelFileELISA()
    }
  })
  
  loadExcelFileELISA <- function() {
    alert$alertContext <- ""
    mess = readfile(
      filename = input$ELISAImport$datapath,
      isFileUploaded = !is.null(input$ELISAImport) && file.exists(input$ELISAImport$datapath),
      type = "Excel",
      allDouble = T,
      colname = F,
      colors = T
    )
    
    if(setequal(names(mess), c("message", "call"))) {
      showAlert("Error", mess[["message"]], "error", 5000)
    } else {
      elisaResult$Initdata = mess$x
      FlagsELISA$EXPcol = mess$fill
      elisaResult$ELISAcell_COLOR = mess$SNtable
      elisaResult$ELISAcell_SN = matrix("", nrow = nrow(elisaResult$ELISAcell_COLOR), ncol = ncol(elisaResult$ELISAcell_COLOR))
      
      removeModal()
      showAlert("Success", "The Excel has been uploaded  with success", "success", 2000)
    }
  }
  
  observe({
    if (!is.null(elisaResult$Initdata) && is.null(elisaResult$TablePlot)) {
      tableExcelColored(session = session,
                        Result = elisaResult, 
                        FlagsExp = FlagsELISA,
                        type = "Initialize")
    }
  })
  
  observeEvent(c(elisaResult$TablePlot,elisaResult$ELISAcell_EXP), {
    if (!is.null(elisaResult$TablePlot)) {
      ELISAtb <- elisaResult$TablePlot
      output$ELISAmatrix <- renderDT(ELISAtb, server = FALSE)
      
      if (!is.null(elisaResult$ELISAcell_EXP) && !is.null(elisaResult$ELISAcell_COLOR)) {
        matTime <- as.matrix(elisaResult$ELISAcell_EXP)
        matExp <- as.matrix(elisaResult$ELISAcell_COLOR)
        
        if (!(all(matTime == "") || all(matExp == ""))) {
          mat <- as.matrix(elisaResult$Initdata)
          elisaV <- expand.grid(seq_len(nrow(mat)), seq_len(ncol(mat))) %>%
            rowwise() %>%
            mutate(values = mat[Var1, Var2])
          elisaT <- expand.grid(seq_len(nrow(matTime)), seq_len(ncol(matTime))) %>%
            rowwise() %>%
            mutate(time = matTime[Var1, Var2])
          elisaE <- expand.grid(seq_len(nrow(matExp)), seq_len(ncol(matExp))) %>%
            rowwise() %>%
            mutate(exp = matExp[Var1, Var2])
          elisaTot <- merge(elisaV, merge(elisaT, elisaE)) %>%
            na.omit() %>%
            filter(time != "", exp != "")
          
          elisaResult$data <- elisaTot
          
          output$ELISAinitplots <- renderPlot({
            ggplot(elisaTot, aes(x = time, y = values, col = exp), alpha = 1.4) +
              geom_point(aes(group = exp)) +
              scale_color_manual(values = FlagsELISA$EXPcol) +
              theme_bw() +
              labs(x = "Times", y = "Values", col = "Exp", fill = "Exp") +
              theme(legend.position = c(0, 1),
                    legend.justification = c(0, 1),
                    legend.direction = "vertical",
                    legend.background = element_rect(linewidth = 0.5,
                                                     linetype = "solid",
                                                     colour = "black"))
          })
        }
      }
    } 
  })
  
  observe({
    color_codes <- FlagsELISA$EXPcol
    color_names <- names(FlagsELISA$EXPcol)
    
    valid_colors <- color_codes != "white"
    color_codes <- color_codes[valid_colors]
    color_names <- color_names[valid_colors]
    
    mid_point <- ceiling(length(color_codes) / 2)
    left_colors <- color_codes[1:mid_point]
    right_colors <- color_codes[(mid_point+1):length(color_codes)]

    left_formatted_data <- get_formatted_data(left_colors, color_names[1:mid_point], elisaResult, elisaResult$ELISAcell_EXP,"ELISA")
    right_formatted_data <- get_formatted_data(right_colors, color_names[(mid_point+1):length(color_codes)], elisaResult, elisaResult$ELISAcell_EXP, "ELISA")
    
    left_data_elisa(left_formatted_data)
    right_data_elisa(right_formatted_data)

    output$leftTableElisa <- renderDataTable(
      left_data_elisa(), 
      escape = FALSE, 
      editable = list(target = "cell", disable = list(columns = 0:3)),
      options = list(
        dom = 't',
        paging = FALSE,
        info = FALSE,
        searching = FALSE, 
        columnDefs = list(
          list(targets = 0, visible = FALSE),
          list(targets = 1, visible = FALSE),
          list(width = '10px', targets = 2),
          list(width = '200px', targets = 3),
          list(width = '80px', targets = 4),
          list(width = '100px', targets = 5),
          list(className = 'dt-head-left dt-body-left', targets = 1)
        )
      )
    )
    
    output$rightTableElisa <- renderDataTable(
      right_data_elisa(), 
      escape = FALSE, 
      editable = list(target = "cell", disable = list(columns = 0:3)),
      options = list(
        dom = 't',
        paging = FALSE,
        info = FALSE,
        searching = FALSE,
        editable = TRUE,
        columnDefs = list(
          list(targets = 0, visible = FALSE),
          list(targets = 1, visible = FALSE),
          list(width = '10px', targets = 2),
          list(width = '200px', targets = 3),
          list(width = '80px', targets = 4),
          list(width = '100px', targets = 5),
          list(className = 'dt-head-left dt-body-left', targets = 1)
        )
      )
    )
  })
  
  observeEvent(input$leftTableElisa_cell_edit, {
    info <- input$leftTableElisa_cell_edit
    data <- left_data_elisa() 
    updatedText <- updateTable("left", "ELISA", info, data, elisaResult, FlagsELISA)
    
    output$ELISASelectedValues <- renderText(updatedText)  
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent(input$rightTableElisa_cell_edit, {
    info <- input$rightTableElisa_cell_edit
    data <- right_data_elisa() 
    updatedText <- updateTable("right", "ELISA", info, data, elisaResult, FlagsELISA)
    
    output$ELISASelectedValues <- renderText(updatedText)  
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent(input$ELISAmatrix_cell_clicked, {
    req(input$ELISAmatrix_cell_clicked)  
    
    cellSelected = as.numeric(input$ELISAmatrix_cell_clicked)
    FlagsELISA$cellCoo = cellCoo = c(cellSelected[1], cellSelected[2] + 1)
    
    allExp <- unique(na.omit(c(elisaResult$ELISAcell_EXP)))  
    selectedExp <- ifelse(is.null(elisaResult$ELISAcell_EXP[cellCoo[1], cellCoo[2]]), "", elisaResult$ELISAcell_EXP[cellCoo[1], cellCoo[2]])
    
    updateSelectizeInput(inputId = "ELISAcell_EXP", 
                         choices = allExp,
                         selected = selectedExp)
    
    allSN <- unique(na.omit(c(elisaResult$ELISAcell_SN)))  
    selectedSN <- ifelse(is.null(elisaResult$ELISAcell_SN[cellCoo[1], cellCoo[2]]), "", elisaResult$ELISAcell_SN[cellCoo[1], cellCoo[2]])
    
    updateSelectizeInput(inputId = "ELISAcell_SN",
                         choices = allSN,
                         selected = selectedSN)
  })
  
  observeEvent(input$ELISAcell_SN, {
    if (!is.null(elisaResult$ELISAcell_COLOR) && !is.null(FlagsELISA$cellCoo) && !anyNA(FlagsELISA$cellCoo)) {
      ELISAtb = elisaResult$TablePlot
      cellCoo = FlagsELISA$cellCoo
      
      value.bef = elisaResult$ELISAcell_COLOR[cellCoo[1], cellCoo[2]] 
      value.now = input$ELISAcell_SN
      
      if (value.now != "" && value.now != value.bef) {
        currentValues <- elisaResult$Initdata[cellCoo[1], cellCoo[2]]
        
        elisaResult$ELISAcell_COLOR[cellCoo[1], cellCoo[2]] = value.now
        elisaResult$ELISAcell_SN[cellCoo[1], cellCoo[2]] = value.now
        ELISAtb$x$data[cellCoo[1], paste0("Col", cellCoo[2])] = value.now
        
        if (!input$ELISAcell_SN %in% FlagsELISA$AllExp) {
          exp = unique(c(FlagsELISA$AllExp, input$ELISAcell_SN))
          FlagsELISA$AllExp = exp
        }
        
        tableExcelColored(session = session,
                          Result = elisaResult, 
                          FlagsExp = FlagsELISA,
                          type = "Update")
        
        output$ELISASelectedValues <- renderText(paste("Updated value", paste(currentValues), ": sample name ", value.now))
        output$ELISAmatrix <- renderDataTable({elisaResult$TablePlot})
      }
    } else return()
  }, ignoreInit = TRUE)
  
  observeEvent(input$ELISAcell_EXP, {
    if (!is.null(elisaResult$ELISAcell_EXP) && !is.null(FlagsELISA$cellCoo) && !anyNA(FlagsELISA$cellCoo)) {
      ELISAtb = elisaResult$TablePlot
      cellCoo = FlagsELISA$cellCoo
      
      value.bef = elisaResult$ELISAcell_EXP[cellCoo[1], cellCoo[2]] 
      value.now = input$ELISAcell_EXP
      
      if (value.now != "" && value.now != value.bef) {
        currentValues <- elisaResult$Initdata[cellCoo[1], cellCoo[2]]
        
        elisaResult$ELISAcell_EXP[cellCoo[1], cellCoo[2]] = value.now
        tableExcelColored(session = session,
                          Result = elisaResult, 
                          FlagsExp = FlagsELISA,
                          type = "Update",
        )
        
        output$ELISASelectedValues <- renderText(paste("Updated value", paste(currentValues), ": time ", value.now))
        output$ELISAmatrix <- renderDataTable({elisaResult$TablePlot})
      }
    } else return()
  }, ignoreInit = TRUE)
  
  ## update Baselines checkBox
  observeEvent(c(FlagsELISA$AllExp,FlagsELISA$BASEselected,FlagsELISA$BLANCHEselected),{
    if(length(FlagsELISA$AllExp) > 1){
      exp = FlagsELISA$AllExp
      exp = exp[exp != ""]
      
      bool.tmp = exp %in% unique(c(FlagsELISA$BLANCHEselected,FlagsELISA$BASEselected))
      if( length(bool.tmp) > 0  )
        exp = exp[!bool.tmp]
      
      updateSelectizeInput(session,"ELISA_standcurve",
                           choices = exp,
                           selected = ifelse(FlagsELISA$STDCselected %in% exp,FlagsELISA$STDCselected,"") 
      )
    }
  })
  observeEvent(c(FlagsELISA$AllExp,FlagsELISA$BASEselected,FlagsELISA$STDCselected),{
    if(length(FlagsELISA$AllExp) > 1){
      exp = FlagsELISA$AllExp
      exp = exp[exp != ""]
      
      bool.tmp = exp %in% unique(c(FlagsELISA$STDCselected,FlagsELISA$BASEselected))
      if( length(bool.tmp) > 0  )
        exp = exp[!bool.tmp]
      
      updateCheckboxGroupInput(session,"ELISA_blanks",
                               choices = exp,
                               selected = FlagsELISA$BLANCHEselected )
    }
  })
  
  observeEvent(c(FlagsELISA$AllExp,FlagsELISA$BLANCHEselected,FlagsELISA$STDCselected),{
    if(length(FlagsELISA$AllExp) > 1){
      exp = FlagsELISA$AllExp
      exp = exp[exp != ""]
      
      bool.tmp = exp %in% unique(c(FlagsELISA$STDCselected,FlagsELISA$BLANCHEselected))
      if( length(bool.tmp) > 0  )
        exp = exp[!bool.tmp]
      
      exp_selec = input$ELISA_baselines
      
      updateCheckboxGroupInput(session,"ELISA_baselines",
                               choices = exp,
                               selected = FlagsELISA$BASEselected )
    }
  })
  
  ## select the baselines, std curves, and blank
  observeEvent(input$ELISA_baselines,{
    FlagsELISA$BASEselected = input$ELISA_baselines
    FlagsELISA$EXPselected = FlagsELISA$AllExp[! FlagsELISA$AllExp %in% c(FlagsELISA$STDCselected,FlagsELISA$BASEselected,FlagsELISA$BLANCHEselected)]
  },ignoreNULL = F)
  observeEvent(input$ELISA_standcurve,{
    FlagsELISA$STDCselected = input$ELISA_standcurve
    FlagsELISA$EXPselected = FlagsELISA$AllExp[! FlagsELISA$AllExp %in% c(FlagsELISA$STDCselected,FlagsELISA$BASEselected,FlagsELISA$BLANCHEselected)]
  },ignoreNULL = F)
  observeEvent(input$ELISA_blanks,{
    FlagsELISA$BLANCHEselected = input$ELISA_blanks
    FlagsELISA$EXPselected = FlagsELISA$AllExp[! FlagsELISA$AllExp %in% c(FlagsELISA$STDCselected,FlagsELISA$BASEselected,FlagsELISA$BLANCHEselected)]
  },ignoreNULL = F)
  
  toListen_elisa <- reactive({
    exp = FlagsELISA$EXPselected
    exp = exp[exp != ""]
    if(length(exp) > 0 )
    {
      Input_baselEXP = lapply(exp,
                              function(i) input[[paste0("elisa_Exp",i)]])
      Input_blEXP = lapply(unique(exp,FlagsELISA$BASELINEselected),
                           function(i) input[[paste0("elisa_blExp",i)]] )
      InputEXP = c(Input_baselEXP,Input_blEXP)
      
      which(sapply(InputEXP, function(x) 
        ifelse(is.null(x), T, ifelse(x == "", T, F) ) ) ) -> indexesEXPnull
      
      if(length(indexesEXPnull) > 0 )
        listReturn = InputEXP[-indexesEXPnull]
      else
        listReturn = InputEXP
    }else{
      listReturn = list()
    }
    
    if(length(listReturn) == 0){
      return(list("Nothing",elisaResult$ELISAcell_EXP,elisaResult$ELISAcell_SN))
    }else{
      return(c(listReturn,list(elisaResult$ELISAcell_EXP,elisaResult$ELISAcell_SN)) )
    }
  })
  
  observeEvent(toListen_elisa(),{
    baselines = FlagsELISA$BASEselected
    baselines = baselines[baselines != ""]
    
    if(toListen_elisa()[[1]] != "Nothing" ){
      exp = FlagsELISA$EXPselected
      exp = exp[exp != ""]
      expNotBlank = unique(c(exp,baselines))
      
      MapBaseline = do.call(rbind,
                            lapply(exp,function(i){
                              if( length(input[[paste0("elisa_Exp",i)]]) > 0 && input[[paste0("elisa_Exp",i)]] != ""){
                                data.frame(Exp = i, Baseline = input[[paste0("elisa_Exp",i)]])
                              }else{
                                data.frame(Exp = i, Baseline = NA)
                              }
                            })
      ) %>% na.omit()
      
      MapBlank = do.call(rbind,
                         lapply(expNotBlank,
                                function(i){
                                  if( length(input[[paste0("elisa_blExp",i)]]) > 0 && input[[paste0("elisa_blExp",i)]] != ""){
                                    data.frame(Exp = i, Blank = input[[paste0("elisa_blExp",i)]])
                                  }else{
                                    data.frame(Exp = i, Blank = NA)
                                  }
                                })
      ) %>% na.omit()
      
      elisaResult$MapBaseline = MapBaseline
      elisaResult$MapBlank = MapBlank
      
      mat = as.matrix(elisaResult$Initdata)
      elisaV = expand.grid(seq_len(nrow(mat)), seq_len(ncol(mat))) %>%
        rowwise() %>%
        mutate(values = mat[Var1, Var2])
      matTime =  as.matrix(elisaResult$ELISAcell_EXP)
      elisaT = expand.grid(seq_len(nrow(matTime)), seq_len(ncol(matTime))) %>%
        rowwise() %>%
        mutate(time = matTime[Var1, Var2])
      matExp =  as.matrix(elisaResult$ELISAcell_SN)
      elisaE = expand.grid(seq_len(nrow(matExp)), seq_len(ncol(matExp))) %>%
        rowwise() %>%
        mutate(exp = matExp[Var1, Var2])
      elisaTot = merge(elisaV,merge(elisaT,elisaE)) %>%
        filter(exp != "")
      
      elisaTotAverage = elisaTot %>%
        #mutate(time = ifelse(exp %in% MapBlank$Blank, 0, time)) %>%
        group_by(time, exp) %>%
        summarize(meanValues = mean(values))
      
      # merging exp with blank for the substraction
      
      elisaTot_bl = right_join( elisaTotAverage,MapBlank, 
                                by= c("exp"= "Blank") )%>%
        rename(BlankValues = meanValues, Blank =  exp, exp = Exp )
      
      elisaTotAverage = merge( elisaTotAverage %>% filter( exp %in%elisaTot_bl$exp ),
                               elisaTot_bl %>% ungroup(),all.x = T, by = c("exp","time") ) 
      elisaTotAverage[is.na(elisaTotAverage[,])] = 0
      elisaTotAverage = elisaTotAverage %>% mutate(meanValues = meanValues - BlankValues )
      
      # merging exp with baseline
      elisaTot_base = merge(MapBaseline, elisaTotAverage,
                            by.y = "exp", by.x = "Baseline",all = T) %>%
        rename(BaseValues = meanValues) %>% select(-Blank,-BlankValues)
      
      elisaTot_base = merge(elisaTotAverage, elisaTot_base, 
                            by.x = c("exp","time"), by.y = c("Exp","time"),
                            all.x = T  )
      
      elisaResult$data = elisaTot
      
      if(length(elisaTot_base[,1]) != 0 && !is.null(elisaResult$Regression) ){
        
        elisamean = elisaTot_base %>%
          rename( MeanExperiment = meanValues,
                  MeanBaseline = BaseValues ) %>%
          dplyr::mutate(Quantification =  elisaResult$Regression$fun(MeanExperiment) ) %>%
          #MeanExperiment/MeanBaseline * 100) %>%
          rename(Experiment = exp,Time = time) 
        
        output$ELISAtables = renderDT(elisamean)
        
        elisaResult$dataFinal = elisamean
        
        output$ELISAplots = renderPlot(
          {
            elisamean %>%
              ggplot( aes(x = Time, y = Quantification,
                          fill= Experiment, group = Experiment ) )+
              geom_bar(position = "dodge",stat = "identity")+
              theme_bw()+
              labs(x = "Time", col = "Experiments",
                   y = "Average quantifications obtained\n from the lm ")
          }
        )
      }else{
        output$ELISAtables = renderDT(data.frame(Error = "No linear model!"))
      }
    }
  })
  
  # here the Exp boxes are updated every time a new experiment is added 
  observeEvent(FlagsELISA$EXPselected,{
    expToselect = FlagsELISA$EXPselected
    baselines =  FlagsELISA$BASEselected
    blanks = FlagsELISA$BLANCHEselected
    
    expToselect = expToselect[expToselect != ""]
    
    # baselines updating
    output$ElisaBaselineSelection <- renderUI({
      select_output_list <- lapply(expToselect[! expToselect %in% baselines],
                                   function(i) {
                                     if(length(input[[paste0("elisa_Exp",i)]])>0)
                                       expsel = input[[paste0("elisa_Exp",i)]]
                                     else 
                                       expsel = ""
                                     
                                     selectInput(inputId = paste0("elisa_Exp",i),
                                                 label = i,
                                                 choices = c("",baselines),
                                                 selected = expsel)
                                   })
      do.call(tagList, select_output_list)
    })
    # blanks updating
    output$ElisaBlankSelection <- renderUI({
      select_output_list <- lapply(unique(c(expToselect,baselines)), function(i) {
        
        if(length(input[[paste0("elisa_blExp",i)]])>0)
          expsel = input[[paste0("elisa_blExp",i)]]
        else 
          expsel = ""
        
        selectInput(inputId = paste0("elisa_blExp",i),
                    label = i,
                    choices = c("",blanks),
                    selected = expsel)
      })
      do.call(tagList, select_output_list)
    })
  })

  observeEvent(input$ELISA_standcurve,{
    elisaResult$data -> data
    if(input$ELISA_standcurve != ""){
      
      standcurve = data %>%
        filter(exp %in% input$ELISA_standcurve) %>%
        # group_by(exp,time) %>%
        # summarise(AverageMeasures = mean(values)) %>%
        # ungroup() %>%
        select(exp,time,values) %>%
        rename(Measures = values) %>%
        mutate(Concentrations = NA )
      
      # If nothing changes w..r.t. the already saved table then I keep the old one!
      if(!is.null(elisaResult$Tablestandcurve) && 
         all.equal(elisaResult$Tablestandcurve %>% select(-Concentrations),
                   standcurve  %>% select(-Concentrations) ))
      {
        standcurve =  elisaResult$Tablestandcurve
      }else{
        elisaResult$Tablestandcurve = standcurve
      }
      
      
      output$ELISA_Table_stdcurve <- DT::renderDataTable({
        DT::datatable( standcurve,
                       selection = 'none',
                       editable = list(target = "cell",
                                       disable = list(columns = 0:2) ),
                       #options = list(lengthChange = FALSE, autoWidth = TRUE),
                       rownames= FALSE
        )
      })
    } 
  })
  
  observeEvent(elisaResult$Tablestandcurve,{
    if(!is.null(elisaResult$Tablestandcurve) && dim(elisaResult$Tablestandcurve)[1]!=0){
      
      output$ELISA_Table_stdcurve <- DT::renderDataTable({
        DT::datatable( 
          elisaResult$Tablestandcurve,
          selection = 'none',
          editable = list(target = "cell",
                          disable = list(columns = 0:2) ),
          #options = list(lengthChange = FALSE, autoWidth = TRUE),
          rownames= FALSE
        )
      })
    } 
  })
  observeEvent(input$ELISA_Table_stdcurve_cell_edit, {
    cells = input$ELISA_Table_stdcurve_cell_edit
    cells$col = cells$col + 1
    elisaResult$Tablestandcurve <- editData( elisaResult$Tablestandcurve ,
                                             cells,
                                             'ELISA_Table_stdcurve')
  })
  observeEvent(input$ELISA_buttonRegression,{
    standcurve = elisaResult$Tablestandcurve
    standcurve$Concentrations = as.numeric(standcurve$Concentrations)
    if(!is.null(standcurve)){
      standcurve = standcurve %>% na.omit()
      
      
      regressionPlot = ggplot(standcurve,aes(Concentrations, Measures)) +
        geom_point() +
        theme_bw()
      
      if(input$regressionType == "Linear"){
        modelStancurve = lm(Measures~Concentrations, data = standcurve)
        
        infoLM = data.frame(x = min(standcurve$Concentrations) + c(1,1),
                            y = max(standcurve$Measures) + c(2,1.75),
                            text = c( paste0("y = ", signif(modelStancurve$coef[[2]], 5), "x + ",signif(modelStancurve$coef[[1]],5 )),
                                      paste0("Adj R2 = ",signif(summary(modelStancurve)$adj.r.squared, 5))) )
        
        fun = paste0("(x - ",modelStancurve$coef[[1]],")/", modelStancurve$coef[[2]])
        
        regressionPlot =  regressionPlot +
          geom_smooth(method='lm', col = "red") +
          geom_text(data= infoLM,
                    aes(x = x, y = y, label =text ),
                    vjust = "inward", hjust = "inward" )
        
      }
      else if(input$regressionType == "Quadratic")
      {
        #this is not implemented
        standcurve$Concentrations2 = standcurve$Concentrations^2
        modelStancurve = lm(Measures~Concentrations+Concentrations2, data = standcurve)
        
        infoLM = data.frame(x = min(standcurve$Concentrations) + c(1,1),
                            y = max(standcurve$Measures) + c(2,1.75),
                            text = c( paste0("y = ", signif(modelStancurve$coef[[3]], 5), "x^2 + ",
                                             signif(modelStancurve$coef[[2]], 5), "x + ",signif(modelStancurve$coef[[1]],5 )),
                                      paste0("Adj R2 = ",signif(summary(modelStancurve)$adj.r.squared, 5))) )
        
        fun = paste0(modelStancurve$coef[[3]],"*x^2 + ",modelStancurve$coef[[2]],"*x + ",modelStancurve$coef[[1]] )
        
        regressionPlot =  regressionPlot  +
          geom_point() +
          stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,col="red")+
          geom_text(data= infoLM,
                    aes(x = x, y = y, label =text ),
                    vjust = "inward", hjust = "inward" )
      }
      else if(input$regressionType == "Hyperbola"){
        
        outNLreg = tryCatch(
          {
            modelStancurve<-nls(
              Measures ~ a*Concentrations/(b+Concentrations), 
              data = standcurve, #%>% group_by(Concentrations) %>% summarise(Measures = mean(Measures)),
              start = list(a = 1,b = 1)
            )
          }, 
          error = function(e){
            return(e)
          })
        
        if(!is.null(outNLreg$mess)){
          modelStancurve = NULL
          regressionPlot = ggplot()+ geom_text(data = data.frame(x = 1,y =1,text = paste0("Error: ",outNLreg$mess)),
                                               aes(x,y,label = text),color = "red")
        }else{
          modelStancurve = outNLreg
          coef = modelStancurve$m$getPars()
          r2 = 1- sum(modelStancurve$m$resid()^2)/(sum(( mean(standcurve$Measures) - modelStancurve$m$predict() )^2))
          
          infoLM = data.frame(x = min(standcurve$Concentrations) + c(1,1),
                              y = max(standcurve$Measures) + c(2,1.75),
                              text = c( paste0("y = ", signif(coef["a"], 5), "x / ( ",
                                               signif(coef["b"], 5), " + x ) "),
                                        paste0("R2 = ",signif(r2, 5))) )
          
          dfHyperbola = data.frame(x = seq(min(standcurve$Concentrations),max(standcurve$Concentrations),length.out = 20)) %>%
            mutate(y = (coef["a"]*x/((coef["b"]+x)) ) )
          
          fun = paste0(coef["b"],"*x/(",coef["a"],"-x)")
          
          regressionPlot =  regressionPlot  +
            geom_point() +
            geom_line(data = dfHyperbola,aes(x = x,y = y),size = 1,col="red" )+
            geom_text(data= infoLM,
                      aes(x = x, y = y, label =text ),
                      vjust = "inward", hjust = "inward" )
        }
      }
      
      elisaResult$Regression = list(data = modelStancurve, plot = regressionPlot, fun = function(x){ eval( parse(text = fun ) ) } )
      
    }else{
      regressionPlot = ggplot()
    }
    output$ELISAregression <- renderPlot(regressionPlot)
  })
  
  # save everytime there is a change in the results
  # ELISAresultListen <- reactive({
  #   reactiveValuesToList(elisaResult)
  # })
  # observeEvent(ELISAresultListen(), {
  #   DataAnalysisModule$elisaResult = reactiveValuesToList(elisaResult)
  # })
  
  ### End ELISA analysis ####  
  #### ENDOCYTOSIS analysis ####
  observeEvent(input$NextEndocQuantif,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "tablesENDOC")
  })
  
  endocResult = reactiveValues(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    ENDOCcell_TIME = NULL,
    ENDOCcell_COLOR = NULL,
    ENDOCcell_EXP = NULL,
    MapBaseline = NULL,
    MapBlank = NULL)
  
  endocResult0 = list(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    ENDOCcell_TIME = NULL,
    ENDOCcell_COLOR = NULL,
    MapBaseline = NULL,
    MapBlank = NULL)
  
  ENDOCresultListen <- reactive({
    reactiveValuesToList(endocResult)
  })
  observeEvent(ENDOCresultListen(), {
    DataAnalysisModule$endocResult = reactiveValuesToList(endocResult)
    DataAnalysisModule$endocResult$Flags = reactiveValuesToList(FlagsENDOC)
  })
  
  ##
  FlagsENDOC <- reactiveValues(cellCoo = NULL,
                               AllExp = "",
                               BASEselected = "",
                               BLANCHEselected = "",
                               EXPselected = "",
                               EXPcol = NULL)
  
  left_data_endoc <- reactiveVal()
  right_data_endoc <- reactiveVal()
  
  observeEvent(input$LoadENDOC_Button,{
    alert$alertContext <- "ENDOC-reset"
    if( !is.null(endocResult$Initdata) ) {
      shinyalert(
        title = "Important message",
        text = "Do you want to update the ENDOC data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else loadExcelFileENDOC()
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "ENDOC-reset") {  
      resetPanel("ENDOC", flags = FlagsENDOC, result = endocResult)
      updateCheckboxGroupInput(session, "ENDOC_baselines", choices = list(), selected = character(0))
      updateCheckboxGroupInput(session, "ENDOC_blanks", choices = list(), selected = character(0))
  
      updateSelectizeInput(session, "ENDOCcell_EXP", choices = character(0), selected = character(0))
      updateSelectizeInput(session, "ENDOCcell_TIME", choices = character(0), selected = character(0))
      
      left_data_endoc <- NULL
      right_data_endoc <- NULL
      loadExcelFileENDOC()
    }
  })
  
  loadExcelFileENDOC <- function() {
    alert$alertContext <- ""
  
    for(nameList in names(endocResult0)) 
      endocResult[[nameList]] <- endocResult0[[nameList]]
  
      mess = readfile(
        filename = input$ENDOCImport$datapath,
        isFileUploaded = !is.null(input$ENDOCImport) && file.exists(input$ENDOCImport$datapath),
        type = "Excel",
        allDouble = T,
        colname = F,
        colors = T,
      )
      
      if(setequal(names(mess), c("message", "call"))) {
        showAlert("Error", mess[["message"]], "error", 5000)
      } else {
        endocResult$Initdata = mess$x
        FlagsENDOC$EXPcol = mess$fill
        endocResult$ENDOCcell_COLOR = mess$SNtable
        endocResult$ENDOCcell_EXP <- matrix("", nrow = nrow(endocResult$ENDOCcell_COLOR), ncol = ncol(endocResult$ENDOCcell_COLOR))
        
        removeModal()
        showAlert("Success", "The Excel has been uploaded  with success", "success", 2000)
      }
  }
  
  observe({
    if (!is.null(endocResult$Initdata) && is.null(endocResult$TablePlot)) {
      tableExcelColored(session = session,
                        Result = endocResult, 
                        FlagsExp = FlagsENDOC,
                        type = "Initialize")
    }
  })
    
  observe({
    color_codes <- FlagsENDOC$EXPcol
    color_names <- names(FlagsENDOC$EXPcol)
    
    valid_colors <- color_codes != "white"
    color_codes <- color_codes[valid_colors]
    color_names <- color_names[valid_colors]
    
    mid_point <- ceiling(length(color_codes) / 2)
    left_colors <- color_codes[1:mid_point]
    right_colors <- color_codes[(mid_point+1):length(color_codes)]
    
    left_formatted_data <- get_formatted_data(left_colors, color_names[1:mid_point], endocResult, endocResult$ENDOCcell_TIME,"ENDOC")
    right_formatted_data <- get_formatted_data(right_colors, color_names[(mid_point+1):length(color_codes)], endocResult, endocResult$ENDOCcell_TIME, "ENDOC")
  
    left_data_endoc(left_formatted_data)
    right_data_endoc(right_formatted_data)
    
    output$leftTableEndoc <- renderDataTable(
      left_data_endoc(), 
      escape = FALSE, 
      editable = list(target = "cell", disable = list(columns = 0:3)),
      options = list(
        dom = 't',
        paging = FALSE,
        info = FALSE,
        searching = FALSE, 
        columnDefs = list(
          list(targets = 0, visible = FALSE),
          list(targets = 1, visible = FALSE),
          list(width = '10px', targets = 2),
          list(width = '100px', targets = 3),
          list(width = '100px', targets = 4),
          list(width = '100px', targets = 5),
          list(className = 'dt-head-left dt-body-left', targets = 1)
        )
      )
    )
    
    output$rightTableEndoc <- renderDataTable(
      right_data_endoc(), 
      escape = FALSE, 
      editable = list(target = "cell", disable = list(columns = 0:3)),
      options = list(
        dom = 't',
        paging = FALSE,
        info = FALSE,
        searching = FALSE,
        editable = TRUE,
        columnDefs = list(
          list(targets = 0, visible = FALSE),
          list(targets = 1, visible = FALSE),
          list(width = '10px', targets = 2),
          list(width = '100px', targets = 3),
          list(width = '100px', targets = 4),
          list(width = '100px', targets = 5),
          list(className = 'dt-head-left dt-body-left', targets = 1)
        )
      )
    )
  })

  observeEvent(input$leftTableEndoc_cell_edit, {
    info <- input$leftTableEndoc_cell_edit
    data <- left_data_endoc() 
    updatedText <- updateTable("left", "ENDOC", info, data, endocResult, FlagsENDOC)
    
    output$ENDOCSelectedValues <- renderText(updatedText)  
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent(input$rightTableEndoc_cell_edit, {
    info <- input$rightTableEndoc_cell_edit
    data <- right_data_endoc() 
    updatedText <- updateTable("right", "ENDOC", info, data, endocResult, FlagsENDOC)
    
    output$ENDOCSelectedValues <- renderText(updatedText)  
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent(input$ENDOCmatrix_cell_clicked, {
    req(input$ENDOCmatrix_cell_clicked)  
    
    cellSelected = as.numeric(input$ENDOCmatrix_cell_clicked)
    FlagsENDOC$cellCoo = cellCoo = c(cellSelected[1], cellSelected[2] + 1)
    
    allTime <- unique(na.omit(c(endocResult$ENDOCcell_TIME)))  
    selectedTime <- ifelse(is.null(endocResult$ENDOCcell_TIME[cellCoo[1], cellCoo[2]]), "", endocResult$ENDOCcell_TIME[cellCoo[1], cellCoo[2]])
    
    updateSelectizeInput(inputId = "ENDOCcell_TIME", 
                         choices = allTime,
                         selected = selectedTime)
    
    allConditions <- unique(na.omit(c(endocResult$ENDOCcell_EXP)))  
    selectedCondition <- ifelse(is.null(endocResult$ENDOCcell_EXP[cellCoo[1], cellCoo[2]]), "", endocResult$ENDOCcell_EXP[cellCoo[1], cellCoo[2]])
    
    updateSelectizeInput(inputId = "ENDOCcell_EXP",
                         choices = allConditions,
                         selected = selectedCondition)
  })
  
  observeEvent(input$ENDOCcell_EXP, {
    if (!is.null(endocResult$ENDOCcell_COLOR) && !is.null(FlagsENDOC$cellCoo) && !anyNA(FlagsENDOC$cellCoo)) {
      ENDOCtb = endocResult$TablePlot
      cellCoo = FlagsENDOC$cellCoo
      
      value.bef = endocResult$ENDOCcell_COLOR[cellCoo[1], cellCoo[2]] 
      value.now = input$ENDOCcell_EXP
      
      if (value.now != "" && value.now != value.bef) {
        currentValues <- endocResult$Initdata[cellCoo[1], cellCoo[2]]
        
        endocResult$ENDOCcell_COLOR[cellCoo[1], cellCoo[2]] = value.now
        endocResult$ENDOCcell_EXP[cellCoo[1], cellCoo[2]] = value.now
        ENDOCtb$x$data[cellCoo[1], paste0("Col", cellCoo[2])] = value.now
        
        if (!input$ENDOCcell_EXP %in% FlagsENDOC$AllExp) {
          exp = unique(c(FlagsENDOC$AllExp, input$ENDOCcell_EXP))
          FlagsENDOC$AllExp = exp
        }
        
        tableExcelColored(session = session,
                           Result = endocResult, 
                           FlagsExp = FlagsENDOC,
                           type = "Update")
        
        output$ENDOCSelectedValues <- renderText(paste("Updated value", paste(currentValues), ": experimental condition ", value.now))
        output$ENDOCmatrix <- renderDataTable({endocResult$TablePlot})
      }
    } else return()
  }, ignoreInit = TRUE)
  
  observeEvent(input$ENDOCcell_TIME, {
    if (!is.null(endocResult$ENDOCcell_TIME) && !is.null(FlagsENDOC$cellCoo) && !anyNA(FlagsENDOC$cellCoo)) {
      ENDOCtb = endocResult$TablePlot
      cellCoo = FlagsENDOC$cellCoo
      
      value.bef = endocResult$ENDOCcell_TIME[cellCoo[1], cellCoo[2]] 
      value.now = input$ENDOCcell_TIME
    
      if (value.now != "" && value.now != value.bef) {
        currentValues <- endocResult$Initdata[cellCoo[1], cellCoo[2]]
        
        endocResult$ENDOCcell_TIME[cellCoo[1], cellCoo[2]] = value.now
        tableExcelColored(session = session,
                          Result = endocResult, 
                          FlagsExp = FlagsENDOC,
                          type = "Update",
                          )
        
        output$ENDOCSelectedValues <- renderText(paste("Updated value", paste(currentValues), ": time ", value.now))
        output$ENDOCmatrix <- renderDataTable({endocResult$TablePlot})
      }
    } else return()
  }, ignoreInit = TRUE)
  
  ## update Baselines checkBox
  observeEvent(c(FlagsENDOC$AllExp,FlagsENDOC$BLANCHEselected),{
    if(length(FlagsENDOC$AllExp) > 1){
      exp = FlagsENDOC$AllExp
      exp = exp[exp != ""]
      
      if(!( length(FlagsENDOC$BLANCHEselected) == 1 && FlagsENDOC$BLANCHEselected == "") )
        exp = exp[!exp %in% FlagsENDOC$BLANCHEselected]
      
      exp_selec = input$ENDOC_baselines
      
      updateCheckboxGroupInput(session,"ENDOC_baselines",
                               choices = exp,
                               selected = exp_selec )
      
      FlagsENDOC$EXPselected = exp
    }
  })
  
  observeEvent(c(FlagsENDOC$AllExp,FlagsENDOC$BASEselected),{
    if(length(FlagsENDOC$AllExp) > 1){
      exp = FlagsENDOC$AllExp
      exp = exp[exp != ""]
      
      if(! (length(FlagsENDOC$BASEselected) == 1 && FlagsENDOC$BASEselected == "") )
        exp = exp[!exp %in% FlagsENDOC$BASEselected]
      
      exp_selec = input$ENDOC_blanks
      
      updateCheckboxGroupInput(session,"ENDOC_blanks",
                               choices = exp,
                               selected = exp_selec )
      
      FlagsENDOC$EXPselected = exp
    }
  })
  
  ## select the baselines and blank
  observeEvent(input$ENDOC_baselines,{
    FlagsENDOC$BASEselected = input$ENDOC_baselines
    FlagsENDOC$EXPselected = FlagsENDOC$AllExp[! FlagsENDOC$AllExp %in% c(FlagsENDOC$BASEselected,FlagsENDOC$BLANCHEselected)]
  },ignoreNULL = F)
  observeEvent(input$ENDOC_blanks,{
    FlagsENDOC$BLANCHEselected = input$ENDOC_blanks
    FlagsENDOC$EXPselected = FlagsENDOC$AllExp[! FlagsENDOC$AllExp %in% c(FlagsENDOC$BASEselected,FlagsENDOC$BLANCHEselected)]
  },ignoreNULL = F)
  
  toListen_endoc <- reactive({
    exp = FlagsENDOC$EXPselected
    exp = exp[exp != ""]
    if(length(exp) > 0 )
    {
      Input_baselEXP = lapply(exp,
                              function(i) input[[paste0("Exp",i)]])
      Input_blEXP = lapply(unique(exp,FlagsENDOC$BASELINEselected),
                           function(i) input[[paste0("blExp",i)]] )
      InputEXP = c(Input_baselEXP,Input_blEXP)
      
      which(sapply(InputEXP, is.null)) -> indexesEXPnull
      if(length(indexesEXPnull) > 0 )
        listReturn = InputEXP[-indexesEXPnull]
      else
        listReturn = InputEXP
    }else{
      listReturn = list()
    }
    
    if(length(listReturn) == 0){
      return(list("Nothing",endocResult$ENDOCcell_TIME,endocResult$ENDOCcell_COLOR))
    }else{
      return(c(listReturn,list(endocResult$ENDOCcell_TIME,endocResult$ENDOCcell_COLOR)) )
    }
  })
  
  observeEvent(toListen_endoc(),{
    baselines = FlagsENDOC$BASEselected
    baselines = baselines[baselines != ""]
    
    if(toListen_endoc()[[1]] != "Nothing"){
      exp = FlagsENDOC$EXPselected
      exp = exp[exp != ""]
      expNotBlank = unique(c(exp,baselines))
      
      MapBaseline = do.call(rbind,
                            lapply(exp,function(i){
                              if( length(input[[paste0("Exp",i)]]) > 0 && input[[paste0("Exp",i)]] != ""){
                                data.frame(Exp = i, Baseline = input[[paste0("Exp",i)]])
                              }else{
                                data.frame(Exp = i, Baseline = NA)
                              }
                            })
      ) %>% na.omit()
      
      MapBlank = do.call(rbind,
                         lapply(expNotBlank,
                                function(i){
                                  if( length(input[[paste0("blExp",i)]]) > 0 && input[[paste0("blExp",i)]] != ""){
                                    data.frame(Exp = i, Blank = input[[paste0("blExp",i)]])
                                  }else{
                                    data.frame(Exp = i, Blank = NA)
                                  }
                                })
      ) %>% na.omit()
      
      if(dim(MapBaseline)[1]!=0 && dim(MapBlank)[1]!=0 ){
        
        endocResult$MapBaseline = MapBaseline
        endocResult$MapBlank = MapBlank
        
        mat = as.matrix(endocResult$Initdata)
        endocV = expand.grid(seq_len(nrow(mat)), seq_len(ncol(mat))) %>%
          rowwise() %>%
          mutate(values = mat[Var1, Var2])
        matTime =  as.matrix(endocResult$ENDOCcell_TIME)
        endocT = expand.grid(seq_len(nrow(matTime)), seq_len(ncol(matTime))) %>%
          rowwise() %>%
          mutate(time = matTime[Var1, Var2])
        matExp =  as.matrix(endocResult$ENDOCcell_COLOR)
        endocE = expand.grid(seq_len(nrow(matExp)), seq_len(ncol(matExp))) %>%
          rowwise() %>%
          mutate(exp = matExp[Var1, Var2])
        endocTot = merge(endocV,merge(endocT,endocE)) %>%
          filter(exp != "")
        
        endocTotAverage = endocTot %>%
          mutate(time = ifelse(exp %in% MapBlank$Blank, 0, time)) %>%
          group_by(time, exp) %>%
          summarize(meanValues = mean(values))
        
        # merging exp with blank for the substraction
        
        endocTot_bl = right_join( endocTotAverage,MapBlank, 
                                  by= c("exp"= "Blank") )%>%
          rename(BlankValues = meanValues, Blank =  exp, exp = Exp )
        
        endocTotAverage = merge( endocTotAverage %>% filter(! exp %in%endocTot_bl$Blank ),
                                 endocTot_bl %>% ungroup() %>%
                                   dplyr::select(-time), all.x = T, by = "exp") %>%
          rename(Exp = exp) 
        endocTotAverage = endocTotAverage %>% mutate(meanValues = meanValues - BlankValues )
        
        # merging exp with baseline
        endocTot_base = merge(MapBaseline, endocTotAverage,
                              by.y = "Exp", by.x = "Baseline") %>%
          rename(BaseValues = meanValues) %>% 
          select(-Blank,-BlankValues)
        
        if(length(unique(endocTot_base$time)) == 1){ 
          # if there is only one point then the baseline is used to normalize every times
          endocTot_base = endocTot_base %>% select(-time)
        }
        endocTot_base = merge(endocTotAverage,
                              endocTot_base )
        # by.x = c("exp","time"),
        # by.y = c("Exp","time")
        #)
        endocResult$data = endocTot
        
        if(length(endocTotAverage[,1]) != 0 ){
          endocmean = endocTot_base %>%
            rename( MeanExperiment = meanValues,
                    MeanBaseline = BaseValues ) %>%
            dplyr::mutate(Quantification = MeanExperiment/MeanBaseline * 100) %>%
            rename(Experiment = Exp,Time = time) 
          
          output$ENDOCtables = renderDT(endocmean)
          
          endocResult$dataFinal = endocmean
          
          output$ENDOCplots = renderPlot(
            {
              pl1 = endocmean %>%
                ggplot( aes(x = Time, y = MeanBaseline,
                            col= Baseline, group = Experiment ) )+
                geom_point( )+
                geom_line()+
                theme_bw()+
                labs(x = "Time", col = "Baselines selected",
                     y = "Average Baseline\n quantifications")
              
              pl2 = endocmean %>%
                mutate(ExperimentBaseline = paste0(Experiment,"/",Baseline)) %>%
                ggplot( aes(x = Time, y = Quantification,
                            col= ExperimentBaseline, group = Experiment ) )+
                geom_point()+
                geom_line()+
                theme_bw()+
                labs(x = "Time", col = "Ratio\n Experiment / Baseline",
                     y = "Ratio of the average quantifications\n (as %)")
              
              pl2/pl1
            }
          )
        }else{
          output$ENDOCtables = renderDT(data.frame(Error = "No baseline is associated with the experiment replicants!"))
        }        
      }
    }
  })
  
  # here the Exp boxes are updated every time a new experiment is added 
  observeEvent(FlagsENDOC$EXPselected,{
    expToselect = FlagsENDOC$EXPselected
    baselines =  FlagsENDOC$BASEselected
    blanks = FlagsENDOC$BLANCHEselected
    
    expToselect = expToselect[expToselect != ""]
    
    # baselines updating
    output$EndocBaselineSelection <- renderUI({
      select_output_list <- lapply(expToselect[! expToselect %in% baselines],
                                   function(i) {
                                     if(length(input[[paste0("Exp",i)]])>0)
                                       expsel = input[[paste0("Exp",i)]]
                                     else 
                                       expsel = ""
                                     
                                     selectInput(inputId = paste0("Exp",i),
                                                 label = i,
                                                 choices = c("",baselines),
                                                 selected = expsel)
                                   })
      do.call(tagList, select_output_list)
    })
    
    # blanks updating
    output$EndocBlankSelection <- renderUI({
      select_output_list <- lapply(unique(c(expToselect,baselines)), function(i) {
        
        if(length(input[[paste0("blExp",i)]])>0)
          expsel = input[[paste0("blExp",i)]]
        else 
          expsel = ""
        
        selectInput(inputId = paste0("blExp",i),
                    label = i,
                    choices = c("",blanks),
                    selected = expsel)
      })
      do.call(tagList, select_output_list)
    })
  })
  
  observeEvent(endocResult$TablePlot, {
    ENDOCtb = endocResult$TablePlot
    output$ENDOCmatrix <- renderDT(
      ENDOCtb,
      server = FALSE
    )
    
    ##### Plot the values selected!
    matTime =  as.matrix(endocResult$ENDOCcell_TIME)
    matExp =  as.matrix(endocResult$ENDOCcell_COLOR)
    
    if( !( all(matTime == "")  || all(matExp == "") ) ){
      mat = as.matrix(endocResult$Initdata)
      endocV = expand.grid(seq_len(nrow(mat)), seq_len(ncol(mat))) %>%
        rowwise() %>%
        mutate(values = mat[Var1, Var2])
      endocT = expand.grid(seq_len(nrow(matTime)), seq_len(ncol(matTime))) %>%
        rowwise() %>%
        mutate(time = matTime[Var1, Var2])
      endocE = expand.grid(seq_len(nrow(matExp)), seq_len(ncol(matExp))) %>%
        rowwise() %>%
        mutate(exp = matExp[Var1, Var2])
      endocTot = merge(endocV,merge(endocT,endocE)) %>%
        na.omit() %>%
        filter(time != "",  exp != "") 
      
      endocResult$data = endocTot
      
      output$ENDOCinitplots <- renderPlot(
        ggplot(endocTot, aes(x = time, y=values, col = exp),alpha = 1.4) +
          #geom_boxplot(aes(fill= exp, group = time),alpha = 0.4) +
          geom_point(aes(group = exp)) +
          scale_color_manual(values = FlagsENDOC$EXPcol) + 
          #scale_fill_manual(values = FlagsENDOC$EXPcol) + 
          theme_bw()+
          labs(x = "Times", y = "Values", col = "Exp",fill = "Exp")+
          theme(legend.position = c(0, 1), 
                legend.justification = c(0, 1),
                legend.direction = "vertical",
                legend.background = element_rect(linewidth = 0.5,
                                                 linetype="solid",
                                                 colour ="black"))
      )
    }
  })

  output$downloadENDOCAnalysis <- downloadHandler(
    filename = function() {
      paste('ENDOCanalysis-', Sys.Date(), '.zip', sep='')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      nomeRDS <- paste0("ENDOC_analysis-", Sys.Date(), ".rds")
      nomeXLSX <- paste0("ENDOC_analysis-", Sys.Date(), ".xlsx")
      
      tempRdsPath <- file.path(tempDir, nomeRDS)
      tempXlsxPath <- file.path(tempDir, nomeXLSX)
      
      results <- endocResult  # Assicurati di accedere a endocResult in un contesto reattivo, se necessario
      saveRDS(results, file = tempRdsPath)
      
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "ENDOC")
      
      zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
      
      } 
  )
  
  ### End ENDOC analysis ####
  
  #### CITOXICITY analysis ####
  
  ### End CITOXICITY analysis ####
  
  #### FACS analysis ####
  facsResult = reactiveValues(
    Initdata= NULL,
    data = NULL,
    dataFinal = NULL,
    depth = NULL,
    depthCount = NULL,
    name = NULL,
    originalName = NULL,
    statistics = NULL,
    cells = NULL
  )
  
  facsResult0 = reactiveValues(
    Initdata= NULL,
    data = NULL,
    dataFinal = NULL,
    depth = NULL,
    depthCount = NULL,
    name = NULL,
    originalName = NULL,
    statistics = NULL,
    cells = NULL
  )
  
  FACSresultListen <- reactive({
    reactiveValuesToList(facsResult)
  })
  observeEvent(FACSresultListen(), {
    DataAnalysisModule$facsResult = reactiveValuesToList(facsResult)
    DataAnalysisModule$facsResult$Flags = reactiveValuesToList(FlagsFACS)
  })
  
  FlagsFACS <- reactiveValues(
    actualLevel = NULL,
    allLevel = NULL,
    actualPath = NULL
  )
  
  observeEvent(input$LoadFACS_Button,{
    alert$alertContext <- "FACS-reset"
    if(!is.null(facsResult$name) ) {
      shinyalert(
        title = "Important message",
        text = "Do you want to update the FACS data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else loadExcelFileFACS()
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "FACS-reset") {  
      resetPanel("FACS", flags = FlagsFACS, result = facsResult)
      updateSelectizeInput(session, "FACScell", choices = character(0), selected = character(0))
      
      loadExcelFileFACS()
    }
  })
  
  loadExcelFileFACS <- function() {
    alert$alertContext <- ""
    
    for(nameList in names(facsResult0)) 
      facsResult[[nameList]] <- facsResult0[[nameList]]
    
    mess = readfile(
      filename = input$FACSImport$datapath,
      isFileUploaded = !is.null(input$FACSImport) && file.exists(input$FACSImport$datapath),
      type = "Excel",
      allDouble = TRUE,
      colname = FALSE,
      colors = FALSE,
    )
    
    if (setequal(names(mess), c("message", "call"))) {
      showAlert("Error", mess[["message"]], "error", 5000)
    } else {
      if (nrow(mess) > 1) {
        data <- mess[-1, , drop = FALSE]  
        facsResult$Initdata <- data
        
        facsResult$depth <- vector("list", nrow(data))
        facsResult$depthCount <- numeric(nrow(data))
        facsResult$name <- vector("list", nrow(data))
        facsResult$originalName <- vector("list", nrow(data))
        facsResult$statistics <- vector("list", nrow(data))
        facsResult$cells <- vector("list", nrow(data))
        
        for (i in 1:nrow(data)) {
          x <- data[i, 1]
          
          if (is.na(x)) {
            facsResult$depth[i] <- NA
            facsResult$depthCount[i] <- 0
            facsResult$name[i] <- data[i, 2]
            facsResult$statistics[i] <- data[i, 3]
            facsResult$cells[i] <- data[i, 4]
          } else {
            x_cleaned <- gsub(" ", "", x)  
            facsResult$depth[i] <- x 
            facsResult$depthCount[i] <- str_count(x_cleaned, fixed(">"))
            facsResult$name[i] <- data[i, 2]
            facsResult$statistics[i] <- data[i, 3]
            facsResult$cells[i] <- data[i, 4]
          }
        }
        
        removeModal()
        
        maxDepth <- max(facsResult$depthCount, na.rm = TRUE)
        
        # Primo: cambia la tab
        updateTabsetPanel(session, "SideTabs", selected = "tablesFACS")
        
        # Secondo: Ascolta il cambio della tab e poi esegue le funzioni
        observeEvent(input$SideTabs, {
          if (input$SideTabs == "tablesFACS") {
            updateSelectizeUI(maxDepth)
            FlagsFACS$actualLevel <- 0
            showAlert("Success", "The Excel has been uploaded with success", "success", 2000)
          }
        }, ignoreInit = TRUE)  # Assicurati di non eseguire al primo caricamento
      }
    }
  }
  
  updateSelectizeUI <- function(maxDepth) {
    output$dynamicSelectize <- renderUI({
      rowContent <- fluidRow(
        lapply(1:maxDepth, function(i) {
          column(
            2, offset = 1,
            tags$div(style = "display: none;", id = paste("div_FACScell", i, sep = ""),
                     selectizeInput(
                       inputId = paste("FACScell", i, sep = "_"),
                       label = paste("Gate", i),
                       choices = c(),
                       options = list(placeholder = 'Select the next gate', create = TRUE)
                     )
            )
          )
        })
      )
      rowContent
    })
  }
  
  observeEvent(FlagsFACS$actualLevel, {
  currentInputId <- paste("FACScell", FlagsFACS$actualLevel, sep = "_")
  
  if (FlagsFACS$actualLevel == 0) {
    valid_indices <- facsResult$depthCount == FlagsFACS$actualLevel
    filtered_cells <- list()
    filtered_names <- list()
    
    if (any(valid_indices)) {
      for (i in which(valid_indices)) {
        filtered_cells <- c(filtered_cells, facsResult$cells[i])
        filtered_names <- c(filtered_names, facsResult$name[i])
      }
    }
    
    data_for_table <- data.frame(
      Name = unlist(filtered_names),
      Start = unlist(filtered_cells),
      stringsAsFactors = FALSE
    )
    
    facsResult$originalName <- unlist(filtered_names)
    FlagsFACS$data <- data_for_table  
    FlagsFACS$actualPath <- ""
    
    output$FACSmatrix <- renderDT({
      datatable(data_for_table, options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(title = "", targets = 2)  
        )
      ))
    })
    
    data_for_name_update <- data.frame(
      Name = unlist(filtered_names),
      New_name = rep("-", length(filtered_names)),  # Colonna con trattini iniziali
      stringsAsFactors = FALSE
    )
    
    output$FACSnameUpdate <- renderDT({
      datatable(data_for_name_update, options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 1, width = '50%'),  # Imposta la larghezza della prima colonna
          list(targets = 2, width = '50%')
        )
      ), editable = list(target = 'cell', columns = 2))  # Permetti l'editing della seconda colonna
    })
  } else {
    selected_item <- input[[currentInputId]]
    path_components <- strsplit(FlagsFACS$actualPath, "/")[[1]]
    if (length(path_components) >= FlagsFACS$actualLevel) {
      path_components <- path_components[1:(FlagsFACS$actualLevel - 1)]
    }
    path_components <- c(path_components, selected_item)
    FlagsFACS$actualPath <- paste(path_components, collapse = "/")
    
    if (ncol(FlagsFACS$data) > length(path_components) + 1) {  
      FlagsFACS$data <- FlagsFACS$data[, c(1, (2:length(path_components) + 1)), drop = FALSE]
    }
    
    new_column_name <- tail(path_components, n = 1)  
    new_column_data <- vector("numeric", length = nrow(FlagsFACS$data))
    
    for (i in seq_len(nrow(FlagsFACS$data))) {
      new_path <- paste0(FlagsFACS$data$Name[i], "/", FlagsFACS$actualPath)
      index <- which(facsResult$name == new_path)
      
      if (length(index) == 1) {
        new_column_data[i] <- sprintf("%.2f%%", facsResult$statistics[index])
      } else {
        new_column_data[i] <- NA
      }
    }
    
    maxLevels <- max(facsResult$depthCount, na.rm = TRUE)
    if (FlagsFACS$actualLevel < maxLevels) {
      lapply((FlagsFACS$actualLevel + 1):maxLevels, function(level) {
        updateSelectizeInput(session, paste("FACScell", level, sep = "_"), choices = list())
      })
    }
    
    FlagsFACS$data[[new_column_name]] <- new_column_data
    
    output$FACSmatrix <- renderDT({
      datatable(FlagsFACS$data, options = list(autoWidth = TRUE))
    })
  }
  loadDrop() 
  }, ignoreInit = TRUE)
  
  observe({
    FlagsFACS$actualPath
    choices <- colnames(FlagsFACS$data)
    choices <- choices[choices != "Name"]
    
    if (length(choices) > 1) {
      choices <- choices[-length(choices)]
    }
    
    updateSelectizeInput(session = session,
                         inputId = "selectBaseGate",
                         choices = choices
    )
  })
  
  escapeRegex <- function(string) {
    gsub("([\\\\^$.*+?()[{\\]|-])", "\\\\\\1", string)
  }
  
  loadDrop <- function() {
    targetLevel <- FlagsFACS$actualLevel + 1
    currentPath <- FlagsFACS$actualPath
    
    escapedPath <- escapeRegex(currentPath)
    regex_path <- paste0(".*", escapedPath, "/[^/]+$")
    valid_indices <- facsResult$depthCount == targetLevel & grepl(regex_path, facsResult$name)
    
    valid_names <- facsResult$name[valid_indices]
    valid_names <- as.character(valid_names)

    if (length(valid_names) > 0) {
      short_names <- sapply(strsplit(valid_names, "/", fixed = TRUE), function(x) tail(x, 1))
    } else {
      short_names <- "no valid names found"
    }
    
    nextInputId <- paste("FACScell", targetLevel, sep = "_")
    nextDivId <- paste("div_FACScell", targetLevel, sep = "")
    
    if (length(short_names) == 0) {
      updateSelectInput(session, nextInputId, choices = list("No choices available" = ""), selected = "")
    } else {
      updateSelectInput(session, nextInputId, choices = setNames(short_names, short_names), selected = character(0))
    }
    
    shinyjs::runjs(paste0('setTimeout(function() { $("#', nextDivId, '").css("display", "block"); }, 200);'))
  }
  
  
  observe({
    if (is.null(facsResult$depthCount) || length(facsResult$depthCount) == 0) {
      numLevels <- 0
    } else {
      numLevels <- max(facsResult$depthCount, na.rm = TRUE)
    }
    
    if (numLevels > 0) {
      if (exists("levelObservers", envir = .GlobalEnv)) {
        lapply(levelObservers, function(observer) {
          observer$destroy()
        })
      }
      
      levelObservers <<- lapply(1:numLevels, function(level) {
        observeEvent(input[[paste("FACScell", level, sep = "_")]], {
          selectedValue <- input[[paste("FACScell", level, sep = "_")]]
          if (!is.null(selectedValue) && nzchar(selectedValue) && selectedValue != "no valid names found") {
            FlagsFACS$actualLevel <- 0
            FlagsFACS$actualLevel <- level
          }
        }, ignoreInit = TRUE)
      })
    } else {
      if (exists("levelObservers", envir = .GlobalEnv)) {
        lapply(levelObservers, function(observer) {
          observer$destroy()
        })
      }
      levelObservers <<- list()  
    }
  })
  
  convert_percent_to_final <- function(data_row, selected_gate) {
    last_percentage_name <- paste(names(data_row)[length(data_row)], selected_gate, sep = "/")
    
    if (selected_gate == "Start") {
      initial_value <- as.numeric(data_row[2])
      percentages <- as.numeric(sub("%", "", data_row[3:length(data_row)])) / 100
      
      final_percentage <- (Reduce(`*`, percentages) * 100 / initial_value) * 100
    } else {
      selected_gate_index <- which(names(data_row) == selected_gate)
      
      if (selected_gate_index == length(data_row) - 1) {
        final_percentage <- as.numeric(sub("%", "", data_row[selected_gate_index + 1]))
      } else {
        percentages <- as.numeric(sub("%", "", data_row[(selected_gate_index + 2):length(data_row)])) / 100
        
        if (length(percentages) == 0) {
          final_percentage <- as.numeric(sub("%", "", data_row[selected_gate_index + 1]))
        } else {
          selected_percentage <- as.numeric(sub("%", "", data_row[selected_gate_index + 1])) / 100
          final_percentage <- selected_percentage * Reduce(`*`, percentages) * 100
        }
      }
    }
    
    data.frame(Name = data_row[1], FinalPercent = sprintf("%.2f%%", final_percentage)) %>%
      setNames(c("Name", last_percentage_name))
  }
  
  observeEvent(input$SaveFACSanalysis, {
    selected_gate <- input$selectBaseGate
    selected_gate_index <- which(names(FlagsFACS$data) == selected_gate)
    
    if (length(selected_gate_index) == 0) {
      showNotification("Selected column not found in the dataset", type = "error")
      return()
    }
    
    current_data <- FlagsFACS$data[, c(1, 2, selected_gate_index:ncol(FlagsFACS$data)), drop = FALSE]
    
    processed_data <- lapply(split(current_data, seq(nrow(current_data))), function(row) {
      convert_percent_to_final(row, selected_gate)
    })
    
    new_data <- do.call(rbind, processed_data)
    
    if (is.null(facsResult$data)) {
      facsResult$data <- new_data
    } else {
      new_column_name <- names(new_data)[2]
      if (names(facsResult$data)[ncol(facsResult$data)] != new_column_name) {
        facsResult$data[new_column_name] <- new_data[, 2]
      }
    }
    
    output$FACSresult <- renderDT({
      datatable(facsResult$data, 
                options = list(
                  autoWidth = TRUE,
                  columnDefs = list(
                    list(visible = FALSE, targets = 0)  
                  )
                )
      )
    })
    facsResult$dataFinal <- facsResult$data
  })
  
  output$downloadFACSanalysis <- downloadHandler(
    filename = function() {
      paste('FACSanalysis-', Sys.Date(), '.zip', sep='')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      
      nomeRDS <- paste0("FACSquant_analysis-", Sys.Date(), ".rds")
      nomeXLSX <- paste0("FACSquant_analysis-", Sys.Date(), ".xlsx")
      
      tempRdsPath <- file.path(tempDir, nomeRDS)
      tempXlsxPath <- file.path(tempDir, nomeXLSX)
      results <- DataAnalysisModule$facsResult
      saveRDS(results, file = tempRdsPath)
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "FACS", PanelStructures)
      
      zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
    },
  )
  
  observeEvent(input$FACSnameUpdate_cell_edit, {
    info <- input$FACSnameUpdate_cell_edit
    
    row <- info$row
    col <- info$col
    new_value <- info$value
    
    if (col == 2) {  
      associated_name <- FlagsFACS$data[row, "Name"]
      
      index <- which(facsResult$name == associated_name)
      
      if (length(index) > 0) {
        if (new_value == "") {
          original_name <- as.character(facsResult$originalName[row])
          for (i in seq_along(facsResult$name)) {
            name_parts <- strsplit(as.character(facsResult$name[i]), "/", fixed = TRUE)[[1]]
            if (name_parts[1] == associated_name) {
              name_parts[1] <- original_name
              facsResult$name[i] <- paste(name_parts, collapse = "/")
            }
          }
        } else {
          for (i in seq_along(facsResult$name)) {
            if (startsWith(as.character(facsResult$name[i]), as.character(associated_name))) {
              suffix <- substr(facsResult$name[i], nchar(associated_name) + 1, nchar(facsResult$name[i]))
              facsResult$name[i] <- paste0(new_value, suffix)
            }
          }
        }
        print(facsResult$name)
      }
    }
  })
  
  observeEvent(facsResult$name, {
    valid_indices <- facsResult$depthCount == FlagsFACS$actualLevel
    filtered_names <- facsResult$name[valid_indices]
    
    FlagsFACS$data$Name <- sapply(filtered_names, function(name) strsplit(name, "/")[[1]][1])
    
    proxy <- dataTableProxy('FACSmatrix')
    replaceData(proxy, FlagsFACS$data, resetPaging = FALSE)
  }, ignoreInit = TRUE)
  
  
  ### End FACS analysis ####
  
  #start statistics
  DataStatisticModule = reactiveValues(WB = list(),
                                       PRCC = list(),
                                       ELISA = list(),
                                       ENDOC = list(),
                                       CYTOTOX = list(),
                                       Flag = F)
  
  DataStatisticresultListen <- reactive({
    reactiveValuesToList(DataStatisticModule)
  })
  
  observeEvent(DataStatisticresultListen(),{
    
    if(DataStatisticModule$Flag){
      AnalysisNames = names(DataStatisticModule)[names(DataStatisticModule) != "Flag"]
      Analysis = rep(F,length(AnalysisNames) )
      names(Analysis) = AnalysisNames
      for(j in AnalysisNames)
        Analysis[j] = all(sapply(DataStatisticModule[[j]], is.null))
      
      AnalysisNames = AnalysisNames[!Analysis]
      
      updateSelectizeInput(inputId = "StatAnalysis",
                           choices = c("",AnalysisNames),
                           selected = "")
      
      DataStatisticModule$Flag = F
    }
  })
  
  observeEvent(input$loadStatAnalysis_file_Button,{
    manageSpinner(TRUE)
    
    result <- readfile(
      filename = input$loadStatAnalysis_file$datapath, 
      type = "RDsMulti",
      isFileUploaded = !is.null(input$loadStatAnalysis_file)
    )
    
    if (!is.null(result$error)) {
      showAlert("Error", result$error, "error", 5000)
      manageSpinner(FALSE)
      return()
    }
    
    datapaths <- input$loadStatAnalysis_file$datapath
    for(dpath in 1:length(datapaths)){
      mess <- readRDS(datapaths[dpath])
      
      if(!(all(names(mess) %in% names(DataAnalysisModule)) ||
           all(names(mess) %in% names(elisaResult)) ||
           all(names(mess) %in% names(wbquantResult)) || 
           all(names(mess) %in% names(pcrResult)) ||
           all(names(mess) %in% names(cytotoxResult)) ||
           all(names(mess) %in% names(endocResult)))){
        showAlert("Error", paste(mess[["message"]],"\n The file must be RDs saved through the Data Analysis module."), "error", 5000)
        manageSpinner(FALSE)
        return()
      }
      
      DataStatisticModule$Flag <- TRUE
      
      if(all(names(mess) %in% names(wbquantResult)) || all(names(mess) %in% names(DataAnalysisModule))){
        DataStatisticModule$WB[[dpath]] <- mess$AdjRelDensitiy %>% mutate(DataSet = dpath)
      } else if(all(names(mess) %in% names(pcrResult)) || all(names(mess) %in% names(DataAnalysisModule))){
        DataAnalysisModule$PRCC[[dpath]]  <- mess
      } else if(all(names(mess) %in% names(endocResult)) || all(names(mess) %in% names(DataAnalysisModule))){
        DataAnalysisModule$ENDOC[[dpath]]  <- mess
      } else if(all(names(mess) %in% names(elisaResult)) || all(names(mess) %in% names(DataAnalysisModule))){
        DataAnalysisModule$ELISA[[dpath]]  <- mess
      } else if(all(names(mess) %in% names(cytotoxResult)) || all(names(mess) %in% names(DataAnalysisModule))){
        DataAnalysisModule$CYTOTOX[[dpath]]  <- mess
      }
    }
    manageSpinner(FALSE)
    showAlert("Success", "The RDs files have been uploaded with success", "success", 2000)
    return(NULL)
  })
  
  observeEvent(input$StatAnalysis, {
    if (input$StatAnalysis != "") {
      DataStatisticModule[[input$StatAnalysis]] -> results
      do.call(rbind, results) -> results
      
      res = resTTest = NULL
      resplot = ggplot()
      
      switch(input$StatAnalysis, 
          "WB" =  {
          res = results %>%
            select(DataSet, SampleName, AdjRelDens) %>%
            mutate(SampleName = gsub(pattern = "^[0-9]. ", x = SampleName, replacement = ""),
                   ColorSet = as.character(DataSet)) 
          
          points = res %>%
            mutate(SampleName = as.factor(SampleName))
          
          stats = points %>%
            group_by(SampleName) %>%
            summarise(Mean = mean(AdjRelDens), Sd = sd(AdjRelDens), .groups = 'drop')
          
          resplot = ggplot(stats, aes(x = SampleName, y = Mean)) + 
            geom_bar(stat="identity", color="black", fill = "#BAE1FF", position=position_dodge()) +
            geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd), width=.2, position=position_dodge(.9)) +
            geom_point(data = points, aes(x = SampleName, y = AdjRelDens, color = ColorSet), position = position_jitter(width = 0.2), size = 3) +
            theme_bw()
        }
      )
      
      output$TabStat = renderDT({stats})
      output$PlotStat = renderPlot({resplot})
      output$TabTTest = renderDT({resTTest})
    }
  })
  
  ### End Statistic ####
  
  
  #----------------------------------------------------------------------------------
  # OTHER ANALYSIS TO DO
  cytotoxResult  = reactiveValues(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    CYTOTOXcell_EXP = NULL,
    CYTOTOXcell_REP = NULL,
    CYTOTOXcell_SN = NULL,
    MapBaseline = NULL)
  
  # DOWNLOAD REPORT E RDS
  output$downloadReport <- downloadHandler(
    filename = function() {
      "Report.html"
    },
    content = function(file) {
      if (checkAnalysis()) {
        showAlert("Error", "no analyzes to download", "error", 5000)
        return(NULL)
      }
      
      manageSpinner(TRUE)
      parmsList = list(ResultList = reactiveValuesToList(DataAnalysisModule))
      rmarkdown::render("inst/shiny/report.Rmd",
                        output_file = file, output_format = "html_document",
                        params = parmsList)
      manageSpinner(FALSE)
      showAlert("Success", "Download completed successfully", "success", 2000)
    }
  )
  
  output$downloadRDSwholeAnalysis <- downloadHandler(
    filename = function() {
      paste('DataIntegrationModuleAnalysis-', Sys.Date(), '.RDs', sep='')
    },
    content = function(file) {
      if (checkAnalysis()) {
        showAlert("Error", "no analyzes to download", "error", 5000)
        return(NULL)
      }
      manageSpinner(TRUE)
      saveRDS(reactiveValuesToList(DataAnalysisModule), file = file)
      manageSpinner(FALSE)
      showAlert("Success", "Download completed successfully", "success", 2000)
    }
  )
  
  checkAnalysis <- function() {
    if (!is.null(wbResult$Normalizer) || !is.null(wbResult$Im) || !is.null(wbResult$Planes) ||
        !is.null(wbResult$TruncatedPanelsValue) || !is.null(wbResult$PanelsValue) ||
        !is.null(wbResult$Plots) || !is.null(wbResult$TruncatedPlots) || !is.null(wbResult$pl) ||
        !identical(wbResult$AUCdf, data.frame(SampleName = "-", Truncation = "-", AUC = "-")))
      return (FALSE)
    if (!is.null(wbquantResult$NormWBanalysis) || !is.null(wbquantResult$NormWBanalysis_filtered) ||
        !is.null(wbquantResult$WBanalysis) || !is.null(wbquantResult$WBanalysis_filtered) ||
        !is.null(wbquantResult$AdjRelDensitiy))
      return (FALSE)
    if (!is.null(pcrResult$Initdata) || !is.null(pcrResult$selectPCRcolumns) ||
        !is.null(pcrResult$data) || !is.null(pcrResult$PCRnorm) || !is.null(pcrResult$NewPCR) ||
        !is.null(pcrResult$plotPRC))
      return (FALSE)
    if (!is.null(elisaResult$Initdata) || !is.null(elisaResult$data) ||
        !is.null(elisaResult$TablePlot) || !is.null(elisaResult$dataFinal) ||
        !is.null(elisaResult$ELISAcell_EXP) || !is.null(elisaResult$ELISAcell_SN) ||
        !is.null(elisaResult$MapBaseline) || !is.null(elisaResult$MapBlank) ||
        !is.null(elisaResult$Tablestandcurve) || !is.null(elisaResult$Regression))
      return (FALSE)
    if (!is.null(cytotoxResult$Initdata) || !is.null(cytotoxResult$data) ||
        !is.null(cytotoxResult$TablePlot) || !is.null(cytotoxResult$dataFinal) ||
        !is.null(cytotoxResult$CYTOTOXcell_EXP) || !is.null(cytotoxResult$CYTOTOXcell_REP) ||
        !is.null(cytotoxResult$CYTOTOXcell_SN) || !is.null(cytotoxResult$MapBaseline))
      return (FALSE)
    if (!is.null(endocResult$Initdata) || !is.null(cytotoxResult$data) ||
        !is.null(endocResult$TablePlot) || !is.null(cytotoxResult$dataFinal) ||
        !is.null(endocResult$ENDOCcell_TIME) || !is.null(cytotoxResult$ENDOCcell_COLOR) ||
        !is.null(endocResult$MapBaseline) || !is.null(cytotoxResult$MapBlank))
      return (FALSE)
  
    return (TRUE)
  }
  
  # END DOWNLOAD
}

