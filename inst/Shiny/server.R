#shiny.maxRequestSize=1000*1024^2
#shiny.launch.browser = .rs.invokeShinyWindowExternal

Sys.setenv("DATAVERSE_SERVER" = "dataverse.harvard.edu")
APIkey_path = system.file("Data",".APIkey", package = "ORCA")

#source(system.file("Shiny","AuxFunctions.R", package = "ORCA"))
source("./inst/Shiny/AuxFunctions.R")

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  alert <- reactiveValues(alertContext = NULL)
  
  DataAnalysisModule <- reactiveValues(wbResult = NULL,
                                       wbquantResult = NULL,
                                       endocResult = NULL,
                                       elisaResult = NULL,
                                       pcrResult = NULL,
                                       cytotoxResult = NULL)
  
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
  
  
  # save everytime there is a change in the results
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
      loadTifFile()
    }
  })
  
  loadTifFile <- function() {
    file <- !is.null(input$imImport) && file.exists(input$imImport$datapath)
    mess = readfile(filename = input$imImport$datapath, type = "tif", file)
    
    
    if(setequal(names(mess), c("message", "call"))) {
      showAlert("Error", mess[["message"]], "error", 5000)
    } else {
      Flags$ShowTif <- TRUE
      wbResult$Im = mess
      updateTabsetPanel(session, "SideTabs", selected = "plane")
      alert$alertContext <- NULL
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
      print(Planes)
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
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "WB")
      
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
  
  # selecting rows
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
  
  # the relative density and adjusted is calculated
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
        else{ # we admit only one SampleName
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
    
    if( !is.null(pcrResult$Initdata) ) {
      # If data is already present, show a modal asking for confirmation to update
      showModal(modalDialog(
        title = "Important message",
        "Do you want to update the RT-PCR data already present?",
        easyClose = TRUE,
        footer= tagList(actionButton("confirmUploadPCR", "Update"),
                        modalButton("Cancel")
        )
      ))
    } else {
      mess = readfile(
        filename = input$PCRImport$datapath,
        type = "Excel", 
        isFileUploaded = !is.null(input$PCRImport) && file.exists(input$PCRImport$datapath),
        colname = TRUE, 
        namesAll = namesAll, 
        allDouble = FALSE, 
        colors = FALSE 
      )
      
      if (!is.null(mess$message) && mess$call == "") {
        showAlert("Error", mess$message, "error", 5000)
      } else {
        validate(
          need(!is.null(input$PCRImport) && file.exists(input$PCRImport$datapath),
               "Please select an RT-qPCR excel file!!")
        )
        pcrResult$Initdata = mess
        
        updateSelectInput(session, "PCR_gene",
                          choices = c("", colnames(pcrResult$Initdata)),
                          selected = "")
        updateSelectInput(session, "PCR_sample",
                          choices = c("", colnames(pcrResult$Initdata)),
                          selected = "")
        updateSelectInput(session, "PCR_value",
                          choices = c("", colnames(pcrResult$Initdata)),
                          selected = "")
        updateSelectInput(session, "PCR_time",
                          choices = c("", colnames(pcrResult$Initdata)),
                          selected = "")
        showAlert("Success", "The RT-qPCR excel has been uploaded with success", "success", 2000)
      }
    }
  })
  
  
  observeEvent(input$confirmUploadPCR,{
    
    pcrResult = pcrResult0
    
    FlagsPCR$norm=F 
    FlagsPCR$baseline = F
    
    output$PCRtables <- renderUI({ NULL })
    output$PCRtablesComp <- renderUI({ NULL })
    output$PCRplot <- renderPlot({ NULL })
    
    output$LoadingError_PCR <- renderText({
      validate(
        need(!is.null(input$PCRImport) && file.exists(input$PCRImport$datapath) ,
             "Please select an RT-PCR excel file!!" )
      )
      
      mess = readfile(
        filename = input$PCRImport$datapath,
        type = "Excel"
      )
      
      validate(
        need(!setequal(names(mess),c("message","call")) ,
             mess[["message"]])
      )
      
      pcrResult$Initdata = mess
      
      updateSelectInput(session,"PCR_gene",
                        choices = c("",colnames(pcrResult$Initdata)),
                        selected = ""
      )
      updateSelectInput(session,"PCR_sample",
                        choices = c("",colnames(pcrResult$Initdata)),
                        selected = ""
      )
      updateSelectInput(session,"PCR_value",
                        choices = c("",colnames(pcrResult$Initdata)),
                        selected = ""
      )
      updateSelectInput(session,"PCR_time",
                        choices = c("",colnames(pcrResult$Initdata)),
                        selected = ""
      )
      
      removeModal()
      
      "The RDs has been uploaded  with success"
    })
  })
  
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
  
  observeEvent(input$StatAnalysis, {
    if (input$StatAnalysis != "") {
      DataStatisticModule[[input$StatAnalysis]] -> results
      do.call(rbind, results) -> results
      
      res = resTTest = NULL
      resplot = ggplot()
      
      if (input$StatAnalysis == "WB") {
        res = results %>%
          select(DataSet, SampleName, AdjRelDens) %>%
          mutate(SampleName = gsub(pattern = "^[0-9]. ", x = SampleName, replacement = ""))
        
        # Usa pivot_longer invece di gather
        dati_punti = res %>% 
          tidyr::pivot_longer(cols = -SampleName, names_to = "DataSet", values_to = "AdjRelDens") %>% 
          mutate(SampleName = as.factor(SampleName))
        
        # Calcola media e sd per ogni SampleName
        stats = dati_punti %>%
          group_by(SampleName) %>%
          summarise(Mean = mean(AdjRelDens), Sd = sd(AdjRelDens))
        
        resplot = ggplot(stats, aes(x = SampleName, y = Mean)) + 
          geom_bar(stat="identity", color="black", fill = "#BAE1FF", position=position_dodge()) +
          geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd), width=.2, position=position_dodge(.9)) +
          geom_point(data = dati_punti, aes(x = SampleName, y = AdjRelDens), position = position_jitter(width = 0.2), color = "red") +
          theme_bw()
      }
      
      output$TabStat = renderDT({stats}) # Assicurati che stats sia il dataframe che vuoi mostrare
      output$PlotStat = renderPlot({resplot})
      output$TabTTest = renderDT({resTTest})
    }
  })
  
  
  ### End Statistic ####
  
  
  #----------------------------------------------------------------------------------
  # OTHER ANALYSIS TO DO
  elisaResult  = reactiveValues(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    ELISAcell_EXP = NULL,
    ELISAcell_SN = NULL,
    MapBaseline = NULL,
    MapBlank = NULL,
    Tablestandcurve = NULL,
    Regression = NULL)
  
  cytotoxResult  = reactiveValues(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    CYTOTOXcell_EXP = NULL,
    CYTOTOXcell_REP = NULL,
    CYTOTOXcell_SN = NULL,
    MapBaseline = NULL)
  
  endocResult = reactiveValues(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    ENDOCcell_TIME = NULL,
    ENDOCcell_SN = NULL,
    MapBaseline = NULL,
    MapBlank = NULL)
  
  
}

