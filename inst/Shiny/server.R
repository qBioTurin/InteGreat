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
      loadExcelFile()
    }
  })
  
  loadExcelFilePCR <- function() {
    alert$alertContext <- ""
    
    mess = readfile(
      filename = input$PCRImport$datapath,
      type = "Excel", 
      isFileUploaded = !is.null(input$PCRImport) && file.exists(input$PCRImport$datapath),
      colname = FALSE, 
      namesAll = namesAll, 
      allDouble = TRUE, 
      colors = TRUE 
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
    ENDOCcell_SN = NULL,
    MapBaseline = NULL,
    MapBlank = NULL)
  
  endocResult0 = list(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    ENDOCcell_TIME = NULL,
    ENDOCcell_SN = NULL,
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
  
  observeEvent(input$LoadENDOC_Button,{
    alert$alertContext <- "ENDOC-reset"
    if( !is.null(endocResult$Initdata) ) {
      shinyalert(
        title = "Important message",
        text = "Do you want to update the WB data already present, by resetting the previous analysis?",
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
  
      updateSelectizeInput(session, "ENDOCcell_SN", choices = character(0), selected = character(0))
      updateSelectizeInput(session, "ENDOCcell_TIME", choices = character(0), selected = character(0))
      
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
        endocResult$ENDOCcell_SN = mess$SNtable
        
        removeModal()
        loadEndocColor()
        showAlert("Success", "The RDs has been uploaded  with success", "success", 2000)
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
  
  loadEndocColor <- function() {
    observe({
      color_names <- names(FlagsENDOC$EXPcol)
      color_codes <- FlagsENDOC$EXPcol
      
      selected_color_index <- which(color_names == "selected" | color_codes == "#49ff00")
      if (length(selected_color_index) > 0) {
        color_names <- color_names[-selected_color_index]
        color_codes <- color_codes[-selected_color_index]
      }
      
      color_styles <- sapply(color_codes, function(color_code) {
        text_color <- if (grepl("^(#FFFFFF|#FFC000|#FFFF00|#D6D6D6|#EBEBEB)", color_code, ignore.case = TRUE)) "#000000" else "#FFFFFF"
        paste0("background-color: ", color_code, "; color: ", text_color, ";")
      })
      
      updatePickerInput(session, "colorDropdown",
                        choices = color_names,
                        options = list(`style` = "width: 100px;"),
                        choicesOpt = list(
                          elementId = color_names,
                          style = color_styles
                        ),
                        selected = NULL
      )
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
      
      get_formatted_data <- function(colors, color_names) {
        if (length(colors) == 0) {
          return(data.frame(Color = character(), Values = character(), Name = character()))
        }
        formatted_data <- vector("list", length(colors))
        for (i in seq_along(colors)) {
          matching_indices <- which(endocResult$ENDOCcell_SN == color_names[i], arr.ind = TRUE)
          if (nrow(matching_indices) > 0) {
            selected_values <- apply(matching_indices, 1, function(idx) {
              endocResult$Initdata[idx["row"], idx["col"]]
            })
            formatted_output <- paste(unlist(selected_values), collapse = " - ")
            formatted_data[[i]] <- data.frame(
              Color = sprintf("<div style='background-color: %s; padding: 10px; margin-right:20px; '></div>", colors[i]),
              Values = formatted_output,
              Name = color_names[i]
            )
          } else {
            formatted_data[[i]] <- data.frame(
              Color = sprintf("<div style='background-color: %s; padding: 10px;'></div>", colors[i]),
              Values = "No matching indices found.",
              Name = color_names[i]
            )
          }
        }
        return(do.call(rbind, formatted_data))
      }
      
      left_data <- get_formatted_data(left_colors, color_names[1:mid_point])
      right_data <- get_formatted_data(right_colors, color_names[(mid_point+1):length(color_codes)])
      
      # Aggiungi qui la stampa della struttura
      print("Left Data Structure:")
      print(str(left_data))
      print("Right Data Structure:")
      print(str(right_data))
      
      output$leftTable <- renderDataTable(
        left_data, escape = FALSE, 
        editable = list(target = "cell", 
                   disable = list(columns = 0:1)),
        options = list(
        dom = 't',
        paging = FALSE,
        info = FALSE,
        searching = FALSE, 
        columnDefs = list(
          list(width = '10px', targets = 0),
          list(width = '10px', targets = 1),
          list(width = '200px', targets = 2),
          list(className = 'dt-head-left dt-body-left', targets = 1)
        )
      ))
      output$rightTable <- renderDataTable(
        left_data, escape = FALSE, 
        editable = list(target = "cell", 
                        disable = list(columns = 0:1)),
        options = list(
        dom = 't',
        paging = FALSE,
        info = FALSE,
        searching = FALSE,
        editable= TRUE,
        columnDefs = list(
          list(width = '10px', targets = 0),
          list(width = '10px', targets = 1),
          list(width = '200px', targets = 2),
          list(className = 'dt-head-left dt-body-left', targets = 1)
        )
      ))
    })
  }
  
  observeEvent(input$colorDropdown, {
    req(input$colorDropdown)
    
    selectedColorName <- input$colorDropdown

    whiteKey <- names(FlagsENDOC$EXPcol)[FlagsENDOC$EXPcol == "white"]

    if (exists("originalColor") && originalColor != "" && any(endocResult$ENDOCcell_SN == "selected")) {
      originalIndices <- which(endocResult$ENDOCcell_SN == "selected", arr.ind = TRUE)

      if (length(originalIndices) > 0) {
        for (idx in 1:nrow(originalIndices)) {
          endocResult$ENDOCcell_SN[originalIndices[idx, "row"], originalIndices[idx, "col"]] <- originalColor
        }
      }
    }
    
    if (selectedColorName != "white" && selectedColorName != "#FFFFFF" && selectedColorName != whiteKey) {
      matchingIndices <- which(endocResult$ENDOCcell_SN == selectedColorName, arr.ind = TRUE)

      if (nrow(matchingIndices) > 0) {
        originalColor <<- selectedColorName  
        selectedValues <- apply(matchingIndices, 1, function(idx) {
          endocResult$Initdata[idx["row"], idx["col"]]
        })

        selectedValuesVector <- unlist(selectedValues)
        formattedOutput <- paste("Selected Values:", paste(selectedValuesVector, collapse = " - "))
        output$ENDOCSelectedValues <- renderText({ formattedOutput })
        
        apply(matchingIndices, 1, function(idx) {
          endocResult$ENDOCcell_SN[idx["row"], idx["col"]] <- "selected"
          updateSelectizeInput(inputId = "ENDOCcell_TIME",
                               selected = ifelse(is.null(endocResult$ENDOCcell_TIME[idx["row"], idx["col"]]),
                                                 "",
                                                 endocResult$ENDOCcell_TIME[idx["row"], idx["col"]])
          )
          updateSelectizeInput(inputId = "ENDOCcell_SN",
                               selected = ifelse(is.null(endocResult$ENDOCcell_SN[idx["row"], idx["col"]]),
                                                 "",
                                                 endocResult$ENDOCcell_SN[idx["row"], idx["col"]])
          )
        })
        endocResult$TablePlot <- NULL
        invalidateLater(1000, session)
      } else {
        output$ENDOCSelectedValues <- renderText("No matching indices found.")
      }
    } else {
      output$ENDOCSelectedValues <- renderText("No selected values")
    }
  })
  
  observeEvent(input$ENDOCcell_TIME, {
    if (!is.null(endocResult$ENDOCcell_TIME)) {
      selectedIndices <- which(endocResult$ENDOCcell_SN == "selected", arr.ind = TRUE)
      
      if (nrow(selectedIndices) > 0) {
        apply(selectedIndices, 1, function(idx) {
          endocResult$ENDOCcell_TIME[idx["row"], idx["col"]] <- input$ENDOCcell_TIME
        })
      }
    }
  })
  
  observeEvent(input$ENDOCcell_SN, {
    req(input$ENDOCcell_SN)  

    selectedIndices <- which(endocResult$ENDOCcell_SN == "selected", arr.ind = TRUE)
    if (nrow(selectedIndices) > 0) {
      updates <- list()  
      
      for (idx in seq_len(nrow(selectedIndices))) {
        cellCoo <- c(selectedIndices[idx, "row"], selectedIndices[idx, "col"])
        value.now <- input$ENDOCcell_SN
        value.bef <- endocResult$ENDOCcell_SN[cellCoo[1], cellCoo[2]]
        
        if (value.now != "" && value.now != value.bef) {
          updates[[length(updates) + 1]] <- list(cellCoo = cellCoo, value = value.now)
        }
      }
      
      if (length(updates) > 0) {
        for (update in updates) {
          cellCoo <- update$cellCoo
          value.now <- update$value
          endocResult$ENDOCcell_SN[cellCoo[1], cellCoo[2]] <- value.now
          endocResult$TablePlot$x$data[cellCoo[1], paste0("Col", cellCoo[2])] <- value.now
          
          if (nzchar(value.now) && !value.now %in% FlagsENDOC$AllExp) {
            exp <- unique(c(FlagsENDOC$AllExp, value.now))
            exp <- exp[exp != ""]
            FlagsENDOC$AllExp <- exp
          }
        }

      }
    } else showAlert("Error", "please, select before a row color", "error", 5000)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  ## update Baselines checkBox
  observeEvent(c(FlagsENDOC$AllExp,FlagsENDOC$BLANCHEselected),{
    if(length(FlagsENDOC$AllExp) >= 1){
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
    if(length(FlagsENDOC$AllExp) >= 1){
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
      return(list("Nothing",endocResult$ENDOCcell_TIME,endocResult$ENDOCcell_SN))
    }else{
      return(c(listReturn,list(endocResult$ENDOCcell_TIME,endocResult$ENDOCcell_SN)) )
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
        matExp =  as.matrix(endocResult$ENDOCcell_SN)
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
    matExp =  as.matrix(endocResult$ENDOCcell_SN)
    
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
                legend.background = element_rect(size=0.5,
                                                 linetype="solid",
                                                 colour ="black"))
      )
    }
  })
  
  ### End ENDOC analysis ####
  
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
        !is.null(endocResult$ENDOCcell_TIME) || !is.null(cytotoxResult$ENDOCcell_SN) ||
        !is.null(endocResult$MapBaseline) || !is.null(cytotoxResult$MapBlank))
      return (FALSE)
  
    return (TRUE)
  }
  
  # END DOWNLOAD
}

