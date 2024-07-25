# shiny.maxRequestSize=1000*1024^2
# shiny.launch.browser = .rs.invokeShinyWindowExternal

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
                                       facsResult = NULL,
                                       bcaResult = NULL,
                                       ifResult = NULL,
                                       facsResult = NULL)
  
  DataStatisticModule = reactiveValues(WB = list(),
                                       PCR = list(),
                                       ELISA = list(),
                                       ENDOC = list(),
                                       CYTOTOX = list(),
                                       IF = list(),
                                       FACS = list(),
                                       Flag = F)
  
  DataIntegrationModule <- reactiveValues(WB = list(),
                                          PCR = list(),
                                          ELISA = list(),
                                          ENDOC = list(),
                                          CYTOTOX = list(),
                                          IF = list(),
                                          FACS = list(),
                                          Stat = list(),
                                          Flag = F)
  
  MapAnalysisNames =c("WB", "WB comparison", "Endocytosis", "ELISA", "RT-qPCR", "Cytotoxicity", "FACS","BCA","IF") 
  names(MapAnalysisNames) =c("wbResult", "wbquantResult", "endocResult", "elisaResult", "pcrResult", "cytotoxResult", "facsResult","bcaResult","ifResult") 
  
  #### DATA INTEGRATION ####
  ### Omics ####
  
  Omics = reactiveValues(DefaultInit = NULL,
                         Default = NULL,
                         UserInit = NULL,
                         User = NULL,
                         Data = NULL, SelectedRows = list(Default = NULL, User = NULL), FinalSelectedRow = NULL)
  
  # default
  protfile = system.file("Data/Proteomic","pmic13104Simplified.xlsx", package = "ORCA")
  if(protfile == ""){
    # In this case it means that we are inside the docker
    protfile = "~/../home/data/Examples/Proteomic/pmic13104Simplified.xlsx"
  }
  Default = readxl::read_excel(protfile)
  Omics$DefaultInit = Omics$Default = Default
  Omics$Data = Default %>%
    dplyr::select(`Protein name`,iBAQ) %>%
    rename(Value = iBAQ, Name = `Protein name`)
  
  output$OmicsDefault_table <- renderDT(
    Omics$Default,
    filter = 'top',
    server = FALSE,
    options = list(lengthChange = FALSE,
                   autoWidth = FALSE,
                   columnDefs = list(list(
                     width = "125px",
                     targets = "_all"
                   )) ),
    rownames= FALSE
  )
  
  # user
  observeEvent(input$LoadOmics_Button,{
    output$LoadingError_Prot <- renderText({
      if(!is.null(input$LoadOmics_Button) && file.exists(input$LoadOmics_Button$datapath) )
      {
        "Please select an excel files!!"
      }
      
      mess = readfile(
        filename = input$LoadOmics_Button$datapath,
        type = "Excel"
      )
      
      validate(
        need(!setequal(names(mess),c("message","call")) ,
             paste(mess[["message"]],"\n all the files must be RDs saved from DataIntegrationModule app." ))
      )
      
      
      colmns = colnames(mess[,sapply(mess, is.numeric)])
      validate(
        need(length(colmns) > 0 ,
             "At least one column should be numeric."
        )
      )
      
      updateSelectInput(inputId = "RescalingColms_User",
                        choices = colmns ) 
      Omics$UserInit = Omics$User = mess
      "The excel file has been uploaded  with success"
    })
  })
  
  observeEvent(input$ResetOmics_Button,{
    Omics$User = NULL
    Omics$Default -> Default
    Omics$Data = Default %>%
      dplyr::select(`Protein name`,iBAQ) %>%
      rename(Value = iBAQ, Name = `Protein name`)
    
  })
  
  output$OmicsUser_table <- renderDT(
    Omics$User,
    filter = 'top',
    server = FALSE,
    selection = list(target = 'column'),
    options = list(lengthChange = FALSE,
                   autoWidth = FALSE,
                   columnDefs = list(list(
                     width = "125px", targets = "_all"
                   )) ),
    rownames= FALSE
  )
  
  observeEvent(input$OmicsUser_table_columns_selected,
               {
                 if(!is.null(Omics$User)){
                   output$LoadingError_Prot <- renderText({
                     
                     validate(
                       need(length(input$OmicsUser_table_columns_selected) == 2 ,
                            "Please select only two columns (one numeric and one character)." )
                     )
                     
                     prot = Omics$User[,input$OmicsUser_table_columns_selected+1]
                     validate(
                       need(all(sapply(prot, class) %in% c("character", "numeric")),
                            "Please select only two columns (one numeric and one character)." )
                     )
                     
                     protType = sapply(prot, class)
                     names(prot)[protType == "numeric"] = "Value"
                     names(prot)[protType == "character"] = "Name"
                     Omics$Data = prot
                   })
                 }
               })
  
  # Observe selected rows in OmicsDefault_table
  observeEvent(input$OmicsDefault_table_rows_selected, {
    if (!is.null(input$OmicsDefault_table_rows_selected)) {
      Omics$SelectedRows$Default <- Omics$Default[input$OmicsDefault_table_rows_selected, ]
    }
  })
  
  # Observe selected rows in OmicsUser_table
  observeEvent(input$OmicsUser_table_rows_selected, {
    if (!is.null(input$OmicsUser_table_rows_selected) && !is.null(Omics$User)) {
      Omics$SelectedRows$User <- Omics$User[input$OmicsUser_table_rows_selected, ]
    }
  })
  
  # Render selected rows tables
  output$SelectedOmicsDefault_table <- renderDT(
    Omics$SelectedRows$Default,
    filter = 'top',
    server = FALSE,
    options = list(lengthChange = FALSE,
                   autoWidth = FALSE,
                   columnDefs = list(list(
                     width = "125px",
                     targets = "_all"
                   ))),
    rownames = FALSE,
    selection = 'single'
  )
  
  output$SelectedOmicsUser_table <- renderDT(
    Omics$SelectedRows$User,
    filter = 'top',
    server = FALSE,
    options = list(lengthChange = FALSE,
                   autoWidth = FALSE,
                   columnDefs = list(list(
                     width = "125px",
                     targets = "_all"
                   ))),
    rownames = FALSE,
    selection = 'single'
  )
  
  # here, only one row can be selected
  observe({
    SelectedOmicsUser_table_rows_selected = input$SelectedOmicsUser_table_rows_selected
    User = Omics$SelectedRows$User
    proxy <- dataTableProxy('SelectedOmicsDefault_table')
    selectRows(proxy, NULL)
    Omics$FinalSelectedRow <-  User[SelectedOmicsUser_table_rows_selected, , drop = FALSE]
  })
  observe({
    SelectedOmicsDefault_table_rows_selected = input$SelectedOmicsDefault_table_rows_selected
    Default = Omics$SelectedRows$Default
    proxy <- dataTableProxy('SelectedOmicsUser_table')
    selectRows(proxy, NULL)
    Omics$FinalSelectedRow <-  Default[SelectedOmicsDefault_table_rows_selected, , drop = FALSE]  
  })
  
  # rescaling
  observeEvent(input$InputDefault_rescaling,{
    
    if(!( is.numeric(input$InputDefault_rescaling) && input$InputDefault_rescaling > 0)){
      showAlert("Error", "Please ensure the rescaling factor is a numeric value greater than 0.", "error", 5000)
      return()
    }
    
    if(!is.null(Omics$DefaultInit)){
      Omics$Default$iBAQ = Omics$DefaultInit$iBAQ/input$InputDefault_rescaling
    }
  })
  
  observeEvent(input$InputUser_rescaling,{
    if(!( length(input$RescalingColms_User)>0)){
      #showAlert("Error", "Please select at least one numeric column.", "error", 5000)
      #return()
      
      if(!( is.numeric(input$InputUser_rescaling) && input$InputUser_rescaling > 0)){
        showAlert("Error", "Please ensure the rescaling factor is a numeric value greater than 0.", "error", 5000)
        return()
      }
      
      colms = input$RescalingColms_User
      if(!is.null(Omics$UserInit)){
        Omics$User[,colms] = Omics$UserInit[,colms]/input$InputUser_rescaling
      }
    }
  })
  
  
  ### end Omics ###
  
  ### integration ####
  
  observeEvent(input$LoadIntG_Button,{
    manageSpinner(TRUE)
    
    result <- readfile(
      filename = input$IntGImport$datapath, 
      type = "RDsMulti",
      isFileUploaded = !is.null(input$LoadIntG_Button)
    )
    
    if (!is.null(result$error)) {
      showAlert("Error", result$error, "error", 5000)
      manageSpinner(FALSE)
      return()
    }
    
    datapaths <- input$IntGImport$datapath
    
    updateMultiValues(datapaths,DataIntegrationModule)
    
    manageSpinner(FALSE)
    showAlert("Success", "The RDs files have been uploaded with success", "success", 2000)
  })
  
  DataIntegrationresultListen <- reactive({
    reactiveValuesToList(DataIntegrationModule)
  })
  
  observeEvent(DataIntegrationresultListen(),{
    if(DataIntegrationModule$Flag){
      AnalysisNames = names(DataIntegrationModule)[names(DataIntegrationModule) != "Flag"]
      Analysis = rep(F,length(AnalysisNames) )
      names(Analysis) = AnalysisNames
      for(j in AnalysisNames)
        Analysis[j] = all(sapply(DataIntegrationModule[[j]], is.null))
      
      AnalysisNames = AnalysisNames[!Analysis]
      
      updateSelectizeInput(inputId = "IntegrAnalysis",
                           choices = c("",AnalysisNames),
                           selected = "")
      
      DataIntegrationModule$Flag = F
    }
  })
  
  IntegAnalysisResults <- reactive({
    Omics$FinalSelectedRow
    
    if (!is.null(input$IntegrAnalysis) && input$IntegrAnalysis != "") {
      results <- DataIntegrationModule[[input$IntegrAnalysis]]
      
      if(length(Omics$FinalSelectedRow$iBAQ) == 0){
        showAlert("Error", "An Omics value has to be selected", "error", 5000)
        return()
      }
      
      switch(input$IntegrAnalysis, 
             "Stat" = {
               results = do.call(rbind, results)
               finalTable = results %>% mutate(
                 OmicsValue = mean(Omics$FinalSelectedRow$iBAQ),
                 MeanScaled = OmicsValue*Mean/100 # I have to divide 100 because Mean is in %
               )%>% select(-sd) 
               
               list(Table = finalTable)
             },
             "WB" =  {
               
             },
             "IF" = {
               
             },
             "FACS" = {
               
             }
      )
    } else {
      NULL
    }
  })
  
  output$tables_IntG = renderTable({IntegAnalysisResults()$Table})
  #### end integration ###
  
  ### Other integration ####
  observeEvent(input$LoadOther_Button,{
    output$LoadingError_Other <- renderText({
      validate(
        need(!is.null(input$OtherImport) && all(file.exists(input$OtherImport$datapath)) ,
             "Please select one excel files!!" )
      )
      
      mess = readfile(
        filename = input$OtherImport$datapath,
        type = "Excel"
      )
      
      validate(
        need(!setequal(names(mess),c("message","call")) ,
             mess[["message"]]
        )
      )
      
      validate(
        need(all(c("character", "numeric") %in% sapply(mess, class)),
             "Please, the excel file must have two columns: one numeric and one character." )
      )
      
      DataIntegrationModule$otherTabs = mess
      "The excel file has been uploaded  with success"
    })
  })
  
  observe({
    if(!is.null(DataIntegrationModule$otherTabs)){
      
      updateSelectizeInput(inputId = "Selectprot_other",
                           choices = c("",unique(Omics$Data$Name)),
                           selected = "",
                           server = T )
    }
  })
  
  output$Other_table <- renderDT({
    DataIntegrationModule$otherTabs
  },
  filter = 'top',
  server = FALSE,
  selection = list(target = 'column'),
  options = list(lengthChange = FALSE,
                 autoWidth = TRUE),
  rownames= FALSE
  )
  
  output$Other_tableMean <- renderDT({
    DataIntegrationModule$otherTabsMean
  },
  filter = 'top',
  server = FALSE,
  options = list(lengthChange = FALSE,
                 autoWidth = TRUE),
  rownames= FALSE
  )
  
  observeEvent(input$Other_table_columns_selected,{
    if(!is.null(DataIntegrationModule$otherTabs) && !is.null(input$Other_table_columns_selected) && length(input$Other_table_columns_selected)>0){
      
      prot = Omics$Data
      data = DataIntegrationModule$otherTabs
      colmns = colnames(data)[input$Other_table_columns_selected + 1]
      colmns = names(sapply(data[,colmns], class)[ sapply(data[,colmns], class) == "character"])
      
      # check at leat one column selected is character
      if(length(colmns) > 0){
        colmnsNumeric = colnames(data)[-which(colnames(data) %in% colmns)]
        colmnsNumeric = colmnsNumeric[which(sapply(data[,colmnsNumeric], class) %in% "numeric" )]
        
        otherTabs = data %>% 
          tidyr::gather(colmnsNumeric, value = "Value", key= "Numeric Columns") %>%
          group_by_at(c("Numeric Columns",colmns)) %>%
          dplyr::summarise(Mean = mean(Value) ) %>%
          ungroup()
        
        DataIntegrationModule$otherTabsMean = otherTabs
      }
    }
  })
  
  observeEvent(input$Selectprot_other,{
    if(!is.null(DataIntegrationModule$otherTabsMean)){
      if(input$Selectprot_other != ""){
        prot = Omics$Data
        as.numeric(prot[prot$Name == input$Selectprot_other,"Value"]) -> pr
        
        DataIntegrationModule$otherTabsMean$`Rel.Prot.`  = DataIntegrationModule$otherTabsMean$Mean * pr
      }else{
        DataIntegrationModule$otherTabsMean$`Rel.Prot.`  = DataIntegrationModule$otherTabsMean$Mean
      }
    }
  })
  #### End other integration
  
  
  ### END DATA INTEGRATION ####
  
  #### DATA ANALYSIS ####
  
  #### BCA analysis ####
  observeEvent(input$NextBCAQuantif,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "tablesBCA")
  })
  #
  
  bcaResult = reactiveValues(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    dataQuant = NULL,
    BCAcell_EXP = "",
    BCAcell_SN = NULL,
    BCAcell_COLOR = NULL,
    MapBaseline = NULL,
    MapBlank = NULL,
    Tablestandcurve = NULL,
    Regression = NULL)
  
  bcaResult0 = list(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    dataQuant = NULL,
    BCAcell_EXP = "",
    BCAcell_SN = NULL,
    BCAcell_COLOR = NULL,
    MapBaseline = NULL,
    MapBlank = NULL,
    Tablestandcurve = NULL,
    Regression = NULL)
  
  left_data_bca <- reactiveVal()
  right_data_bca <- reactiveVal()
  
  # save everytime there is a change in the results
  BCAresultListen <- reactive({
    reactiveValuesToList(bcaResult)
  })
  observeEvent(BCAresultListen(), {
    DataAnalysisModule$bcaResult = reactiveValuesToList(bcaResult)
    DataAnalysisModule$bcaResult$Flags = reactiveValuesToList(FlagsBCA)
  })
  
  ##
  FlagsBCA <- reactiveValues(cellCoo = NULL,
                             AllExp = "",
                             BASEselected = "",
                             STDCselected = "",
                             BLANCHEselected = "",
                             EXPselected = "",
                             EXPcol = NULL)
  
  observeEvent(input$LoadBCA_Button,{
    alert$alertContext <- "BCA-reset"
    if( !is.null(bcaResult$Initdata) ) {
      shinyalert(
        title = "Important message",
        text = "Do you want to update the BCA data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else loadExcelFileBCA()
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "BCA-reset") {  
      resetPanel("BCA", flags = FlagsBCA, result = bcaResult)
      
      loadExcelFileBCA()
    }
  })
  
  loadExcelFileBCA <- function() {
    alert$alertContext <- ""
    mess = readfile(
      filename = input$BCAImport$datapath,
      isFileUploaded = !is.null(input$BCAImport) && file.exists(input$BCAImport$datapath),
      type = "Excel",
      allDouble = T,
      colname = F,
      colors = T
    )
    
    if(setequal(names(mess), c("message", "call"))) {
      showAlert("Error", mess[["message"]], "error", 5000)
    } else {
      bcaResult$Initdata = mess$x
      FlagsBCA$EXPcol = mess$fill
      bcaResult$BCAcell_COLOR = mess$SNtable
      bcaResult$BCAcell_SN = matrix("", nrow = nrow(bcaResult$BCAcell_COLOR), ncol = ncol(bcaResult$BCAcell_COLOR))
      
      removeModal()
      showAlert("Success", "The Excel has been uploaded  with success", "success", 2000)
    }
  }
  
  observe({
    if (!is.null(bcaResult$Initdata) && is.null(bcaResult$TablePlot)) {
      tableExcelColored(session = session,
                        Result = bcaResult, 
                        FlagsExp = FlagsBCA,
                        type = "Initialize")
    }
  })
  
  observeEvent(c(bcaResult$TablePlot,bcaResult$BCAcell_EXP), {
    if (!is.null(bcaResult$TablePlot)) {
      BCAtb <- bcaResult$TablePlot
      output$BCAmatrix <- renderDT(BCAtb, server = FALSE)
      
      if (!is.null(bcaResult$BCAcell_EXP) && !is.null(bcaResult$BCAcell_COLOR)) {
        matTime <- as.matrix(bcaResult$BCAcell_EXP)
        matExp <- as.matrix(bcaResult$BCAcell_COLOR)
        
        #if (!(all(matTime == "") || all(matExp == ""))) {
        mat <- as.matrix(bcaResult$Initdata)
        bcaV <- expand.grid(seq_len(nrow(mat)), seq_len(ncol(mat))) %>%
          rowwise() %>%
          mutate(values = mat[Var1, Var2])
        bcaT <- expand.grid(seq_len(nrow(matTime)), seq_len(ncol(matTime))) %>%
          rowwise() %>%
          mutate(time = matTime[Var1, Var2])
        bcaE <- expand.grid(seq_len(nrow(matExp)), seq_len(ncol(matExp))) %>%
          rowwise() %>%
          mutate(exp = matExp[Var1, Var2])
        bcaTot <- merge(bcaV, merge(bcaT, bcaE)) %>%
          na.omit()
        #filter(time != "", exp != "")
        
        bcaResult$data <- bcaTot
        #}
      }
    } 
  })
  
  observe({
    color_codes <- FlagsBCA$EXPcol
    color_names <- names(FlagsBCA$EXPcol)
    
    valid_colors <- color_codes != "white"
    color_codes <- color_codes[valid_colors]
    color_names <- color_names[valid_colors]
    
    mid_point <- ceiling(length(color_codes) / 2)
    left_colors <- color_codes[1:mid_point]
    right_colors <- color_codes[(mid_point+1):length(color_codes)]
    
    left_formatted_data <- get_formatted_data(left_colors, color_names[1:mid_point], bcaResult, bcaResult$BCAcell_EXP,"BCA")
    right_formatted_data <- get_formatted_data(right_colors, color_names[(mid_point+1):length(color_codes)], bcaResult, bcaResult$BCAcell_EXP, "BCA")
    
    left_data_bca(left_formatted_data)
    right_data_bca(right_formatted_data)
    
    output$leftTableBCA <- renderDataTable(
      left_data_bca(), 
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
    
    output$rightTableBCA <- renderDataTable(
      right_data_bca(), 
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
  
  observeEvent(input$leftTableBCA_cell_edit, {
    info <- input$leftTableBCA_cell_edit
    data <- left_data_bca() 
    updatedText <- updateTable("left", "BCA", info, data, bcaResult, FlagsBCA)
    
    output$BCASelectedValues <- renderText(updatedText)  
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  observeEvent(input$rightTableBCA_cell_edit, {
    info <- input$rightTableBCA_cell_edit
    data <- right_data_bca() 
    updatedText <- updateTable("right", "BCA", info, data, bcaResult, FlagsBCA)
    
    output$BCASelectedValues <- renderText(updatedText)  
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent(input$BCAmatrix_cell_clicked, {
    req(input$BCAmatrix_cell_clicked)  
    
    cellSelected = as.numeric(input$BCAmatrix_cell_clicked)
    FlagsBCA$cellCoo = cellCoo = c(cellSelected[1], cellSelected[2] + 1)
    
    allExp <- unique(na.omit(c(bcaResult$BCAcell_EXP)))  
    selectedExp <- ifelse(is.null(bcaResult$BCAcell_EXP[cellCoo[1], cellCoo[2]]), "", bcaResult$BCAcell_EXP[cellCoo[1], cellCoo[2]])
    
    updateSelectizeInput(inputId = "BCAcell_EXP", 
                         choices = c("",allExp),
                         selected = selectedExp)
    
    allSN <- unique(na.omit(c(bcaResult$BCAcell_SN)))  
    selectedSN <- ifelse(is.null(bcaResult$BCAcell_SN[cellCoo[1], cellCoo[2]]), "", bcaResult$BCAcell_SN[cellCoo[1], cellCoo[2]])
    
    updateSelectizeInput(inputId = "BCAcell_SN",
                         choices = c("",allSN),
                         selected = selectedSN)
  })
  
  observeEvent(input$BCAcell_SN, {
    if (!is.null(bcaResult$BCAcell_COLOR) && !is.null(FlagsBCA$cellCoo) && !anyNA(FlagsBCA$cellCoo)) {
      BCAtb = bcaResult$TablePlot
      cellCoo = FlagsBCA$cellCoo
      
      value.bef = bcaResult$BCAcell_COLOR[cellCoo[1], cellCoo[2]] 
      value.now = input$BCAcell_SN
      
      if (value.now != "" && value.now != value.bef) {
        currentValues <- bcaResult$Initdata[cellCoo[1], cellCoo[2]]
        
        bcaResult$BCAcell_COLOR[cellCoo[1], cellCoo[2]] = value.now
        bcaResult$BCAcell_SN[cellCoo[1], cellCoo[2]] = value.now
        BCAtb$x$data[cellCoo[1], paste0("Col", cellCoo[2])] = value.now
        
        if (!input$BCAcell_SN %in% FlagsBCA$AllExp) {
          exp = unique(c(FlagsBCA$AllExp, input$BCAcell_SN))
          FlagsBCA$AllExp = exp
        }
        
        tableExcelColored(session = session,
                          Result = bcaResult, 
                          FlagsExp = FlagsBCA,
                          type = "Update")
        
        output$BCASelectedValues <- renderText(paste("Updated value", paste(currentValues), ": sample name ", value.now))
        output$BCAmatrix <- renderDataTable({bcaResult$TablePlot})
      }
    } else return()
  }, ignoreInit = TRUE)
  observeEvent(input$BCAcell_EXP, {
    if (!is.null(bcaResult$BCAcell_EXP) && !is.null(FlagsBCA$cellCoo) && !anyNA(FlagsBCA$cellCoo)) {
      BCAtb = bcaResult$TablePlot
      cellCoo = FlagsBCA$cellCoo
      
      value.bef = bcaResult$BCAcell_EXP[cellCoo[1], cellCoo[2]] 
      value.now = input$BCAcell_EXP
      
      if (value.now != "" && value.now != value.bef) {
        currentValues <- bcaResult$Initdata[cellCoo[1], cellCoo[2]]
        
        bcaResult$BCAcell_EXP[cellCoo[1], cellCoo[2]] = value.now
        tableExcelColored(session = session,
                          Result = bcaResult, 
                          FlagsExp = FlagsBCA,
                          type = "Update"
        )
        
        output$BCASelectedValues <- renderText(paste("Updated value", paste(currentValues), ": Exp Condition ", value.now))
        output$BCAmatrix <- renderDataTable({bcaResult$TablePlot})
      }
    } else return()
  }, ignoreInit = TRUE)
  
  ## update checkBox
  observeEvent(FlagsBCA$AllExp,{
    if(length(FlagsBCA$AllExp) > 1){
      exp = FlagsBCA$AllExp
      exp = exp[exp != ""]
      updateSelectizeInput(session,"BCA_standcurve",
                           choices = exp,
                           selected = ifelse(FlagsBCA$STDCselected %in% exp,FlagsBCA$STDCselected,"") 
      )
    }
  })
  
  
  ## select the baselines, std curves, and blank
  observeEvent(input$BCA_standcurve,{
    FlagsBCA$STDCselected = input$BCA_standcurve
    FlagsBCA$EXPselected = FlagsBCA$AllExp[! FlagsBCA$AllExp %in% c(FlagsBCA$STDCselected,FlagsBCA$BASEselected,FlagsBCA$BLANCHEselected)]
  },ignoreNULL = F)
  # observeEvent(input$BCA_blanks,{
  #   FlagsBCA$BLANCHEselected = input$BCA_blanks
  #   FlagsBCA$EXPselected = FlagsBCA$AllExp[! FlagsBCA$AllExp %in% c(FlagsBCA$STDCselected,FlagsBCA$BASEselected,FlagsBCA$BLANCHEselected)]
  # },ignoreNULL = F)
  
  observe({
    BCA_blanks = input$BCA_blanks
    standcurve = bcaResult$Tablestandcurve
    
    isolate({
      if( !is.null(BCA_blanks) && !is.null(standcurve) && BCA_blanks == "yes" ){
        
        exp = FlagsBCA$EXPselected
        stcd = FlagsBCA$STDCselected
        exp = exp[exp != ""]
        expNotBlank = unique(c(stcd,exp))
        
        BlankV = standcurve %>% filter(Concentrations == min(Concentrations)) %>% summarize(M = mean(Measures)) %>% pull(M)
        
        standcurve$BlankValues = BlankV
        standcurve$MeasuresWithoutBlank =  standcurve$Measures - BlankV
        
        bcaResult$MapBlank = BlankV
        bcaResult$Tablestandcurve = standcurve
        
      }else{
        bcaResult$MapBlank = 0
      }
    })
  })
  
  ## update the data with the blank and baseline
  observe({
    if(!is.null(bcaResult$Initdata)){
      bcaResult$MapBlank -> MapBlank
      stcd = input$BCA_standcurve
      bcaResult$Regression -> Regression
      bcaTot = bcaResult$data
      isolate({
        if(!is.null(MapBlank) && !is.null(bcaTot)){
          
          bcaTotAverage = bcaTot %>%
            filter(exp != stcd) %>%
            #mutate(time = ifelse(exp %in% MapBlank$Blank, 0, time)) %>%
            group_by(time, exp) %>%
            summarize(AverageValues = mean(values), 
                      BlankValues = MapBlank,
                      AverageValuesWithoutBlank = AverageValues - BlankValues) %>%
            ungroup()
          
          if( !is.null(Regression) ){  
            bcamean = bcaTotAverage %>%
              dplyr::mutate(Quantification =  Regression$fun(AverageValuesWithoutBlank) ) %>%
              rename(SampleName = exp,ExpCondition = time) 
            
            output$BCAtables = renderDT(bcamean)
            
            bcaResult$dataQuant = bcamean
            
            bcameanNew = bcamean %>%
              select(SampleName,ExpCondition,Quantification) %>%
              rename(Ug = Quantification)
            bcaResult$dataFinal = bcameanNew
            output$BCAtablesUG = renderDT(bcameanNew)
            
          }else{
            output$BCAtables = renderDT(data.frame(Error = "No linear model!"))
          }
        }
      })
    }
  })
  
  observe({
    input$BCA_standcurve -> BCA_standcurve
    bcaResult$data -> data
    
    isolate({
      if(BCA_standcurve != ""){
        
        standcurve = data %>%
          filter(exp %in% BCA_standcurve) %>%
          # group_by(exp,time) %>%
          # summarise(AverageMeasures = mean(values)) %>%
          # ungroup() %>%
          select(exp,time,values) %>%
          rename(Measures = values, Concentrations = time) %>%
          filter(Concentrations != "")
        
        # If nothing changes w..r.t. the already saved table then I keep the old one!
        if(is.null(bcaResult$Tablestandcurve) || 
           !identical(bcaResult$Tablestandcurve,standcurve) )
        {
          bcaResult$Tablestandcurve = standcurve
        }
      } 
    })
  })
  
  ## update table standCurve
  observeEvent(bcaResult$Tablestandcurve,{
    output$BCA_Table_stdcurve <- DT::renderDataTable({
      DT::datatable( bcaResult$Tablestandcurve %>% 
                       select(exp, Concentrations,  Measures,  BlankValues, MeasuresWithoutBlank) %>%
                       rename(SampleName = exp),
                     selection = 'none',
                     # editable = list(target = "cell",
                     #                 disable = list(columns = 0:2) ),
                     rownames= FALSE,
                     options = list(scrollX = TRUE,
                                    searching = FALSE,
                                    dom = 't' # Only display the table
                     )
      )
    })
  })
  
  observeEvent(input$BCA_buttonRegression,{
    standcurve = bcaResult$Tablestandcurve
    
    if(!is.null(standcurve)){
      standcurve$Concentrations = as.numeric(standcurve$Concentrations)
      if(!all(is.na(standcurve$BlankValues)) ){
        standcurve = standcurve %>%
          select(Concentrations, MeasuresWithoutBlank) %>%
          rename( Measures = MeasuresWithoutBlank ) %>% na.omit()
        y="Measures Without Blank"
      }else{
        y="Measures"
        standcurve = standcurve %>% 
          select(Concentrations, Measures) %>% na.omit()
      }
      
      regressionPlot = ggplot(standcurve,aes(Concentrations, Measures)) +
        geom_point() +
        theme_bw()
      
      #if(input$regressionType == "Linear"){
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
                  vjust = "inward", hjust = "inward" )+
        labs(x = "Concentrations", y = y) 
      
      # }
      # else if(input$regressionType == "Quadratic")
      # {
      #   #this is not implemented
      #   standcurve$Concentrations2 = standcurve$Concentrations^2
      #   modelStancurve = lm(Measures~Concentrations+Concentrations2, data = standcurve)
      #   
      #   infoLM = data.frame(x = min(standcurve$Concentrations) + c(1,1),
      #                       y = max(standcurve$Measures) + c(2,1.75),
      #                       text = c( paste0("y = ", signif(modelStancurve$coef[[3]], 5), "x^2 + ",
      #                                        signif(modelStancurve$coef[[2]], 5), "x + ",signif(modelStancurve$coef[[1]],5 )),
      #                                 paste0("Adj R2 = ",signif(summary(modelStancurve)$adj.r.squared, 5))) )
      #   
      #   fun = paste0(modelStancurve$coef[[3]],"*x^2 + ",modelStancurve$coef[[2]],"*x + ",modelStancurve$coef[[1]] )
      #   
      #   regressionPlot =  regressionPlot  +
      #     geom_point() +
      #     stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,col="red")+
      #     geom_text(data= infoLM,
      #               aes(x = x, y = y, label =text ),
      #               vjust = "inward", hjust = "inward" )
      # }
      # else if(input$regressionType == "Hyperbola"){
      #   
      #   outNLreg = tryCatch(
      #     {
      #       modelStancurve<-nls(
      #         Measures ~ a*Concentrations/(b+Concentrations), 
      #         data = standcurve, #%>% group_by(Concentrations) %>% summarise(Measures = mean(Measures)),
      #         start = list(a = 1,b = 1)
      #       )
      #     }, 
      #     error = function(e){
      #       return(e)
      #     })
      #   
      #   if(!is.null(outNLreg$mess)){
      #     modelStancurve = NULL
      #     regressionPlot = ggplot()+ geom_text(data = data.frame(x = 1,y =1,text = paste0("Error: ",outNLreg$mess)),
      #                                          aes(x,y,label = text),color = "red")
      #   }else{
      #     modelStancurve = outNLreg
      #     coef = modelStancurve$m$getPars()
      #     r2 = 1- sum(modelStancurve$m$resid()^2)/(sum(( mean(standcurve$Measures) - modelStancurve$m$predict() )^2))
      #     
      #     infoLM = data.frame(x = min(standcurve$Concentrations) + c(1,1),
      #                         y = max(standcurve$Measures) + c(2,1.75),
      #                         text = c( paste0("y = ", signif(coef["a"], 5), "x / ( ",
      #                                          signif(coef["b"], 5), " + x ) "),
      #                                   paste0("R2 = ",signif(r2, 5))) )
      #     
      #     dfHyperbola = data.frame(x = seq(min(standcurve$Concentrations),max(standcurve$Concentrations),length.out = 20)) %>%
      #       mutate(y = (coef["a"]*x/((coef["b"]+x)) ) )
      #     
      #     fun = paste0(coef["b"],"*x/(",coef["a"],"-x)")
      #     
      #     regressionPlot =  regressionPlot  +
      #       geom_point() +
      #       geom_line(data = dfHyperbola,aes(x = x,y = y),size = 1,col="red" )+
      #       geom_text(data= infoLM,
      #                 aes(x = x, y = y, label =text ),
      #                 vjust = "inward", hjust = "inward" )
      #   }
      # }
      
      bcaResult$Regression = list(data = modelStancurve, plot = regressionPlot, fun = function(x){ eval( parse(text = fun ) ) } )
      
    }else{
      regressionPlot = ggplot()
    }
    output$BCAregression <- renderPlot(regressionPlot)
  })
  
  observeEvent(input$confirmBCA_UGvalue,{
    
    bcaResult$dataFinal -> bcamean
    
    if(!is.null(bcamean)){
      BCA_UGvalue = as.numeric(input$BCA_UGvalue)
      if(is.na(BCA_UGvalue)){
        showAlert("Error", "The value must be numeric!", "error", 5000)
        return()
      }else{
        bcamean[,paste0("Ug/",BCA_UGvalue)] =  bcamean$Ug/BCA_UGvalue
        
        bcaResult$dataFinal = bcamean
        
        output$BCAtablesUG = renderDT(bcamean)
        
      }
    }else{
      showAlert("Error", "The quantification step is missing!", "error", 5000)
      return()
    }
  })
  
  output$downloadBCAAnalysis <- downloadHandler(
    filename = function() {
      paste('BCAanalysis-', Sys.Date(), '.zip', sep='')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      nomeRDS <- paste0("BCA_analysis-", Sys.Date(), ".rds")
      nomeXLSX <- paste0("BCA_analysis-", Sys.Date(), ".xlsx")
      
      tempRdsPath <- file.path(tempDir, nomeRDS)
      tempXlsxPath <- file.path(tempDir, nomeXLSX)
      
      results <- DataAnalysisModule$bcaResult
      saveRDS(results, file = tempRdsPath)
      
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "BCA")
      
      utils::zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
      
    } 
  )
  
  
  # save everytime there is a change in the results
  BCAresultListen <- reactive({
    reactiveValuesToList(bcaResult)
  })
  observeEvent(BCAresultListen(), {
    DataAnalysisModule$bcaResult = reactiveValuesToList(bcaResult)
    DataAnalysisModule$bcaResult$Flags = reactiveValuesToList(FlagsBCA)
  })
  
  ### End BCA analysis ####  
  
  #### IF analysis ####
  
  ifResult = reactiveValues(
    Initdata = NULL,
    selectIFcolumns = NULL,
    data = NULL,
    StatData = NULL,
    FinalData = NULL,
    SubStatData = NULL,
    resplot = NULL,
    TTestData = NULL
  )
  ifResult0 = reactiveValues(
    Initdata = NULL,
    selectIFcolumns = NULL,
    data = NULL,
    StatData = NULL,
    FinalData = NULL,
    SubStatData = NULL,
    resplot = NULL,
    TTestData = NULL
  )
  
  # save everytime there is a change in the results
  IFresultListen <- reactive({
    reactiveValuesToList(ifResult)
  })
  observeEvent(IFresultListen(), {
    DataAnalysisModule$ifResult = reactiveValuesToList(ifResult)
    DataAnalysisModule$ifResult$Flags = reactiveValuesToList(FlagsIF)
  })
  
  ## next buttons
  observeEvent(input$NextQuantif,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "tablesIF")
  })
  observeEvent(input$NextifPlots,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "plotsIF")
  })
  
  FlagsIF <- reactiveValues(norm=F, 
                            baseline = F)
  
  
  observeEvent(input$LoadIF_Button,{
    alert$alertContext <- "IF-reset"
    if( !is.null(ifResult$Initdata) ) {
      shinyalert(
        title = "Important message",
        text = "Do you want to update the WB data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else loadExcelFileIF()
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "IF-reset") {  
      resetPanel("IF", flags = FlagsIF, result = ifResult)
      loadExcelFileIF()
    }
  })
  
  loadExcelFileIF <- function() {
    alert$alertContext <- ""
    
    mess = readfile(
      filename = input$IFImport$datapath,
      type = "Excel", 
      isFileUploaded = !is.null(input$IFImport) && file.exists(input$IFImport$datapath),
    )
    
    
    if (!is.null(mess$message) && mess$call == "") {
      showAlert("Error", mess$message, "error", 5000)
      return()
    }
    
    ifResult$Initdata = mess
    
    updateSelectInput(session,
                      "IF_expcond",
                      choices = c("", colnames(ifResult$Initdata)),
                      selected = "")
    
    showAlert("Success", "The IF excel has been uploaded with success", "success", 2000)
    # change pannel
    updateTabsetPanel(session = session, "SideTabs",
                      selected = "tablesIF")
  }
  
  observeEvent(input$IF_expcond,{
    tryCatch({
      if( !is.null(ifResult$Initdata) && input$IF_expcond != ""){
        selectIFcolumns = input$IF_expcond
        selectIFcolumns = selectIFcolumns[selectIFcolumns!= ""]
        
        IFdata = ifResult$Initdata
        colNames = colnames(IFdata)
        colNames[ colNames == selectIFcolumns] = "ExpCond"
        colnames(IFdata) = colNames
        
        IFdataCalc = IFdata %>%
          group_by(ExpCond) %>%
          mutate(nRow = 1:n()) %>%
          ungroup() %>%
          tidyr::gather(-ExpCond,-nRow, value = "Values", key = "Vars") %>%
          group_by(ExpCond,nRow) %>%
          mutate(Tot = sum(Values), Perc =  Values/Tot*100)
        
        IFinalData = IFdataCalc  %>%
          tidyr::pivot_wider( names_from = Vars, names_glue = "{Vars}_{.value}", 
                              values_from = c(Values, Perc)) %>% ungroup() %>% select(-nRow)
        
        statisticData = IFdataCalc %>% 
          group_by(ExpCond,Vars) %>%
          summarise(MeanValues = mean(Values),
                    sdValues = sd(Values),
                    MeanPerc = mean(Perc),
                    sdPerc = sd(Perc),
                    MeanTot = mean(Tot),
                    sdTot = sd(Tot)) %>%
          tidyr::pivot_wider( names_from = Vars, names_glue = "{Vars}_{.value}", 
                              values_from = c(MeanValues, sdValues,MeanPerc,sdPerc)) %>% ungroup()
        
        ifResult$StatData = statisticData
        ifResult$FinalData = IFinalData
        
        updateSelectInput("IF_TTestvariable",session = session, choices = unique(IFdataCalc$Vars))
        
        output$IFtable = renderDT({
          DT::datatable( IFinalData,
                         selection = 'none',
                         # editable = list(target = "cell",
                         #                 disable = list(columns = 0:2) ),
                         rownames= FALSE,
                         options = list(scrollX = TRUE,
                                        searching = FALSE,
                                        dom = 't' # Only display the table
                         )
          )
        })
        
        output$IFtable_stat = renderDT({
          DT::datatable(statisticData,
                        selection = 'none',
                        # editable = list(target = "cell",
                        #                 disable = list(columns = 0:2) ),
                        rownames= FALSE,
                        options = list(scrollX = TRUE,
                                       searching = FALSE,
                                       dom = 't' # Only display the table
                        )
          )
        })
      }
      else{
        
        updateSelectInput("IF_TTestvariable",session = session, choices = "",selected = "")
        
        output$IFtable = renderDT({ NULL })
        
        output$IFtable_stat = renderDT({ NULL })
      }
    }, error = function(e) {
      showAlert("Error",paste("An error occurred:", e$message), "error", 2000)
    })
  })
  
  observeEvent(input$IF_TTestvariable,{
    ifResult$FinalData -> IFinalData
    input$IF_TTestvariable -> varSel
    
    if(varSel != "" && !is.null(IFinalData)){
      IFinalData[,c("ExpCond", paste0(varSel,"_Perc"))] -> SubData
      colnames(SubData) = c("ExpCond", "Values")
      
      SubDataStat = SubData %>% group_by(ExpCond) %>% summarise(Mean = mean(Values), sd = sd(Values))
      
      resTTest = testStat.function(SubData, varSel)
      
      resplot <- ggplot(SubDataStat, aes(x = ExpCond, y = Mean)) + 
        geom_bar(stat="identity", color="black", fill = "#BAE1FF", position=position_dodge()) +
        geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(.9)) +
        geom_point(data = SubData, aes(x = ExpCond, y = Values, color = ExpCond),
                   position = position_jitter(width = 0.2), size = 3) +
        theme_bw()+
        labs(color = "Experimental Condition", y = "Percentages (%)" , x = "")
      
      output$IFsummarise_plot = renderPlot({resplot})
      
      output$IFsummariseMean = renderDT({
        DT::datatable(SubDataStat,
                      selection = 'none',
                      rownames= FALSE,
                      options = list(scrollX = TRUE,
                                     searching = FALSE,
                                     dom = 't' # Only display the table
                      )
        )
      })
      
      output$IFtable_ttest = renderDT({
        DT::datatable(resTTest,
                      selection = 'none',
                      rownames= FALSE,
                      options = list(scrollX = TRUE,
                                     searching = FALSE,
                                     dom = 't' # Only display the table
                      )
        )
      })
      
      ifResult$SubStatData = SubData
      ifResult$TTestData = resTTest
      ifResult$resplot = resplot
    }else{
      output$IFsummarise_plot = renderPlot({NULL})
      output$IFsummariseMean = renderDT({ NULL })
    }
    
  })
  
  
  observe({
    DataAnalysisModule$ifResult = reactiveValuesToList(ifResult)
  })
  
  output$downloadIFAnalysis <- downloadHandler(
    filename = function() {
      paste('IFanalysis-', Sys.Date(), '.zip', sep='')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      
      nomeRDS <- paste0("IFanalysis-", Sys.Date(), ".rds")
      nomeXLSX <- paste0("IFanalysis-", Sys.Date(), ".xlsx")
      
      tempRdsPath <- file.path(tempDir, nomeRDS)
      tempXlsxPath <- file.path(tempDir, nomeXLSX)
      
      results <- DataAnalysisModule$ifResult
      saveRDS(results, file = tempRdsPath)
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "IF")
      
      utils::zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
    },
  )
  
  #### END IF analysis ####
  
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
      
      utils::zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
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
      
      maxPanelsValue = max(wbResult$PanelsValue$Values) 
      wbResult$AUCdf -> AUCdf
      #AUCdf.new <- AUCdf[length(AUCdf$Truncation),]
      
      lastTrunc = AUCdf %>% 
        group_by(SampleName) %>%
        filter(SampleName == IDlane, row_number()==n() ) %>%
        ungroup() 
      
      # if(length(lastTrunc$Truncation) > 0 & lastTrunc$Truncation != "-")
      #   AUCdf.new$Truncation <- lastTrunc$Truncation
      # else
      #   
      #   
      # AUCdf.new$SampleName <- IDlane
      
      AUCdf.new = lastTrunc
      
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
                          min = round(min(PanelsValue$Values[PanelsValue$ID == IDlane]),digits = 3),
                          max = round(max(PanelsValue$Values[PanelsValue$ID == IDlane]),digits = 3),
                          value = round(min(PanelsValue$Values[PanelsValue$ID == IDlane],digits = 3) ) )
        
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
                                 AdjRelDensity = NULL
  )
  wbquantResult0 = list(NormWBanalysis = NULL,
                        NormWBanalysis_filtered = NULL,
                        WBanalysis = NULL,
                        WBanalysis_filtered = NULL,
                        RelDensitiy = NULL,
                        AdjRelDensity = NULL
  )
  FlagsWBquant = reactiveValues(BothUploaded = F)
  
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
                                 AdjRelDensity = NULL
  )
  
  loadWBanalysi4Quant = function(datapath,wbResult){
    manageSpinner(TRUE)
    names(wbResult) -> namesAll
    
    mess = readfile(
      filename = datapath,
      type = "RDs",
      namesAll = c(namesAll, "Flags")
    )
    
    if( setequal(names(mess),c("message","call")) ) {
      showAlert("Error", mess[["message"]], "error", 5000)
      return(mess[["message"]]) 
    }
    
    #wbquantResult$NormWBanalysis = mess
    manageSpinner(FALSE)
    
    showAlert("Success", "The RDS has been uploaded with success", "success", 2000)
    return(mess)
  }
  
  observeEvent(input$actionB_loadingNormWB, {
    alert$alertContext <- "WBNormQuant-reset"
    if(!is.null(wbquantResult$NormWBanalysis) ) { 
      shinyalert(
        title = "Important message",
        text = "Do you want to update the WB data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else{
      wbquantResult$NormWBanalysis = loadWBanalysi4Quant(input$NormWBImport$datapath,wbResult)
    }
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "WBNormQuant-reset") {  
      
      wbquantResult$NormWBanalysis = NULL
      FlagsWBquant$BothUploaded = F
      updateSelectInput("IdLaneNorm_RelDens",
                        session = session,
                        choices = "Nothing selected",
                        selected = "Nothing selected")
      
      output$AUC_RelDens <- renderDT(NULL)
      output$AUC_AdjRelDens <- renderDT(NULL)
      
      wbquantResult$NormWBanalysis = loadWBanalysi4Quant(input$NormWBImport$datapath,wbResult)
    }
  })
  
  observeEvent(input$actionB_loadingWB, {
    alert$alertContext <- "WBQuant-reset"
    if(!is.null(wbquantResult$NormWBanalysis) ) { 
      shinyalert(
        title = "Important message",
        text = "Do you want to update the WB data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else{
      wbquantResult$WBanalysis = loadWBanalysi4Quant(input$WBImport$datapath,wbResult)
      wbquantResult$WBanalysis_filtered = NULL
    }
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "WBQuant-reset") {  
      
      wbquantResult$WBanalysis = NULL
      FlagsWBquant$BothUploaded = F
      
      output$AUC_RelDens <- renderDT(NULL)
      output$AUC_AdjRelDens <- renderDT(NULL)
      
      wbquantResult$WBanalysis = loadWBanalysi4Quant(input$WBImport$datapath,wbResult)
      wbquantResult$WBanalysis_filtered = NULL    }
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
        
        choices = AUCdf[indexesWB,]$SampleName
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
  
  observe({
    FlagsWBquant$BothUploaded
    input$AUC_WB_rows_selected
    input$AUC_WBnorm_rows_selected
    
    input$IdLaneNorm_RelDens -> IdLaneNorm_RelDens
    isolate({
      table =  data.frame(SampleName = "-",
                          Truncation = "-", 
                          Truncation_Norm = "-",
                          AUC = "-", 
                          RelDens = "-",
                          AUC_Norm = "-",
                          RelDens_Norm = "-")
      tableAdjRel = data.frame(SampleName = "-",
                               AdjRelDensity = "-")
      barPlotAdjRelDens = ggplot()
      if(!is.null(wbquantResult$WBanalysis_filtered) & !is.null(wbquantResult$NormWBanalysis_filtered) &&
         IdLaneNorm_RelDens != "" && IdLaneNorm_RelDens != "Nothing selected"){
        IdLaneNorm_RelDens = input$IdLaneNorm_RelDens
        
        tbWBnormDen = wbquantResult$NormWBanalysis_filtered %>%
          filter(SampleName ==IdLaneNorm_RelDens[1]) %>%
          pull(AUC)
        
        tbWBnorm = wbquantResult$NormWBanalysis_filtered %>% 
          rename(AUC_Norm = AUC,
                 Truncation_Norm = Truncation) %>%
          mutate(RelDens_Norm = AUC_Norm/tbWBnormDen)
        
        tbWBDen = wbquantResult$WBanalysis_filtered %>%
          filter(SampleName ==IdLaneNorm_RelDens[1]) %>% 
          pull(AUC)
        
        tbWB = wbquantResult$WBanalysis_filtered %>%
          mutate(RelDens = AUC/tbWBDen)
        
        
        if(!is.null(tbWBnorm) & !is.null(tbWB) ){
          if(!all(table(tbWB$SampleName)==1) && !all(table(tbWBnorm$SampleName)==1) ){
            showAlert("Error",  "Only one sample name per panel is allowed", "error", 5000)
            return()
          }
          else{ 
            table = merge(tbWB,tbWBnorm,all =T) %>%
              dplyr::select(SampleName, Truncation, AUC, RelDens, AUC_Norm, RelDens_Norm) 
            tableAdjRel = table %>% na.omit() %>% mutate(AdjRelDensity = RelDens/RelDens_Norm) %>% select(SampleName,AdjRelDensity)
            barPlotAdjRelDens = tableAdjRel  %>%
              ggplot() +
              geom_bar(aes(x = SampleName,
                           y = AdjRelDensity,
                           fill = SampleName ),
                       stat = "identity" ) +
              theme_bw()
          }
        }
      }
      
      wbquantResult$RelDensitiy = table
      wbquantResult$AdjRelDensity = tableAdjRel
      
      output$AUC_RelDens <- renderDT(
        table,
        filter = 'top',
        server = FALSE,
        options = list(lengthChange = FALSE, autoWidth = TRUE),
        rownames= FALSE
      )
      
      output$AUC_AdjRelDens <- renderDT(
        tableAdjRel ,
        server = FALSE,
        options = list(lengthChange = FALSE, autoWidth = TRUE),
        rownames= FALSE
      )
      output$plot_AdjRelDens <- renderPlot({
        barPlotAdjRelDens
      })
      
    })
  })
  
  output$downloadWBquantAnalysis <- downloadHandler(
    filename = function() {
      paste0('WBquantanalysis-', Sys.Date(), '.zip')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      tempRDSPath <- file.path(tempDir, paste0("QuantificAnalysis-", Sys.Date(), ".rds"))
      tempExcelPath <- file.path(tempDir, paste0("QuantificAnalysis-", Sys.Date(), ".xlsx"))
      
      resultsRDS <- DataAnalysisModule$wbquantResult
      saveRDS(resultsRDS, tempRDSPath)
      
      resultsExcel <- DataAnalysisModule$wbquantResult
      saveExcel(filename = tempExcelPath, ResultList=resultsExcel, analysis = "WB comparison")
      
      utils::zip(file, files = c(tempRDSPath, tempExcelPath), flags = "-j")
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
    plotPCR = NULL,
    NewPCR = NULL,
    AllGenesFoldChangePlot = NULL,
    AllGenesFoldChangeTable = NULL)
  
  pcrResult0 = list(
    Initdata = NULL,
    selectPCRcolumns = NULL,
    data = NULL,
    PCRnorm = NULL,
    BaselineExp = NULL,
    plotPCR = NULL,
    NewPCR = NULL,
    AllGenesFoldChangePlot = NULL,
    AllGenesFoldChangeTable = NULL)
  
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
      type = "ExcelMulti", 
      isFileUploaded = !is.null(input$PCRImport) && all(file.exists(input$PCRImport$datapath)),
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
        group_by(Sample,Gene,Time) %>% mutate(Value = as.numeric(Value)) %>%
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
        dplyr::select(-Sample,  -HousekGeneMean, -HousekGeneSd)
      
      PCRstep4 = merge(BaselinePCR,PCRstep3,all.y = T,by=c("Gene","Time","HousekGene") )
      
      PCRstep5 = PCRstep4 %>%
        group_by(Sample,Gene,Time,HousekGene) %>%
        dplyr::summarize(
          ddCt = dCt - BaselinedCt,
          Q = 2^{-ddCt},
          Sd = Sd,
          Mean = Mean)%>%
        ungroup()
      
      AllGenes = unique(PCR$Gene)
      pcrResult$NewPCR = PCRstep5
      
      PCRstep5 = PCRstep5 %>% mutate(GeneH = paste(Gene, ", Housekeeping: ",HousekGene))
      
      # plot1 = lapply(unique(PCRstep5$GeneH),function(g){
      #     ggplot(data = PCRstep5 %>% filter(GeneH == g),
      #                  aes(x= as.factor(Time), y = ddCt, col = Sample)) + 
      #     facet_wrap(~GeneH, ncol = 1) +
      #     geom_jitter(width = 0.1, height = 0,size = 2)+
      #     theme_bw()+
      #     labs(x = "Time", y = "DDCT")
      #   })
      # 
      # plot2 = lapply(unique(PCRstep5$GeneH),function(g){
      #   ggplot(data = PCRstep5 %>% filter(GeneH == g),
      #          aes(x= as.factor(Time), y = Q, col = Sample)) + 
      #     facet_wrap(~GeneH, ncol = 1) +
      #     geom_jitter(width = 0.1, height = 0,size = 2)+
      #     theme_bw()+
      #     labs(x = "Time", y = "2^(-DDCT)")
      # })
      # 
      # plot1All = ggplot(data = PCRstep5,
      #          aes(x= as.factor(Time), y = ddCt, col = Sample)) + 
      #     facet_wrap(~GeneH, ncol = 1) +
      #     geom_jitter(width = 0.1, height = 0,size = 2)+
      #     theme_bw()+
      #     labs(x = "Time", y = "DDCT")
      # 
      # plot2All = ggplot(data = PCRstep5,
      #          aes(x= as.factor(Time), y = Q, col = Sample)) + 
      #     facet_wrap(~GeneH, ncol = 1) +
      #     geom_jitter(width = 0.1, height = 0,size = 2)+
      #     theme_bw()+
      #     labs(x = "Time", y = "2^(-DDCT)")
    }
  })
  # pcr plot
  
  observeEvent(pcrResult$NewPCR,{
    NewPCR = req(pcrResult$NewPCR)
    
    NewPCRFiltered = NewPCR %>% filter(!is.na(ddCt) | Sample == pcrResult$BaselineExp )
    updateSelectizeInput(inputId = "HousKgene_plot", choices = c("",unique(NewPCR$HousekGene)),selected = "" )
    updateSelectizeInput(inputId = "Gene_plot",choices = c("",unique(NewPCR$Gene)),selected = "" )
    
    updateSliderInput(session = session, inputId = "CutFoldChange_slider",
                      min = min(-NewPCRFiltered$ddCt,na.rm = T), max = max(-NewPCRFiltered$ddCt,na.rm = T), 
                      value = min(-NewPCRFiltered$ddCt,na.rm = T))
  })
  
  observe({
    NewPCR = req(pcrResult$NewPCR)
    cut = req(input$CutFoldChange_slider)
    
    PCRstep5 = NewPCR %>%
      mutate(GeneH = paste(Gene, ", Housekeeping: ",HousekGene) ) %>% 
      filter(!is.na(ddCt), Sample != pcrResult$BaselineExp ) %>% 
      mutate(Cut = if_else(-ddCt >= cut, "Greater", "Smaller"))
    
    pl = ggplot(PCRstep5,aes(x = Gene, y = -ddCt)) +
      geom_point(aes(color = Cut), size = 3) +
      geom_hline(yintercept = cut, color = "red", linetype = "dashed")+
      facet_wrap(~HousekGene, ncol = 1)+
      theme_bw()+
      labs(x="",y= "Log2(Q)", title = "")+
      scale_color_manual(values = c("Greater" = "#0072B2", "Smaller" = "#56B4E9")) +
      theme(axis.text.x=element_blank(), 
            axis.ticks.x=element_blank()) 
    
    output$FoldchangeAllGenesPlot = renderPlot({  pl  })
    
    table  =  PCRstep5 %>% rename(DDCT = ddCt, `2^(-DDCT)` = Q) %>%
      filter(-DDCT >= cut,!is.na(DDCT), Sample != pcrResult$BaselineExp )
    
    output$AllGenesTable = renderTable({ table })
    
    pcrResult$AllGenesFoldChangePlot = pl
    pcrResult$AllGenesFoldChangeTable = table
    
  })
  
  observe({
    input$Gene_plot -> gene
    input$HousKgene_plot -> Hgene
    isolate({
      if(Hgene != "" && gene != ""){
        
        PCRstep5 = pcrResult$NewPCR %>% 
          filter(HousekGene == Hgene, Gene == gene) %>%
          mutate(GeneH = paste(Gene, ", Housekeeping: ",HousekGene))
        
        plot1 = 
          ggplot(data = PCRstep5,
                 aes(x= as.factor(Time), y = ddCt, col = Sample)) + 
          facet_wrap(~GeneH, ncol = 1) +
          geom_jitter(width = 0.1, height = 0,size = 2)+
          theme_bw()+
          labs(x = "Time", y = "DDCT")
        
        plot2 = ggplot(data = PCRstep5 ,
                       aes(x= as.factor(Time), y = Q, col = Sample)) + 
          facet_wrap(~GeneH, ncol = 1) +
          geom_jitter(width = 0.1, height = 0,size = 2)+
          theme_bw()+
          labs(x = "Time", y = "2^(-DDCT)")
        
        output$SingleGenePlot = renderPlot({plot1|plot2})
        output$SingleGeneTable = renderTable({
          PCRstep5 %>% rename(DDCT = ddCt, `2^(-DDCT)` = Q)
        })
      }
    })
  })
  
  observeEvent(input$SavePCRplot,{
    input$Gene_plot -> gene
    input$HousKgene_plot -> Hgene
    isolate({
      if(Hgene != "" && gene != ""){
        
        plotList = pcrResult$plotPCR 
        
        if(is.null(plotList)) plotList = list()
        
        PCRstep5 = pcrResult$NewPCR %>% 
          filter(HousekGene == Hgene, Gene == gene) %>%
          mutate(GeneH = paste(Gene, ", Housekeeping: ",HousekGene))
        
        plot1 = 
          ggplot(data = PCRstep5,
                 aes(x= as.factor(Time), y = ddCt, col = Sample)) + 
          facet_wrap(~GeneH, ncol = 1) +
          geom_jitter(width = 0.1, height = 0,size = 2)+
          theme_bw()+
          labs(x = "Time", y = "DDCT")
        
        plot2 = ggplot(data = PCRstep5 ,
                       aes(x= as.factor(Time), y = Q, col = Sample)) + 
          facet_wrap(~GeneH, ncol = 1) +
          geom_jitter(width = 0.1, height = 0,size = 2)+
          theme_bw()+
          labs(x = "Time", y = "2^(-DDCT)")
        
        if(is.null(plotList[[paste0(gene,"_H",Hgene)]])){
          
          plotList[[paste0(gene,"_H",Hgene)]] = list(plot = ggplot(), table = NULL)
          
          PCRplotUI =  renderUI({
            plot_output_list <- lapply(names(plotList), function(i) {
              plotOutput(paste0("PCRplot_", i), width = "100%" )
            })
            do.call(tagList, list(plot_output_list))
          })
          PCRtableUI =  renderUI({
            plot_output_list <- lapply(names(plotList), function(i) {
              tablename <- paste("PCRtable_", i, sep="")
              tableOutput(tablename)
            })
            do.call(tagList, list(plot_output_list))
          })
        }
        # Dynamically generate plot output UI
        output$PCRplot <- renderUI({PCRplotUI})
        output$PCRtables <- renderUI({PCRtableUI})
        
        plotList[[paste0(gene,"_H",Hgene)]]$plot = plot1+plot2
        plotList[[paste0(gene,"_H",Hgene)]]$table = PCRstep5 %>% 
          rename(DDCT = ddCt, `2^(-DDCT)` = Q)
        
        output[[paste0("PCRplot_", gene,"_H",Hgene)]] <- renderPlot({
          plotList[[paste0(gene,"_H",Hgene)]]$plot
        })
        output[[paste0("PCRtable_", gene,"_H",Hgene)]] <- renderTable({ plotList[[paste0(gene,"_H",Hgene)]]$table})
        
        pcrResult$plotPCR = plotList
        
        df = do.call(rbind, lapply(plotList, `[[`, 2)) %>% filter(Sample != input$PCRbaseline)
        
        output$PointGenePlot = renderPlot({
          ggplot(df,aes(x = Gene, y = -DDCT)) +
            geom_point() +
            facet_wrap(~HousekGene, ncol = 1)+
            theme_bw()+
            labs(x="",y= "Log2(Q)", title = "")+
            theme(axis.text.x=element_blank(), 
                  axis.ticks.x=element_blank())
        })
        
      }
    })
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
      
      utils::zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
    },
  )
  
  #### END PCR analysis ####
  
  #### ELISA analysis ####
  observeEvent(input$NextElisaQuantif,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "tablesELISA")
  })
  
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
    DataAnalysisModule$elisaResult  = reactiveValuesToList(elisaResult)
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
    tableExcelColored(session = session,
                      Result = elisaResult, 
                      FlagsExp = FlagsELISA,
                      type = "Update")
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent(input$rightTableElisa_cell_edit, {
    info <- input$rightTableElisa_cell_edit
    data <- right_data_elisa() 
    updatedText <- updateTable("right", "ELISA", info, data, elisaResult, FlagsELISA)
    
    output$ELISASelectedValues <- renderText(updatedText)  
    tableExcelColored(session = session,
                      Result = elisaResult, 
                      FlagsExp = FlagsELISA,
                      type = "Update",
    )
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
        
        var1_val <- elisaResult$data$Var1[cellCoo[1]]
        var2_val <- elisaResult$data$Var2[cellCoo[1]]
        elisaResult$data$exp[elisaResult$data$Var1 == var1_val & elisaResult$data$Var2 == var2_val] <- value.now
        
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
      ELISAtb <- elisaResult$TablePlot
      cellCoo <- FlagsELISA$cellCoo
      
      value.bef <- elisaResult$ELISAcell_EXP[cellCoo[1], cellCoo[2]] 
      value.now <- input$ELISAcell_EXP
      
      if (value.now != "" && value.now != value.bef) {
        currentValues <- elisaResult$Initdata[cellCoo[1], cellCoo[2]]
        
        elisaResult$ELISAcell_EXP[cellCoo[1], cellCoo[2]] <- value.now
        
        # Aggiorna data
        elisaResult$data$time[elisaResult$data$Var1 == cellCoo[1] & elisaResult$data$Var2 == cellCoo[2]] <- value.now
        
        elisaResult$Tablestandcurve <- elisaResult$data %>% filter(exp == input$ELISA_standcurve) %>% select(exp, time, values) %>% rename(Concentrations = time,Measures = values)
        
        
        tableExcelColored(session = session,
                          Result = elisaResult, 
                          FlagsExp = FlagsELISA,
                          type = "Update")
        
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
  
  observe({
    elisaResult$data -> data
    
    if(input$ELISA_standcurve != ""){
      print("primo if")
      standcurve <- data %>%
        filter(exp %in% input$ELISA_standcurve) %>%
        select(exp, time, values) %>%
        rename(Concentrations = time, Measures = values)
      
      if(!is.null(elisaResult$Tablestandcurve) && 
         all.equal(elisaResult$Tablestandcurve %>% select(-Concentrations),
                   standcurve %>% select(-Concentrations) ))
      {
        print("secondo if")
        standcurve =  elisaResult$Tablestandcurve
      }else{
        elisaResult$Tablestandcurve = standcurve
      }
      
      print("stampa")
      output$ELISA_Table_stdcurve <- DT::renderDataTable({
        DT::datatable(standcurve,
                      selection = 'none',
                      editable = list(target = "cell",
                                      disable = list(columns = 0:1)),
                      rownames = FALSE
        )
      })
      print("fine stampa")
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
  
  output$downloadElisaAnalysis <- downloadHandler(
    filename = function() {
      paste('ELISAanalysis-', Sys.Date(), '.zip', sep='')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      nomeRDS <- paste0("ELISA_analysis-", Sys.Date(), ".rds")
      nomeXLSX <- paste0("ELISA_analysis-", Sys.Date(), ".xlsx")
      
      tempRdsPath <- file.path(tempDir, nomeRDS)
      tempXlsxPath <- file.path(tempDir, nomeXLSX)
      
      results <- DataAnalysisModule$elisaResult
      saveRDS(results, file = tempRdsPath)
      
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "ELISA")
      
      zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
    } 
  )
  
  # save everytime there is a change in the results
  ELISAresultListen <- reactive({
    reactiveValuesToList(elisaResult)
  })
  observeEvent(ELISAresultListen(), {
    DataAnalysisModule$elisaResult = reactiveValuesToList(elisaResult)
  })
  
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
  
  observeEvent(input$rightTableEndoc_cell_edit, {
    info <- input$rightTableEndoc_cell_edit
    data <- right_data_endoc() 
    updatedText <- updateTable("right", "ENDOC", info, data, endocResult, FlagsENDOC)
    
    tableExcelColored(session = session,
                      Result = endocResult, 
                      FlagsExp = FlagsENDOC,
                      type = "Update")
    output$ENDOCSelectedValues <- renderText(updatedText)  
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  observeEvent(input$leftTableEndoc_cell_edit, {
    info <- input$leftTableEndoc_cell_edit
    data <- left_data_endoc() 
    updatedText <- updateTable("left", "ENDOC", info, data, endocResult, FlagsENDOC)
    
    output$ENDOCSelectedValues <- renderText(updatedText)  
    tableExcelColored(session = session,
                      Result = endocResult, 
                      FlagsExp = FlagsENDOC,
                      type = "Update")
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
      
      results <- DataAnalysisModule$endocResult
      saveRDS(results, file = tempRdsPath)
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "ENDOC")
      
      utils::zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
      
    } 
  )
  
  ### End ENDOC analysis ####
  
  #### CITOXICITY analysis ####
  # next buttons
  observeEvent(input$NextCytotoxQuantif,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "tablesCYTOTOX")
  })
  #
  
  cytotoxResult = reactiveValues(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    CYTOTOXcell_EXP = NULL,
    CYTOTOXcell_REP = NULL,
    CYTOTOXcell_SN = NULL,
    CYTOTOXcell_COLOR = NULL,
    MapBaseline = NULL)
  
  cytotoxResult0 = list(
    Initdata= NULL,
    data = NULL,
    TablePlot = NULL,
    dataFinal = NULL,
    CYTOTOXcell_EXP = NULL,
    CYTOTOXcell_REP = NULL,
    CYTOTOXcell_SN = NULL,
    MapBaseline = NULL)
  
  # save everytime there is a change in the results
  CYTOTOXresultListen <- reactive({
    reactiveValuesToList(cytotoxResult)
  })
  observeEvent(CYTOTOXresultListen(), {
    DataAnalysisModule$cytotoxResult = reactiveValuesToList(cytotoxResult)
    DataAnalysisModule$cytotoxResult$Flags = reactiveValuesToList(FlagsCYTOTOX)
  })
  
  left_data_cytotox <- reactiveVal()
  right_data_cytotox <- reactiveVal()
  
  FlagsCYTOTOX <- reactiveValues(cellCoo = NULL,
                                 AllExp = "",
                                 BASEselected = "",
                                 EXPselected = "",
                                 EXPcol = NULL)
  
  observeEvent(input$LoadCYTOTOX_Button,{
    alert$alertContext <- "CYTOTOX-reset"
    if(!is.null(cytotoxResult$Initdata) ) {
      shinyalert(
        title = "Important message",
        text = "Do you want to update the CYTOTOX data already present, by resetting the previous analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Update",
        cancelButtonText = "Cancel",
      )
    } else loadExcelFileCYTOTOX()
  })
  
  observeEvent(input$shinyalert, {
    removeModal()
    if (input$shinyalert && alert$alertContext == "CYTOTOX-reset") {  
      resetPanel("CYTOTOX", flags = FlagsCYTOTOX, result = cytotoxResult)
      left_data_cytotox <- NULL
      
      loadExcelFileCYTOTOX()
    }
  })
  
  loadExcelFileCYTOTOX <- function() {
    alert$alertContext <- ""
    
    mess = readfile(
      filename = input$CYTOTOXImport$datapath,
      isFileUploaded = !is.null(input$CYTOTOXImport) && file.exists(input$CYTOTOXImport$datapath),
      type = "Excel",
      allDouble = T,
      colname = F,
      colors = T
    )
    
    cytotoxResult$Initdata = mess$x
    FlagsCYTOTOX$EXPcol = mess$fill
    cytotoxResult$CYTOTOXcell_COLOR = mess$SNtable
    cytotoxResult$CYTOTOXcell_SN <- matrix("", nrow = nrow(cytotoxResult$CYTOTOXcell_COLOR), ncol = ncol(cytotoxResult$CYTOTOXcell_COLOR))
    
    showAlert("Success", "The CYTOTOX excel has been uploaded with success", "success", 2000)
  }
  
  observe({
    if( !is.null(cytotoxResult$Initdata) && is.null(cytotoxResult$TablePlot) ){
      tableExcelColored(session = session,
                        Result = cytotoxResult, 
                        FlagsExp = FlagsCYTOTOX,
                        type = "Initialize")
      
      output$CYTOTOXmatrix <-renderDataTable({cytotoxResult$TablePlot})
    }
  })
  
  observe({
    color_codes <- FlagsCYTOTOX$EXPcol
    color_names <- names(FlagsCYTOTOX$EXPcol)
    
    valid_colors <- color_codes != "white"
    color_codes <- color_codes[valid_colors]
    color_names <- color_names[valid_colors]
    
    mid_point <- ceiling(length(color_codes) / 2)
    left_colors <- color_codes[1:length(color_codes)]
    right_colors <- color_codes[(mid_point+1):length(color_codes)]
    
    left_formatted_data <- get_formatted_data(left_colors, color_names[1:length(color_codes)], cytotoxResult, cytotoxResult$CYTOTOXcell_EXP, "CYTOTOX")
    #right_formatted_data <- get_formatted_data(right_colors, color_names[(mid_point+1):length(color_codes)], cytotoxResult, cytotoxResult$CYTOTOXcell_EXP, "CYTOTOX")
    
    left_data_cytotox(left_formatted_data)
    #right_data_cytotox(right_formatted_data)
    
    output$leftTableCytotox <- renderDataTable(
      left_data_cytotox(), 
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
          list(width = '180px', targets = 3),
          list(width = '150px', targets = 4),
          list(width = '200px', targets = 5),
          list(width = '200px', targets = 6),
          list(className = 'dt-head-left dt-body-left', targets = 1)
        )
      )
    )
  })
  
  observeEvent(input$leftTableCytotox_cell_edit, {
    info <- input$leftTableCytotox_cell_edit
    data <- left_data_cytotox() 
    updatedText <- updateTable("left", "CYTOTOX", info, data, cytotoxResult, FlagsCYTOTOX)
    
    output$CYTOTOXSelectedValues <- renderText(updatedText)  
    
    tableExcelColored(session = session,
                      Result = cytotoxResult, 
                      FlagsExp = FlagsCYTOTOX,
                      type = "Update",
    )
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent(input$CYTOTOXmatrix_cell_clicked,{
    req(input$CYTOTOXmatrix_cell_clicked)  
    
    cellSelected= as.numeric(input$CYTOTOXmatrix_cell_clicked)
    FlagsCYTOTOX$cellCoo = cellCoo = c(cellSelected[1],cellSelected[2] + 1)
    
    allExp <- unique(na.omit(c(cytotoxResult$CYTOTOXcell_EXP)))  
    selectedExp <- ifelse(is.null(cytotoxResult$CYTOTOXcell_EXP[cellCoo[1], cellCoo[2]]), "", cytotoxResult$CYTOTOXcell_EXP[cellCoo[1], cellCoo[2]])
    
    updateSelectizeInput(inputId = "CYTOTOXcell_EXP",
                         choices = allExp,
                         selected = selectedExp)
    
    allSN <- unique(na.omit(c(cytotoxResult$CYTOTOXcell_SN)))  
    selectedSN <- ifelse(is.null(cytotoxResult$CYTOTOXcell_SN[cellCoo[1], cellCoo[2]]), "", cytotoxResult$CYTOTOXcell_SN[cellCoo[1], cellCoo[2]])
    
    updateSelectizeInput(inputId = "CYTOTOXcell_SN",
                         choices = allSN,
                         selected = selectedSN)
    
    allREP <- unique(na.omit(c(cytotoxResult$CYTOTOXcell_REP)))  
    selectedREP <- ifelse(is.null(cytotoxResult$CYTOTOXcell_REP[cellCoo[1], cellCoo[2]]), "", cytotoxResult$CYTOTOXcell_REP[cellCoo[1], cellCoo[2]])
    
    updateSelectizeInput(inputId = "CYTOTOXcell_REP",
                         choices = allREP,
                         selected = selectedREP)
  })
  
  observeEvent(input$CYTOTOXcell_EXP,{
    if (!is.null(cytotoxResult$CYTOTOXcell_EXP) && !is.null(FlagsCYTOTOX$cellCoo) && !anyNA(FlagsCYTOTOX$cellCoo)) {
      CYTOTOXtb = cytotoxResult$TablePlot
      cellCoo = FlagsCYTOTOX$cellCoo
      
      value.bef = cytotoxResult$CYTOTOXcell_EXP[cellCoo[1], cellCoo[2]] 
      value.now = input$CYTOTOXcell_EXP
      
      if (value.now != "" && value.now != value.bef) {
        currentValues <- cytotoxResult$Initdata[cellCoo[1], cellCoo[2]]
        
        cytotoxResult$CYTOTOXcell_EXP[cellCoo[1], cellCoo[2]] = value.now
        tableExcelColored(session = session,
                          Result = cytotoxResult, 
                          FlagsExp = FlagsCYTOTOX,
                          type = "Update",
        )
        
        output$CYTOTOXSelectedValues <- renderText(paste("Updated value", paste(currentValues), ": time ", value.now))
        output$CYTOTOXmatrix <- renderDataTable({cytotoxResult$TablePlot})
      }
    }  else return()
  }, ignoreInit = TRUE)
  
  observeEvent(input$CYTOTOXcell_REP,{
    if (!is.null(cytotoxResult$CYTOTOXcell_REP) && !is.null(FlagsCYTOTOX$cellCoo) && !anyNA(FlagsCYTOTOX$cellCoo)) {
      CYTOTOXtb = cytotoxResult$TablePlot
      cellCoo = FlagsCYTOTOX$cellCoo
      
      value.bef = cytotoxResult$CYTOTOXcell_REP[cellCoo[1], cellCoo[2]] 
      value.now = input$CYTOTOXcell_REP
      
      if (value.now != "" && value.now != value.bef) {
        currentValues <- cytotoxResult$Initdata[cellCoo[1], cellCoo[2]]
        cytotoxResult$CYTOTOXcell_REP[cellCoo[1], cellCoo[2]] = value.now
        
        tableExcelColored(session = session,
                          Result = cytotoxResult, 
                          FlagsExp = FlagsCYTOTOX,
                          type = "Update",
        )
        
        output$CYTOTOXSelectedValues <- renderText(paste("Updated value", paste(currentValues), ": time ", value.now))
        output$CYTOTOXmatrix <- renderDataTable({cytotoxResult$TablePlot})    
      }
    }else return()
  }, ignoreInit = TRUE)
  
  observeEvent(input$CYTOTOXcell_SN,{
    if (!is.null(cytotoxResult$CYTOTOXcell_COLOR) && !is.null(FlagsCYTOTOX$cellCoo) && !anyNA(FlagsCYTOTOX$cellCoo)) {
      CYTOTOXtb = cytotoxResult$TablePlot
      cellCoo = FlagsCYTOTOX$cellCoo
      
      value.bef = cytotoxResult$CYTOTOXcell_SN[cellCoo[1],cellCoo[2]] 
      value.now = input$CYTOTOXcell_SN
      
      if (value.now != "" && value.now != value.bef) {
        currentValues <- cytotoxResult$Initdata[cellCoo[1], cellCoo[2]]
        
        cytotoxResult$CYTOTOXcell_COLOR[cellCoo[1], cellCoo[2]] = value.now
        cytotoxResult$CYTOTOXcell_SN[cellCoo[1],cellCoo[2]] = value.now
        CYTOTOXtb$x$data[cellCoo[1],paste0("Col",cellCoo[2])] = value.now
        
        if(! input$CYTOTOXcell_SN %in% FlagsCYTOTOX$AllExp){
          exp = unique(c(FlagsCYTOTOX$AllExp,input$CYTOTOXcell_SN))
          FlagsCYTOTOX$AllExp = exp
        }
        
        ## updating table and colors definition depending on the cell fill 
        tableExcelColored(session = session,
                          Result = cytotoxResult, 
                          FlagsExp = FlagsCYTOTOX,
                          type = "Update")
        
        output$CYTOTOXSelectedValues <- renderText(paste("Updated value", paste(currentValues), ": sample name ", value.now))
        output$CYTOTOXmatrix <-renderDataTable({cytotoxResult$TablePlot})
      }
    } else return()
  }, ignoreInit = TRUE)
  
  ## update Baselines checkBox
  observeEvent(FlagsCYTOTOX$AllExp,{
    if(length(FlagsCYTOTOX$AllExp) > 1){
      exp = FlagsCYTOTOX$AllExp
      exp = exp[exp != ""]
      
      exp_selec = input$CYTOTOX_baselines
      
      updateSelectizeInput(session,"CYTOTOX_baselines",
                           choices = exp,
                           selected = FlagsCYTOTOX$BASEselected )
    }
  })
  
  ## select the baselines
  observeEvent(input$CYTOTOX_baselines,{
    FlagsCYTOTOX$BASEselected = input$CYTOTOX_baselines
    FlagsCYTOTOX$EXPselected = FlagsCYTOTOX$AllExp[! FlagsCYTOTOX$AllExp %in% FlagsCYTOTOX$BASEselected]
  },ignoreNULL = F)
  
  toListen_cytotox <- reactive({
    return( list(cytotoxResult$CYTOTOXcell_EXP,cytotoxResult$CYTOTOXcell_SN, cytotoxResult$CYTOTOXcell_REP, FlagsCYTOTOX$BASEselected) )
  })
  observeEvent(toListen_cytotox(),{
    baselines = FlagsCYTOTOX$BASEselected
    baselines = baselines[baselines != ""]
    if(length(baselines) > 0 )
    {
      CYTOTOXcell_value = data.frame(
        row = c(t(row(cytotoxResult$Initdata))),
        col = c(t(col(cytotoxResult$Initdata))),
        Val = c(t(cytotoxResult$Initdata))
      ) 
      CYTOTOXcell_SN = data.frame(
        row = c(t(row(cytotoxResult$CYTOTOXcell_SN))),
        col = c(t(col(cytotoxResult$CYTOTOXcell_SN))),
        SN = c(t(cytotoxResult$CYTOTOXcell_SN))
      ) 
      CYTOTOXcell_EXP = data.frame(
        row = c(t(row(cytotoxResult$CYTOTOXcell_EXP))),
        col = c(t(col(cytotoxResult$CYTOTOXcell_EXP))),
        EXP = c(t(cytotoxResult$CYTOTOXcell_EXP))
      ) 
      CYTOTOXcell_REP = data.frame(
        row = c(t(row(cytotoxResult$CYTOTOXcell_REP))),
        col = c(t(col(cytotoxResult$CYTOTOXcell_REP))),
        REP = c(t(cytotoxResult$CYTOTOXcell_REP))
      ) 
      
      CYTOTOXcell = merge(merge(merge(CYTOTOXcell_REP,CYTOTOXcell_EXP) , CYTOTOXcell_SN), CYTOTOXcell_value) %>%
        group_by(SN,EXP,REP) %>%
        summarise(MeanV = mean(Val,na.rm = T)) %>%
        ungroup()
      
      CYTOTOXcell_base = CYTOTOXcell %>%
        filter(SN == baselines) %>%
        rename(MeanBaseV = MeanV, Baseline = SN) 
      # if no exp conditions is used for the baseline then we repeat  exp given for the other SN
      # in thius way the same baseline is used for each SN
      if(all(CYTOTOXcell_base$EXP == "") ) {
        CYTOTOXcell_base = do.call("rbind", lapply(unique(CYTOTOXcell$EXP), function(x) CYTOTOXcell_base %>% mutate(EXP = x) ) )
      }
      
      CYTOTOXcell = CYTOTOXcell %>%  filter(SN != baselines)
      CYTOTOXcell = merge(CYTOTOXcell_base,CYTOTOXcell,by= c("REP","EXP"),all.y = T)
      CYTOTOXcell = CYTOTOXcell %>% mutate(Res = (MeanV-MeanBaseV)/(100-MeanBaseV)*100)
      
      cytotoxResult$data =  CYTOTOXcell
      cytotoxResult$dataFinal =  CYTOTOXcell %>%
        select(-Baseline, - MeanBaseV, - MeanV ) %>%
        tidyr::spread(key= "REP", value = "Res") %>%
        rename(`Sample Name` = SN,
               `Experimental Condition` = EXP)
      
      output$CYTOTOXtables = renderDT({
        datatable(
          cytotoxResult$dataFinal,
          rownames= FALSE,
          options = list(
            scrollX = TRUE,
            lengthChange = FALSE,
            dom = 't')
        )
      })
      
      output$CYTOTOXplots = renderPlot({
        pl1 = CYTOTOXcell %>%
          group_by(EXP,SN) %>%
          summarize(Mean = mean(MeanV),SD = sd(MeanV)) %>%
          ggplot(aes(x = as.factor(EXP),y = Mean ,fill = SN, col = SN)) + 
          geom_bar(stat="identity", color="black", position=position_dodge()) +
          geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                        position=position_dodge(.9)) +
          theme_bw()+ 
          theme(legend.position = "bottom")+
          labs(title = "Sample Name mean values with standard deviation bars",
               col="Sample Name",fill="Sample Name",
               x = "Experimental condition", y= "Mean Values")
        
        pl2 = CYTOTOXcell %>% ggplot() +
          geom_boxplot(aes(x = as.factor(EXP),y = Res,fill = SN, col = SN),alpha = 0.4) +
          theme_bw() + theme(legend.position = "bottom") +
          labs(x = "Experimental condition", y= "% Values w.r.t \nthe baseline cell death",
               col="Sample Name",fill="Sample Name")
        pl1+pl2
      })
      
    }
  })
  
  output$downloadCYTOTOXAnalysis <- downloadHandler(
    filename = function() {
      paste('CYTOTOXanalysis-', Sys.Date(), '.zip', sep='')
    },
    content = function(file) {
      manageSpinner(TRUE)
      
      tempDir <- tempdir()
      nomeRDS <- paste0("CYTOTOX_analysis-", Sys.Date(), ".rds")
      nomeXLSX <- paste0("CYTOTOX_analysis-", Sys.Date(), ".xlsx")
      
      tempRdsPath <- file.path(tempDir, nomeRDS)
      tempXlsxPath <- file.path(tempDir, nomeXLSX)
      
      results <- DataAnalysisModule$cytotoxResult
      saveRDS(results, file = tempRdsPath)
      
      saveExcel(filename = tempXlsxPath, ResultList=results, analysis = "CYTOTOX")
      
      zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
    } 
  )
  
  ### End CITOXICITY analysis ####
  
  #### FACS analysis ####
  facsResult = reactiveValues(
    Initdata= NULL,
    data = NULL,
    dataFinal = NULL,
    depth = NULL,
    depthCount = NULL,
    originalName = NULL,
    originalColumnName = NULL,
    name = NULL,
    columnName = NULL,
    statistics = NULL,
    cells = NULL,
    ExpConditionDF = NULL,
    Barplot = NULL,
    StatDF = NULL
  )
  
  facsResult0 = reactiveValues(
    Initdata= NULL,
    data = NULL,
    dataFinal = NULL,
    depth = NULL,
    depthCount = NULL,
    originalName = NULL,
    originalColumnName = NULL,
    name = NULL,
    columnName = NULL,
    statistics = NULL,
    cells = NULL,
    ExpConditionDF = NULL,
    Barplot = NULL,
    StatDF = NULL
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
    actualPath = NULL,
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
        facsResult$statistics <- vector("list", nrow(data))
        facsResult$cells <- vector("list", nrow(data))
        facsResult$originalName <- vector("list", nrow(data))
        
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
        showAlert("Success", "The Excel has been uploaded with success", "success", 2000)
        
        output$dynamicSelectize <- renderUI({
          updateSelectizeUI(maxDepth)
        })
        
        updateTabsetPanel(session, "SideTabs", selected = "tablesFACS")
        
      }
    }
  }
  
  observeEvent(input$SideTabs, {
    if (input$SideTabs == "tablesFACS") {
      FlagsFACS$actualLevel <- 0
    }
  }, ignoreInit = TRUE)  
  
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
      
      ## update the table for the barplot in which an expcondition is associated to each name
      ExpConditionDF <- facsResult$ExpConditionDF <- data.frame(Name = unique(unlist(filtered_names)), ExpCondition = "")
      
      output$FACSexpcond_tab = renderDT({
        DT::datatable( ExpConditionDF,
                       selection = 'none',
                       editable = list(target = "cell",
                                       disable = list(columns = 0) ),
                       rownames= FALSE,
                       options = list(scrollX = TRUE,
                                      searching = FALSE,
                                      dom = 't' # Only display the table
                       )
        )
      })
      ######
      data_for_table <- data.frame(
        Name = unlist(filtered_names),
        Start = "100%",
        stringsAsFactors = FALSE
      )
      
      facsResult$originalName <- unlist(filtered_names)
      FlagsFACS$data <- data_for_table  
      FlagsFACS$actualPath <- ""
      
      output$FACSmatrix <- renderDT({
        datatable(data_for_table, options = list(
          pageLength = 10,
          autoWidth = TRUE
        ))
      })
      
      data_for_name_update <- data.frame(
        Name = unlist(filtered_names),
        New_name = rep("-", length(filtered_names)),  
        stringsAsFactors = FALSE
      )
      
      output$FACSnameUpdate <- renderDT({
        datatable(data_for_name_update, options = list(
          pageLength = 10,
          autoWidth = TRUE,
          columnDefs = list(
            list(targets = 1, width = '50%'),  
            list(targets = 2, width = '50%')
          )
        ), editable = list(target = 'cell', columns = 2))  
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
    loadDrop(facsResult, FlagsFACS, session) 
  }, ignoreInit = TRUE)
  
  observe({
    FlagsFACS$actualPath
    choices <- colnames(FlagsFACS$data)
    
    if (length(choices) > 2) {
      choices <- choices[3:length(choices) - 1]  
    } else {
      choices <- character(0) 
    }
    
    updateSelectizeInput(session = session,
                         inputId = "selectBaseGate",
                         choices = choices
    )
  })
  
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
    
    data.frame(Name = data_row[1], FinalPercent = sprintf("%.2f%%", final_percentage)) %>%
      setNames(c("Name", last_percentage_name))
  }
  
  observeEvent(input$SaveFACSanalysis, {
    if (is.null(input$selectBaseGate) || input$selectBaseGate == "") {
      showAlert("Error", "No gate selected. Please select a gate", "error", 5000)
      return()
    }
    
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
    
    current_names <- FlagsFACS$data$Name
    if (!is.null(facsResult$dataFinal)) {
      facsResult$dataFinal$Name <- current_names
    }
    
    new_column_name <- names(new_data)[2]
    if (is.null(facsResult$dataFinal)) {
      facsResult$dataFinal <- new_data
    } else {
      if (new_column_name %in% names(facsResult$dataFinal)) {
        facsResult$dataFinal[[new_column_name]] <- new_data[, 2]
      } else {
        facsResult$dataFinal <- cbind(facsResult$dataFinal, new_data[, 2])
        names(facsResult$dataFinal)[ncol(facsResult$dataFinal)] <- new_column_name
      }
    }
    
    existing_column_names <- names(facsResult$dataFinal)
    for (col_name in existing_column_names) {
      if (col_name != new_column_name && all(facsResult$dataFinal[[col_name]] == facsResult$dataFinal[[new_column_name]])) {
        facsResult$dataFinal[[col_name]] <- NULL
      }
    }
    
    if (!is.null(facsResult$originalColumnName)) {
      facsResult$originalColumnName <- unique(c(facsResult$originalColumnName, new_column_name))
    } else {
      facsResult$originalColumnName <- c(new_column_name)
    }
    
    output$FACSresult <- renderDT({
      datatable(facsResult$dataFinal, 
                options = list(
                  autoWidth = TRUE,
                  columnDefs = list(
                    list(visible = FALSE, targets = 0)  
                  )
                )
      )
    })
    
    column_names <- colnames(facsResult$dataFinal)
    column_names <- column_names[column_names != "Name"]
    
    data_for_column_update <- data.frame(
      Name = column_names,
      New_name = rep("-", length(column_names)),
      stringsAsFactors = FALSE
    )
    
    facsResult$columnName <- data_for_column_update
    
    output$FACScolumnNameUpdate <- renderDT({
      datatable(data_for_column_update, options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 1, width = '50%'),  
          list(targets = 2, width = '50%')
        )
      ), editable = list(target = 'cell', columns = 2))
    })
  })
  
  observeEvent(input$FACSexpcond_tab_cell_edit, {
    info <- input$FACSexpcond_tab_cell_edit
    str(info)  # Debugging line to show what 'info' contains
    info$col = 2
    # Update the data in the reactive value
    facsResult$ExpConditionDF <- editData(facsResult$ExpConditionDF, info)
  })
  
  observe({
    facsResult$ExpConditionDF -> ExpConditionDF
    facsResult$dataFinal -> dataFinal
    
    isolate({
      if(!is.null(dataFinal) && !is.null(ExpConditionDF)  ){
        
        output$FACSexpcond_tab = renderDT({
          DT::datatable( ExpConditionDF,
                         selection = 'none',
                         editable = list(target = "cell",
                                         disable = list(columns = 0) ),
                         rownames= FALSE,
                         options = list(scrollX = TRUE,
                                        searching = FALSE,
                                        dom = 't' # Only display the table
                         )
          )
        })
        
        dataFinalExpCond = merge(dataFinal,ExpConditionDF) %>%
          tidyr::gather(-Name,-ExpCondition, value = "Percetages", key = "Gate" ) %>%
          mutate(Percetages = as.numeric(gsub(replacement = "",x = Percetages,pattern = "%")))
        
        pl = ggplot(dataFinalExpCond, aes(x =  ExpCondition , y =  Percetages, color = Name, fill = Name)) + 
          geom_bar(stat = "identity", position = "dodge") +
          facet_wrap(~Gate) + theme_bw() + labs(y = "Percetages %", x = "Experimetal condition")
        
        plot1 = lapply(unique(dataFinalExpCond$Gate),function(g){
          ggplot(dataFinalExpCond %>% filter(Gate == g),
                 aes(x =  ExpCondition , y =  Percetages, color = Name, fill = Name)) + 
            geom_bar(stat = "identity", position = "dodge") +
            facet_wrap(~Gate) +
            theme_bw() +
            labs(y = "Percetages %", x = "Experimetal condition")
        })
        
        
        # Dynamically generate plot output UI
        output$FACSexpcond_plot <- renderUI({
          plot_output_list_1 <- lapply(1:length(plot1), function(i) {
            plotOutput(paste0("FACSplot_", i), width = "100%" )
          })
          do.call(tagList, list(plot_output_list_1))
        })
        
        lapply(1:length(plot1), function(i) {
          output[[paste0("FACSplot_", i)]] <- renderPlot({
            plot1[[i]]
          })
        })
        
        facsResult$Barplot = pl
      }
    })
  })
  
  observeEvent(input$FACSstatButton, {
    facsResult$ExpConditionDF -> ExpConditionDF
    facsResult$dataFinal -> dataFinal
    isolate({
      if(!is.null(dataFinal) && !is.null(ExpConditionDF)  ){
        dataFinalExpCond = merge(dataFinal,ExpConditionDF) %>%
          tidyr::gather(-Name,-ExpCondition, value = "Percetages", key = "Gate" ) %>%
          mutate(Percetages = as.numeric(gsub(replacement = "",x = Percetages,pattern = "%")))
        
        ### Ttest
        if(length(unique(dataFinalExpCond$ExpCondition)) > 1){
          StatDF = do.call( rbind, 
                            lapply(unique(dataFinalExpCond$Gate), function(gate){
                              subdata = dataFinalExpCond %>% filter(Gate == gate) %>% select(ExpCondition,Percetages)
                              testStat.function(tibble(subdata), var = gate)
                            }) )
        }
      }
      
    })
    
    facsResult$StatDF = StatDF
    output$FACSstat_tab = renderDT({
      DT::datatable( StatDF,
                     selection = 'none',
                     rownames= FALSE,
                     options = list(scrollX = TRUE,
                                    searching = FALSE,
                                    dom = 't' # Only display the table
                     )
      )
    })
    
  })
  
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
        
      }
    }
  })
  
  observeEvent(facsResult$name, {
    valid_indices <- facsResult$depthCount == FlagsFACS$actualLevel
    filtered_names <- facsResult$name[valid_indices]
    
    unique_names <- unique(sapply(filtered_names, function(name) strsplit(name, "/")[[1]][1]))
    FlagsFACS$data$Name <- unique_names    
    
    proxy <- dataTableProxy('FACSmatrix')
    replaceData(proxy, FlagsFACS$data, resetPaging = FALSE)
    
    if(!is.null(facsResult$dataFinal)){
      facsResult$dataFinal$Name <- unique_names
      proxy <- dataTableProxy('FACSresult')
      replaceData(proxy, facsResult$dataFinal, resetPaging = FALSE)
    }
    
    if(!is.null(facsResult$ExpConditionDF)){
      facsResult$ExpConditionDF$Name <- unique_names
      #proxy <- dataTableProxy('FACSexpcond_tab')
      #replaceData(proxy, facsResult$ExpConditionDF, resetPaging = FALSE)
    }
    
  }, ignoreInit = TRUE)
  
  observeEvent(input$FACScolumnNameUpdate_cell_edit, {
    info <- input$FACScolumnNameUpdate_cell_edit
    
    row <- info$row
    col <- info$col
    new_value <- info$value
    
    old_name <- facsResult$columnName[row, "Name"]
    original_name <- facsResult$originalColumnName[row]
    
    if (new_value == "") {
      old_value <- original_name
      new_value <- original_name
    }
    
    if (!is.null(facsResult$dataFinal)) {
      colnames(facsResult$dataFinal)[colnames(facsResult$dataFinal) == old_name] <- new_value
    }
    
    if (new_value == original_name) {
      facsResult$columnName[row, "Name"] <- original_name
      facsResult$columnName[row, "New_name"] <- "-"
    } else {
      facsResult$columnName[row, "Name"] <- new_value
      facsResult$columnName[row, "New_name"] <- "-"
    }
    
    output$FACSresult <- renderDT({
      datatable(facsResult$dataFinal, 
                options = list(
                  autoWidth = TRUE,
                  columnDefs = list(
                    list(visible = FALSE, targets = 0)  
                  )
                )
      )
    })
    
    column_names <- colnames(facsResult$dataFinal)
    column_names <- column_names[column_names != "Name"]
    
    data_for_column_update <- data.frame(
      Name = facsResult$originalColumnName,
      New_name = sapply(1:nrow(facsResult$columnName), function(i) {
        if (facsResult$columnName$Name[i] == facsResult$originalColumnName[i]) {
          return("-")
        } else {
          return(facsResult$columnName$Name[i])
        }
      }),
      stringsAsFactors = FALSE
    )
    
    #facsResult$columnName <- data_for_column_update
    
    output$FACScolumnNameUpdate <- renderDT({
      datatable(data_for_column_update, options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 1, width = '50%'),  
          list(targets = 2, width = '50%')
        )
      ), editable = list(target = 'cell', columns = 2))
    })
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
      
      utils::zip(file, files = c(tempRdsPath, tempXlsxPath), flags = "-j")
      manageSpinner(FALSE)
    }
  )
  
  ### End FACS analysis ####
  
  ### Start Statistic  ####
  
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
      names(mess) -> namesRes
      if("Flags"%in% namesRes) namesRes = namesRes[ namesRes != "Flags"]
      
      if(all(namesRes %in% names(wbquantResult)) || all(namesRes %in% names(DataAnalysisModule))){
        DataStatisticModule$WB[[dpath]] <- mess$AdjRelDensity %>% mutate(DataSet = dpath)
      } else if(all(namesRes %in% names(pcrResult)) || all(namesRes %in% names(DataAnalysisModule))){
        DataStatisticModule$PCR[[dpath]]  <- mess
      } else if(all(namesRes %in% names(endocResult)) || all(namesRes %in% names(DataAnalysisModule))){
        DataStatisticModule$ENDOC[[dpath]]  <- mess
      } else if(all(namesRes %in% names(elisaResult)) || all(namesRes %in% names(DataAnalysisModule))){
        DataStatisticModule$ELISA[[dpath]]  <- mess
      } else if(all(namesRes %in% names(cytotoxResult)) || all(namesRes %in% names(DataAnalysisModule))){
        DataStatisticModule$CYTOTOX[[dpath]]  <- mess
      } else if(all(namesRes %in% names(ifResult)) || all(namesRes %in% names(DataAnalysisModule))){
        DataStatisticModule$IF[[dpath]]  <- mess
      } else if(all(namesRes %in% names(facsResult)) || all(namesRes %in% names(DataAnalysisModule))){
        DataStatisticModule$FACS[[dpath]]  <- mess
      }else{
        showAlert("Error", paste(mess[["message"]],"\n The file must be RDs saved through the Data Analysis module."), "error", 5000)
        manageSpinner(FALSE)
        return()
      }
      DataStatisticModule$Flag <- TRUE
    }
    manageSpinner(FALSE)
    showAlert("Success", "The RDs files have been uploaded with success", "success", 2000)
    return(NULL)
  })
  
  StatisticalAnalysisResults <- reactive({
    if (!is.null(input$StatAnalysis) && input$StatAnalysis != "") {
      results <- DataStatisticModule[[input$StatAnalysis]]
      resTTest <- NULL
      resMainTest <- NULL
      resPairwise <- NULL
      steps <- NULL
      resplot <- NULL
      stats <- NULL
      
      switch(input$StatAnalysis, 
             "WB" = {
               do.call(rbind, results) -> results
               points <- results %>%
                 select(DataSet, SampleName, AdjRelDensity) %>%
                 mutate(SampleName = gsub(pattern = "^[0-9]. ", x = SampleName, replacement = ""),
                        SampleName = trimws(SampleName),
                        ColorSet = as.character(DataSet))
               
               stats <- points %>%
                 group_by(SampleName) %>%
                 summarise(Mean = mean(AdjRelDensity), sd = sd(AdjRelDensity), .groups = 'drop')
               
               res <- testStat.function(points %>% select(SampleName, AdjRelDensity) %>% as_tibble())
               resTTest <- res$resTTest
               resMainTest <- res$test
               resPairwise <- res$pairwise
               steps <- res$steps
               
               main_test_pvalue <- if (!is.null(resMainTest)) resMainTest$pValue else NA
               main_test_txt <- if (is.null(resPairwise)) {
                 signif(main_test_pvalue, digits = 3) %>% paste0("P value from main test = ", .)
               } else {
                 ""
               }
               
               resplot <- ggplot(stats, aes(x = SampleName, y = Mean)) + 
                 geom_bar(stat="identity", color="black", fill = "#BAE1FF", position=position_dodge()) +
                 geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(.9)) +
                 geom_point(data = points, aes(x = SampleName, y = AdjRelDensity, color = ColorSet),
                            position = position_jitter(width = 0.2), size = 3) +
                 theme_bw() +
                 labs(title = "Results", subtitle = main_test_txt, y = "AdjRelDensity", x = "SampleName", color = "Sample Name") +
                 annotate("text", x = Inf, y = Inf, label = "ns: p > 0.05\n*: p <= 0.05\n**: p <= 0.01\n ***: p <= 0.001", 
                          hjust = 1.1, vjust = 1.5, size = 5, color = "black")
             },
             "IF" = {
               resultsNew <- do.call(rbind,
                                     lapply(1:length(results),
                                            function(l){
                                              d = results[[l]]$SubStatData
                                              d$File = l
                                              d
                                            } 
                                     ) 
               )
               
               stats <- resultsNew %>%
                 group_by(ExpCond) %>% 
                 summarise(Mean = mean(Values), sd = sd(Values))
               
               res <- testStat.function(resultsNew[,c("ExpCond", "Values")])
               resTTest <- res$resTTest
               resMainTest <- res$test
               resPairwise <- res$pairwise
               steps <- res$steps
               
               main_test_pvalue <- if (!is.null(resMainTest)) resMainTest$pValue else NA
               main_test_txt <- if (is.null(resPairwise)) {
                 signif(main_test_pvalue, digits = 3) %>% paste0("P value from main test = ", .)
               } else {
                 ""
               }
               
               resplot <- ggplot(stats, aes(x = ExpCond, y = Mean)) + 
                 geom_bar(stat="identity", color="black", fill = "#BAE1FF", position=position_dodge()) +
                 geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(.9)) +
                 geom_point(data = resultsNew, aes(x = ExpCond, y = Values, color = as.factor(File)),
                            position = position_jitter(width = 0.2), size = 3) +
                 theme_bw()+
                 labs(title = "Results", subtitle = main_test_txt, color = "File") +
                 annotate("text", x = Inf, y = Inf, label = "ns: p > 0.05\n*: p <= 0.05\n**: p <= 0.01\n ***: p <= 0.001", 
                          hjust = 1.1, vjust = 1.5, size = 5, color = "black")
             },
             "PCR" = {
               resultsNew <- do.call(rbind,
                                     lapply(1:length(results),
                                            function(l){
                                              d = results[[l]]
                                              d$File = l
                                              d
                                            } 
                                     ) 
               )
               
               resultsNew <- resultsNew %>% mutate(GeneH = paste(Gene, ", Housekeeping: ", HousekGene))
               stats <- resultsNew %>% group_by(GeneH) %>% summarise(Mean = mean(Q), sd = sd(Q))
               
               res <- testStat.function(resultsNew[, c("GeneH", "Q")])
               resTTest <- res$resTTest
               resMainTest <- res$test
               resPairwise <- res$pairwise
               steps <- res$steps
               
               main_test_pvalue <- if (!is.null(resMainTest)) resMainTest$pValue else NA
               main_test_txt <- if (is.null(resPairwise)) {
                 signif(main_test_pvalue, digits = 3) %>% paste0("P value from main test = ", .)
               } else {
                 ""
               }
               
               resplot <- ggplot(stats, aes(x = GeneH, y = Mean)) + 
                 geom_bar(stat="identity", color="black", fill = "#BAE1FF", position=position_dodge()) +
                 geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(.9)) +
                 geom_point(data = resultsNew, aes(x = GeneH, y = Q, color = as.factor(File)),
                            position = position_jitter(width = 0.2), size = 3) +
                 theme_bw() +
                 labs(title = "Results", subtitle = main_test_txt, y = "Q", x = "GeneH", color = "File") +
                 annotate("text", x = Inf, y = Inf, label = "ns: p > 0.05\n*: p <= 0.05\n**: p <= 0.01\n ***: p <= 0.001", 
                          hjust = 1.1, vjust = 1.5, size = 5, color = "black")
             },
             "FACS" = {
               resultsNew <- do.call(rbind,
                                     lapply(1:length(results), function(l) {
                                       d <- results[[l]]$dataFinal
                                       d <- d %>%
                                         tidyr::gather(-Name, key = "Gate", value = "Perc") %>%
                                         mutate(Perc = as.numeric(gsub("%", "", Perc)))
                                       d$File <- l
                                       d <- merge(d, results[[l]]$ExpConditionDF)
                                       d
                                     })
               )
               
               stats <- resultsNew %>%
                 group_by(ExpCondition, Gate) %>%
                 summarise(Mean = mean(Perc), sd = sd(Perc), .groups = 'drop')
               
               resList <- lapply(unique(resultsNew$Gate), function(g) {
                 gate_results <- resultsNew %>% filter(Gate == g)
                 testStat.function(gate_results[, c("ExpCondition", "Perc")])
               })
               
               resTTest <- do.call(rbind, lapply(resList, function(res) res$resTTest))
               resMainTest <- do.call(rbind, lapply(resList, function(res) res$test))
               resPairwise <- do.call(rbind, lapply(resList, function(res) res$pairwise))
               steps <- unlist(lapply(resList, function(res) res$steps))
               
               main_test_pvalue <- if (!is.null(resMainTest)) resMainTest$pValue else NA
               main_test_txt <- if (is.null(resPairwise)) {
                 signif(main_test_pvalue, digits = 3) %>% paste0("P value from main test = ", .)
               } else {
                 ""
               }
               
               resplot <- ggplot(stats, aes(x = ExpCondition, y = Mean)) + 
                 geom_bar(stat = "identity", color = "black", fill = "#BAE1FF", position = position_dodge()) +
                 geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width = .2, position = position_dodge(.9)) +
                 geom_point(data = resultsNew, aes(x = ExpCondition, y = Perc, color = as.factor(File)),
                            position = position_jitter(width = 0.2), size = 3) +
                 theme_bw() +
                 labs(title = "Results", subtitle = main_test_txt, color = "File") +
                 facet_wrap(~ Gate) +
                 annotate("text", x = Inf, y = Inf, label = "ns: p > 0.05\n*: p <= 0.05\n**: p <= 0.01\n ***: p <= 0.001", 
                          hjust = 1.1, vjust = 1.5, size = 5, color = "black")
             }
      )
      
      if (!is.null(resPairwise) && main_test_pvalue < 0.05) {
        resPairwise <- resPairwise %>%
          separate(Condition, into = c("group1", "group2"), sep = " vs ") %>%
          mutate(group1 = trimws(group1), 
                 group2 = trimws(group2),
                 stars = cut(pValue, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns")))
        
        all_combinations <- combn(unique(stats$SampleName), 2, simplify = FALSE)
        
        post_hoc_pairs <- lapply(all_combinations, function(pair) {
          pair <- unlist(pair)
          c(trimws(pair[1]), trimws(pair[2]))
        })
        
        annotations <- sapply(post_hoc_pairs, function(pair) {
          idx <- which((resPairwise$group1 == pair[1] & resPairwise$group2 == pair[2]) |
                         (resPairwise$group1 == pair[2] & resPairwise$group2 == pair[1]))
          if (length(idx) > 0) {
            as.character(resPairwise$stars[idx])
          } else {
            "ns"
          }
        })
        
        y_max <- max(stats$Mean + stats$sd, na.rm = TRUE)
        y_pos <- seq(y_max + 0.1, y_max + 0.1 + 0.1 * length(post_hoc_pairs), length.out = length(post_hoc_pairs))
        
        if (length(post_hoc_pairs) > 0 && length(annotations) > 0) {
          resplot <- resplot + 
            ggsignif::geom_signif(
              comparisons = post_hoc_pairs,
              y_position = y_pos,
              annotations = annotations,
              tip_length = 0.01,
              textsize = 8 / 1.5,
              vjust = 0.8
            )
        } else {
          print("No valid post-hoc comparisons or annotations to plot.")
        }
      } else {
        print("No pairwise comparisons to plot.")
      }
      
      output$analysis_output <- renderUI({
        steps_formatted <- gsub("Step", "<br><br>Step", steps)
        steps_formatted <- gsub("Group", "<br>&nbsp;&nbsp;&nbsp;&nbsp;Group", steps_formatted)
        HTML(steps_formatted)
      })
      
      
      return(list(Table = stats, TableTest = resTTest, Plot = resplot))
    } else {
      return(list(Table = NULL, TableTest = NULL, Plot = NULL))
    }
  })
  
  
  # Render statistical analysis results
  output$TabStat <- renderDT({
    StatisticalAnalysisResults()$Table
  },
  options = list(
    searching = FALSE,
    dom = 't' # Only display the table
  )
  )
  
  output$PlotStat <- renderPlot({
    StatisticalAnalysisResults()$Plot
  })
  
  output$TabTTest <- renderDT({
    StatisticalAnalysisResults()$TableTest
  },
  options = list(
    searching = FALSE,
    dom = 't' # Only display the table
  ))
  
  ### End Statistic ####
  
  ### Start DATAVERSE ####
  
  observeEvent(input$APIkey,{
    pathORCA <- system.file("Data", package = "ORCA")
    
    if(input$APIkey != "") # the last key used is saved
      write(input$APIkey,file = paste0(pathORCA,"/.APIkey"))
    
  })
  
  DataAnalysisModule0 <- list(wbResult = wbResult0,
                              wbquantResult = wbquantResult0,
                              endocResult = endocResult0,
                              elisaResult = elisaResult0,
                              pcrResult = pcrResult0,
                              cytotoxResult = cytotoxResult0,
                              facsResult = facsResult0,
                              bcaResult = bcaResult0,
                              ifResult = ifResult0,
                              facsResult = facsResult0)
  
  observe({
    listDataAnalysisModule = reactiveValuesToList(DataAnalysisModule)
    namesAnalysis = sapply(names(listDataAnalysisModule), function(x){
      # if(x == "wbquantResult"){
      #   if(! identical(DataAnalysisModule[[x]]$NormWBanalysis,DataAnalysisModule0[[x]]$NormWBanalysis) )
      #     return(x)
      #   else
      #     return("")
      # }
      # else 
      if(!identical(DataAnalysisModule[[x]],DataAnalysisModule0[[x]]) )
        return(x)
      else
        return("")
    })
    namesAnalysis = namesAnalysis[namesAnalysis!= ""]
    if(length(namesAnalysis)>0)
      updateSelectizeInput(inputId = "selectAnalysis_DV",
                           choices = MapAnalysisNames[namesAnalysis])
  })
  
  observeEvent(input$DataverseUpload_Button,{
    
    if(input$selectAnalysis_DV != ""){
      
      if(input$Title_DV == "" && input$Author_DV == "" && input$Description_DV == "" &&
         input$AuthorAff_DV== "" && input$ContactN_DV == "" && input$ContactEmail_DV == "")
        output$LoadingError_DATAVERSE = renderText("Error: missing information")
      else{
        # creation of a temporary folder
        tempfolder = paste0(tempdir(check = T),"/ORCA")
        system(paste0("mkdir ",tempfolder))
        system(paste0("mkdir ",tempfolder,"/dataset"))
        # create the metadata
        result <- rjson::fromJSON(file = system.file("Docker/Dataverse","metadata.json", package = "ORCA") )
        result$dataset_title = input$Title_DV
        result$dataset_description = paste0(input$Description_DV,"\n This dataset has been obtained using the ORCA application, specifically the module: ", input$selectAnalysis_DV)
        result$author_name = input$Author_DV
        result$author_affiliation= input$AuthorAff_DV
        result$dataset_contact_name = input$ContactN_DV
        
        # check of it is a correct email
        result$dataset_contact_email = input$ContactEmail_DV
        
        # result$subject = as.vector(result$subject)
        write(rjson::toJSON(result), file=paste0(tempfolder,"/metadata.json") )
        
        # move the file in the temporary folder
        filenameTMP = paste0(tempfolder,"/dataset/",
                             gsub(pattern = " ", 
                                  x = input$selectAnalysis_DV,
                                  replacement = ""),".xlsx")
        
        saveExcel(filename = filenameTMP,
                  ResultList= DataAnalysisModule[[ names(MapAnalysisNames[MapAnalysisNames == input$selectAnalysis_DV])]],
                  analysis = input$selectAnalysis_DV, PanelStructures)
        
        saveRDS(DataAnalysisModule[[ names(MapAnalysisNames[MapAnalysisNames == input$selectAnalysis_DV])]] ,
                file = paste0(tempfolder,"/dataset/ORCA_",
                              gsub(pattern = " ", 
                                   x = input$selectAnalysis_DV,
                                   replacement = ""),".RDs"))
        # docker run
        ORCA::docker.run(params = paste0("--volume ", tempfolder,
                                         ":/Results/ -d qbioturin/dataversemanagement createDataset --dataverse root --json /Results/metadata.json --datafile_dir /Results/dataset") 
        )
        
      }
      
    }
  })
  
  #initiate_sword_dataset()
  #add_dataset_file()
  #publish_dataset()
  
  #### End DATAVERSE ####
  
  ### LOAD analysis ####
  UploadDataAnalysisModuleAllFalse  = reactiveValues(FlagALL = F,
                                                     FlagUpdate = F,
                                                     FlagWB = F,
                                                     FlagWBquant = F,
                                                     FlagPCR = F,
                                                     FlagELISA = F,
                                                     FlagCYTOTOX = F,
                                                     FlagENDOC = F,
                                                     FlagBCA = F,
                                                     FlagFACS = F)
  UploadDataAnalysisModule = reactiveValues(FlagALL = F,
                                            FlagUpdate = F,
                                            FlagWB = F,
                                            FlagWBquant = F,
                                            FlagPCR = F,
                                            FlagELISA = F,
                                            FlagCYTOTOX = F,
                                            FlagENDOC = F,
                                            FlagBCA = F,
                                            FlagFACS = F)
  
  # general upload in the app
  observeEvent(input$loadAnalysis_Button,{
    
    manageSpinner(TRUE)
    if( is.null(input$loadAnalysis_file) || !all(file.exists(input$loadAnalysis_file$datapath) ) ) {
      showAlert("Error", "Please select one RDs file generated throught the Data Analysis module.", "error", 5000)
      return()
    }
    mess = readRDS(input$loadAnalysis_file$datapath)
    
    messNames = names(mess)
    if("Flags"%in% messNames) messNames = messNames[ messNames != "Flags"]
    
    if(!(all(messNames %in% names(DataAnalysisModule)) ||
         all(messNames %in% names(elisaResult)) ||
         all(messNames %in% names(wbResult)) ||
         all(messNames %in% names(wbquantResult)) ||
         all(messNames %in% names(pcrResult)) ||
         all(messNames %in% names(cytotoxResult)) ||
         all(messNames %in% names(endocResult)) ||
         all(messNames %in% names(bcaResult)) ||
         all(messNames %in% names(facsResult)) ) ) {
      showAlert("Error", "The file must be RDs saved throught the Data Analysis module.", "error", 5000)
      return()
    }
    
    if(all(messNames %in% names(DataAnalysisModule)) ){
      DataAnalysisModule <- mess
      UploadDataAnalysisModule$FlagALL = T
    }else if( all(messNames %in% names(bcaResult)) ){
      DataAnalysisModule$bcaResult <- mess
      UploadDataAnalysisModule$FlagBCA = T
    }else if( all(messNames %in% names(wbResult)) ){
      DataAnalysisModule$wbResult <- mess
      UploadDataAnalysisModule$FlagWB = T
    }else if( all(messNames %in% names(wbquantResult)) ){
      DataAnalysisModule$wbquantResult <- mess
      UploadDataAnalysisModule$FlagWBquant = T
    }else if( all(messNames %in% names(pcrResult)) ){
      DataAnalysisModule$pcrResult <- mess
      UploadDataAnalysisModule$FlagPCR = T
    }else if(all(messNames %in% names(endocResult)) ){
      DataAnalysisModule$endocResult <- mess
      UploadDataAnalysisModule$FlagENDOC = T
    }else if(all(messNames %in% names(elisaResult)) ){
      DataAnalysisModule$elisaResult <- mess
      UploadDataAnalysisModule$FlagELISA = T
    }else if(all(messNames %in% names(cytotoxResult)) ){
      DataAnalysisModule$cytotoxResult <- mess
      UploadDataAnalysisModule$FlagCYTOTOX = T
    }else if( all(messNames %in% names(ifResult)) ){
      DataAnalysisModule$ifResult <- mess
      UploadDataAnalysisModule$FlagIF = T
    }else if( all(messNames %in% names(facsResult)) ){
      DataAnalysisModule$facsResult <- mess
      UploadDataAnalysisModule$FlagFACS = T
    }
    
    UploadDataAnalysisModule$FlagUpdate = T
    manageSpinner(FALSE)
    showAlert("Success", "The RDs file has been uploaded  with success." , "success", 5000)
    
  })
  
  observeEvent(UploadDataAnalysisModule$FlagUpdate,{
    if(UploadDataAnalysisModule$FlagUpdate){
      
      if(UploadDataAnalysisModule$FlagWB || UploadDataAnalysisModule$FlagALL){
        UploadRDs(Flag = "WB",
                  session = session,
                  output = output,
                  DataAnalysisModule = DataAnalysisModule,
                  Result = wbResult, 
                  FlagsExp = Flags,
                  PanelStructures = PanelStructures)
      } else if(UploadDataAnalysisModule$FlagWBquant || UploadDataAnalysisModule$FlagALL){
        UploadRDs(Flag = "WBquant",
                  session = session,
                  output = output,
                  DataAnalysisModule = DataAnalysisModule,
                  Result = wbquantResult, 
                  FlagsExp = FlagsWBquant,
                  PanelStructures = PanelStructures)
      } else if(UploadDataAnalysisModule$FlagBCA || UploadDataAnalysisModule$FlagALL){
        UploadRDs(Flag = "BCA",
                  session = session,
                  output = output,
                  DataAnalysisModule = DataAnalysisModule,
                  Result = bcaResult, 
                  FlagsExp = FlagsBCA)
      }
      else if(UploadDataAnalysisModule$FlagPCR || UploadDataAnalysisModule$FlagALL){
        UploadRDs(Flag = "PCR",
                  session = session,
                  output = output,
                  DataAnalysisModule = DataAnalysisModule,
                  Result = pcrResult, 
                  FlagsExp = FlagsPCR)
      }
      else if(UploadDataAnalysisModule$FlagENDOC || UploadDataAnalysisModule$FlagALL){
        UploadRDs(Flag = "ENDOC",
                  session = session,
                  output = output,
                  DataAnalysisModule = DataAnalysisModule,
                  Result = endocResult, 
                  FlagsExp = FlagsENDOC)
      }
      else if(UploadDataAnalysisModule$FlagELISA || UploadDataAnalysisModule$FlagALL){
        UploadRDs(Flag = "ELISA",
                  session = session,
                  output = output,
                  DataAnalysisModule = DataAnalysisModule,
                  Result = elisaResult, 
                  FlagsExp = FlagsELISA)
        
      }
      else if(UploadDataAnalysisModule$FlagCYTOTOX || UploadDataAnalysisModule$FlagALL){
        UploadRDs(Flag = "CYTOTX",
                  session = session,
                  output = output,
                  DataAnalysisModule = DataAnalysisModule,
                  Result = cytotoxResult, 
                  FlagsExp = FlagsCYTOTOX)
        
      }
      else if(UploadDataAnalysisModule$FlagIF || UploadDataAnalysisModule$FlagALL){
        UploadRDs(Flag = "IF",
                  session = session,
                  output = output,
                  DataAnalysisModule = DataAnalysisModule,
                  Result = ifResult, 
                  FlagsExp = FlagsIF)
        
      }
      else if(UploadDataAnalysisModule$FlagFACS || UploadDataAnalysisModule$FlagALL){
        UploadRDs(Flag = "FACS",
                  session = session,
                  output = output,
                  DataAnalysisModule = DataAnalysisModule,
                  Result = facsResult, 
                  FlagsExp = FlagsFACS)
        
      }
      UploadDataAnalysisModule = UploadDataAnalysisModuleAllFalse
    }
    
  })
  
  #### 
  
  #----------------------------------------------------------------------------------
  # OTHER ANALYSIS TO DO
  
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
      reportpath = system.file("Shiny","report.Rmd", package = "ORCA")
      if(reportpath == ""){
        # In this case it means that we are inside the docker
        reportpath = "~/../home/Shiny/report.Rmd"
      }
      rmarkdown::render(reportpath,
                        output_file = file,
                        output_format = "html_document",
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
  
  # Download statistical analysis button
  output$downloadStatisticalAnalysis <- downloadHandler(
    filename = function() {
      paste("Statistical_Analysis_", input$StatAnalysis, ".RDs", sep = "")
    },
    content = function(file) {
      results <- StatisticalAnalysisResults()
      if (!is.null(results$Table)) {
        manageSpinner(TRUE)
        saveRDS(results, file = file)
        manageSpinner(FALSE)
        showAlert("Success", "Download completed successfully", "success", 2000)
      }
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
        !is.null(wbquantResult$AdjRelDensity))
      return (FALSE)
    if (!is.null(pcrResult$Initdata) || !is.null(pcrResult$selectPCRcolumns) ||
        !is.null(pcrResult$data) || !is.null(pcrResult$PCRnorm) || !is.null(pcrResult$NewPCR) ||
        !is.null(pcrResult$plotPCR))
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
    if (!is.null(endocResult$Initdata) || !is.null(endocResult$data) ||
        !is.null(endocResult$TablePlot) || !is.null(endocResult$dataFinal) ||
        !is.null(endocResult$ENDOCcell_TIME) || !is.null(endocResult$ENDOCcell_SN) ||
        !is.null(endocResult$MapBaseline) || !is.null(endocResult$MapBlank))
      return (FALSE)
    if (!is.null(ifResult$Initdata) || !is.null(ifResult$data) ||
        !is.null(ifResult$dataFinal))
      return (FALSE)
    if (!is.null(facsResult$Initdata) || !is.null(facsResult$data) ||
        !is.null(facsResult$dataFinal))
      return (FALSE)
    
    return (TRUE)
  }
  
  # END DOWNLOAD
}