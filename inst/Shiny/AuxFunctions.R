# function to invoke to send shiny error or success messages
showAlert <- function(title, text, type = "info", time) {
  shinyalert(title = title, text = text, type = type, timer = time)
}

manageSpinner <- function(isDownloading) {
  if(isDownloading == TRUE) {
    show_modal_spinner()
  } else {
    remove_modal_spinner()
  }
}

areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}
resetPanel <- function(type, flags = NULL, panelStructures = NULL, numberOfPlanes = NULL, planeSelected = NULL, result, output = NULL, panelData = NULL) {
  switch(type,
         "WB" = {
           flags$ShowTif <- FALSE
           flags$LanesCut <- FALSE
           flags$CutTab <- "V"
           flags$IDlane <- 0
           
           panelStructures$data <- panelData  
           
           numberOfPlanes$N <- 0
           planeSelected$First <- NULL
           
           result$Normalizer <- NULL
           result$Im <- NULL
           result$Planes <- NULL
           result$TruncatedPanelsValue <- NULL
           result$PanelsValue <- NULL
           result$Plots <- NULL
           result$TruncatedPlots <- NULL
           result$pl <- NULL
           result$AUCdf <- data.frame(SampleName = "-", Truncation = "-", AUC = "-")
           
           output$DataPlot <- renderPlot({}) 
           output$AUC <- renderDT({data.frame()})  
           output$AUC_RelDens <- renderDT({data.frame()})
           output$AUC_AdjRelDens <- renderDT({data.frame()})
           output$plot_AdjRelDens <- renderPlot({})
         },
        "PCR" = {
           flags$norm <- FALSE
           flags$baseline <- FALSE
           
           result$Initdata <- NULL
           result$selectPCRcolumns <- NULL
           result$data <- NULL
           result$PCRnorm <- NULL
           result$BaselineExp <- NULL
           result$plotPRC <- NULL
           result$NewPCR <- NULL
        },
        "ENDOC" = {
           result$Initdata <- NULL
           result$data <- NULL
           result$TablePlot <- NULL
           result$dataFinal <- NULL
           result$ENDOCcell_TIME <- NULL
           result$ENDOCcell_COLOR <- NULL
           result$MapBaseline <- NULL
           result$MapBlank <- NULL

           flags$cellCoo <- NULL
           flags$AllExp <- NULL
           flags$BASEselected <- NULL
           flags$BLANCHEselected <- NULL
           flags$EXPselected <- NULL
           flags$EXPcol <- NULL
        },
        "ELISA" = {
          FlagsELISA <- reactiveValues(cellCoo = NULL,
                                       AllExp = "",
                                       BASEselected = "",
                                       STDCselected = "",
                                       BLANCHEselected = "",
                                       EXPselected = "",
                                       EXPcol = NULL)
          
           result$Initdata <- NULL
           result$data <- NULL
           result$TablePlot <- NULL
           result$dataFinal <- NULL
           result$ELISAcell_EXP <- NULL
           result$ELISAcell_SN <- NULL
           result$ELISAcell_COLOR <- NULL
           result$MapBaseline <- NULL
           result$MapBlank <- NULL
           result$Tablestandcurve <- NULL
           result$Regression <- NULL
           
           flags$cellCoo <- NULL
           flags$AllExp <- NULL
           flags$BASESelected <- ""
           flags$STDCselected <- ""
           flags$BLANCHEselected
           flags$EXPselected <- ""
           flags$EXPcol <- NULL
        },
        "CYTOTOX" = {
           result$Initdata <- NULL
           result$data <- NULL
           result$TablePlot <- NULL
           result$dataFinal <- NULL
           result$CYTOTOXcell_EXP <- NULL
           result$CYTOTOXcell_REP <- NULL
           result$CYTOTOXcell_SN <- NULL
           result <- MapBaseline <- NULL
           
           flags$cellCoo <- NULL
           flags$AllExp <- NULL
           flags$BASEselected <- NULL
           flags$EXPselected <- NULL
           flags$EXPcol <- NULL
        },
        "FACS" = {
           result$Initdata <- NULL
           result$data <- NULL
           result$dataFinal <- NULL
           result$depth <- NULL
           result$depthCount <- NULL
           result$name <- NULL
           result$statistics <- NULL
           result$cells <- NULL
          
           flags$actualLevel <- NULL
           flags$allLevel <- NULL
           flags$actualPath <- NULL
        },
        error = function(cond) {
           showAlert("Error", "an error occured", "error", 5000)
        }
  )
}


# function called when you need to read a file
readfile <- function(filename, type, isFileUploaded, colname = TRUE, namesAll = namesAll, allDouble = FALSE, colors = FALSE) {
  out <- tryCatch({
    switch(type, "tif" = {
      if(!isFileUploaded || !file.exists(filename)) {
        return (list(message = "Please, select a TIF file",call = ""))
      } else if(tolower(tools::file_ext(filename)) != "tif" && tolower(tools::file_ext(filename)) != "tiff") {
        return(list(message = "Please, upload a file with a .tif extension.", call = ""))
      } 
      else {loadImage(filename)}
    },
    "RDs" = {
      if(is.null(filename) || !file.exists(filename)) {
        return(list(message = "Please, select a RDS File!", call = ""))
      }  else if(tolower(tools::file_ext(filename)) != "rds") {
        return(list(message = "Please, upload a file with a .rds extension.", call = ""))
      } 
      
      x = readRDS(filename)
      
      if(!all(names(x) %in% namesAll)) {
        return(list(message = "The RDs file must be generated from Data Analysis module.", call = ""))
      }
      
      if(is.null(x$AUCdf)) {
        return(list(message = "The WB analysis must contain the AUC table", call = ""))
      }
      
      x
    },
    "RDsMulti" = {
      result <- list(data = list(), error = NULL)
      filenames <- filename 
      
      if (is.null(filenames) || length(filenames) == 0) {
        result$error <- "Please select one or more .rds files."
        return(result)
      }
      
      for (filename in filenames) {
        if (!file.exists(filename)) {
          result$error <- paste("The file", filename, "does not exist.")
          return(result)
        }
        
        if (tolower(tools::file_ext(filename)) != "rds") {
          result$error <- paste("The file", filename, "is not a .rds file")
          return(result)
        }
        
      }
      return(result)
    },
    "Excel" = {
      if(!isFileUploaded || !file.exists(filename)) {
        return (list(message = "Please, select an Excel file", call = ""))
      } else if(tolower(tools::file_ext(filename)) != "xls" && tolower(tools::file_ext(filename)) != "xlsx") {
        return(list(message = "Please, upload a file with a .xls or .xlsx extension.", call = ""))
      } 
      
      x = readxl::read_excel(filename, col_names = colname)
      x <- x[, !is.na(colnames(x))]
      
      if (allDouble) {
        xstr = which(sapply(x, function(col) !is.numeric(col)))
        if (length(xstr) > 0) {
          for (i in which(sapply(x, is.numeric))) {
            x[[i]] = as.double(x[[i]])
          }
        }
      }
      if (colors) {
        wb = loadWorkbook(filename)
        sheetName = wb$sheet_names[1]
        
        l = lapply(wb$styleObjects, function(x) {
          if (x$sheet == sheetName) {
            if (all(areColors(paste0("#", unname(x$style$fill$fillFg))))) {
              color = paste0("#", unname(x$style$fill$fillFg))
              if (grep(color, pattern = "^#FF") && all(areColors(gsub(replacement = "#", x = color, pattern = "^#FF")))) {
                color = gsub(replacement = "#", x = color, pattern = "^#FF")
              }
            } else {
              color = randomcoloR::randomColor(1)
            }
            if (color == "#FFFFFF") { 
              color = "white"
            }
            df = data.frame(row = x$rows, col = x$cols,
                            fill = ifelse(length(x$style$fill$fillFg) > 0 && color != "white", 
                                          color,
                                          "white"))
          }
        })
        l = do.call(rbind, l[lengths(l) > 0])
        SN = table(l$fill)
        l$SN = paste0("Color ", match(l$fill, names(SN)))
        
        tb.SN = matrix("", nrow = max(l$row), ncol = max(l$col))
        
        for (j in 1:ncol(tb.SN)) {
          fill_col = dplyr::filter(l, col == j)
          if (!all(fill_col$fill == "white")) {
            tb.SN[fill_col$row, j] = fill_col$SN
          }
        }
        
        tb.SN <- tb.SN[, colSums(tb.SN != "") > 0]
        
        col = dplyr::select(l, fill, SN) %>% dplyr::distinct()
        vectcol = col$fill
        names(vectcol) = col$SN
        
        return(list(x = x, SNtable = tb.SN, fill = vectcol))
      }
      
      return(x)
    },
    
    {list(message = "Unsupported file type.", call = "")} #default switch case
    )
  }, 
  error = function(cond) {
    list(message = "An error occurred.", call = "")
  })
  return(out)
}

# load the image
LoadImage = function(pathImage){
  im <- OpenImageR::readImage(pathImage,as.is = T, convert=TRUE)
  
  if( length(dim(im)) != 2 && dim(im)[3] == 4 ){ # if the tiff image is already gray
    im = im[,,-4]
    im = rgb_2gray(im)
  }else if( length(dim(im)) != 2 && dim(im)[3] == 3 ){
    im = rgb_2gray(im)
  }
  
  
  JpegImage = tempfile(fileext = "jpeg")
  
  # unlink()
  jpeg(file=JpegImage,
       width = dim(im)[2],
       height = dim(im)[1],
       units = "px")
  imageShow(im)
  dev.off()
  
  img <- jpeg::readJPEG(JpegImage, native = F)
  
  return(list(WB = im, RGB = img))
}

AUCfunction<-function(AUCdf,PanelsValue,bind=T,session = session,SName="1",AUCdf.new=NULL){
  if(is.null(AUCdf.new)){
    
    if(length(AUCdf[,1])==1 & AUCdf$AUC[1] == "-")
    {
      AUCdf.new <- AUCdf
    }else{
      AUCdf2 = AUCdf %>% filter(SampleName == SName)
      
      if(length(AUCdf2[,1]) == 0){
        AUCdf.new = AUCdf2
        AUCdf.new[1,] = rep(NA,length(names(AUCdf2)))
        AUCdf.new$Truncation = "-"
        AUCdf.new$SampleName = SName
      }else{
        AUCdf.new <- AUCdf2[length(AUCdf2$Truncation),]
      }
    }
  }
  
  PanelsValue.Lane <- PanelsValue[which(PanelsValue$ID == SName),]
  id <- order(PanelsValue.Lane$Y)
  AUCdf.new$AUC <- round(
    sum(diff(PanelsValue.Lane$Y[id])*rollmean(PanelsValue.Lane$Values[id],2)),
    digits = 4)
  AUCdf.new$SampleName <- paste(SName)
  
  if(length(AUCdf[,1])==1 & AUCdf$AUC[1] == "-")
  {
    A<-AUCdf.new
  }else{
    A<-rbind(AUCdf,AUCdf.new) 
  }
  return(unique(A)) 
}

saveExcel <- function(filename, ResultList, analysis, PanelStructures = NULL) {
  if (file.exists(filename)) {
    file.remove(filename)
  }
  
  switch(analysis, 
      "WB" = {
        wb <- createWorkbook("WB")
        
        addWorksheet(wb, "WBimage")
        ListIm <- ResultList[["Im"]] 
        im <- ListIm$RGB
        
        plot(c(1, dim(im)[2]), c(1, dim(im)[1]), type='n', ann=FALSE)
        rasterImage(im, 1, 1, dim(im)[2], dim(im)[1])
        insertPlot(wb = wb, sheet="WBimage")
        
        addWorksheet(wb, "WBimage and protein bands")
        plot(c(1, dim(im)[2]), c(1, dim(im)[1]), type='n', ann=FALSE)
        rasterImage(im, 1, 1, dim(im)[2], dim(im)[1])
        if (!is.null(PanelStructures$data)) {
          for (i in seq_len(nrow(PanelStructures$data))) {
            with(PanelStructures$data, {
              rect(xmin[i], ymin[i], xmax[i], ymax[i], border = "red")
            })
          }
        }
        insertPlot(wb = wb, sheet = "WBimage and protein bands")
        
        startRow <- 22
        writeDataTable(wb, PanelStructures$data, sheet = "WBimage and protein bands", startRow = startRow, startCol = 1)
        
        addWorksheet(wb, "Plot")
        print(ResultList[["Plots"]])
        insertPlot(wb = wb, sheet="Plot")
        
        addWorksheet(wb, "Truncated Plot")
        print(ResultList[["TruncatedPlots"]])
        insertPlot(wb = wb, sheet="Truncated Plot")
        
        addWorksheet(wb, "AUC")
        finaldata = ResultList[["AUCdf"]]
        writeDataTable(wb, finaldata, sheet="AUC")
        
        saveWorkbook(wb, filename)
        return(1)
      },
      "WB comparison" = {
      ## Create a new workbook
      wb <- createWorkbook("WB comparison")
      
      ## initial data
      addWorksheet(wb,"Normalizer WB")
      writeDataTable(wb, sheet = "Normalizer WB", ResultList[["NormWBanalysis_filtered"]])
      
      addWorksheet(wb,"WB")
      writeDataTable(wb, sheet = "WB", ResultList[["WBanalysis_filtered"]])
      
      ### Analysis
      
      addWorksheet(wb,"RelDensitiy")
      writeDataTable(wb, sheet = "RelDensitiy", ResultList[["RelDensitiy"]])
      
      addWorksheet(wb,"AdjRelDensitiy")
      writeDataTable(wb, sheet = "AdjRelDensitiy", ResultList[["AdjRelDensitiy"]])
      
      addWorksheet(wb,"Barplot AdjRelDensitiy")
      if(!is.null( ResultList[["AdjRelDensitiy"]])){
        print(
          ResultList[["AdjRelDensitiy"]] %>% 
            mutate(Normalizer = paste0("Sample: ",SampleName ),
                   WB = paste0("Sample: ",SampleName))  %>%
            ggplot() +
            geom_bar(aes(x = SampleName,
                         y = AdjRelDens,
                         fill = Normalizer ),
                     stat = "identity" ) +
            theme_bw()
        )
      }
      insertPlot(wb, sheet = "Barplot AdjRelDensitiy")
    }, 
    "RT-qPCR" = {
      wb <- createWorkbook("RTqPCR")
      
      addWorksheet(wb,"Table")
      writeDataTable(wb, sheet = "Table", ResultList[["Initdata"]])
      
      addWorksheet(wb,"Norm PRC")
      writeDataTable(wb,ResultList[["NewPCR"]], sheet="Norm PRC")
      
      insertPlot(wb = wb,  sheet="Norm PRC",
                 startCol=dim(ResultList[["NewPCR"]])[2]+ 2)
      
    },
    "ENDOC" = {
      wb <- createWorkbook("ENDOC")  
      
      if (!is.null(ResultList$Initdata) && is.data.frame(ResultList$Initdata)) {
        addWorksheet(wb, "TablePlot")  # Aggiunge un nuovo foglio al workbook
        writeDataTable(wb, ResultList$Initdata, sheet = "TablePlot")  # Scrive i dati nel foglio
        print("Initdata scritto nel foglio Excel")
      } else {
        print("Errore: Initdata non disponibile o non è un data.frame")
      }
      
      if (!is.null(ResultList$dataFinal) && is.data.frame(ResultList$dataFinal)) {
        addWorksheet(wb, "Results Analysis")  
        writeDataTable(wb, ResultList$dataFinal, sheet = "Results Analysis")  
      }
    },
    "FACS" = {
      wb <- createWorkbook("FACS")
      
      print(ResultList$Initdata)
      print(ResultList$data)
      print(ResultList$dataFinal)
      
      if (!is.null(ResultList$Initdata) && is.data.frame(ResultList$Initdata)) {
        addWorksheet(wb, "Init data") 
        writeDataTable(wb, ResultList$Initdata, sheet = "Init data")  
      }
      if (!is.null(ResultList$dataFinal) && is.data.frame(ResultList$dataFinal)) {
        addWorksheet(wb, "Final data") 
        writeDataTable(wb, ResultList$dataFinal, sheet = "Final data")  
      }
    }
  )
  
  saveWorkbook(wb, filename)  # Salva il workbook
  return(1)  # Restituisce 1 per indicare il successo
}

tableExcelColored = function(session, output, Result, FlagsExp, type){
  switch(type,
      "Initialize" = {
      ExpDataTable = Result$Initdata

      if(is.null(FlagsExp$EXPcol)){
        ExpDataTable.colors = matrix("",nrow = nrow(ExpDataTable),ncol=ncol(ExpDataTable))
      }else{
        ExpDataTable.colors = Result[[grep(x=names(Result), pattern = "cell_COLOR", value = T)]]
        
        if(is.null(ExpDataTable.colors) ) {
          ExpDataTable.colors = matrix("", nrow = nrow(ExpDataTable), ncol = ncol(ExpDataTable))
        }else if( nrow(ExpDataTable.colors) != nrow(ExpDataTable) ||  ncol(ExpDataTable.colors) != ncol(ExpDataTable) ){
          ExpDataTable.colors = ExpDataTable.colors[1:nrow(ExpDataTable),1:ncol(ExpDataTable)]
        }
      }
      
      completeExpDataTable = cbind(ExpDataTable,ExpDataTable.colors)
      
      cols.keep <- paste0('V',1:length(ExpDataTable[1,])) 
      cols.color <- paste0('Col',1:length(ExpDataTable[1,]))
      colnames(completeExpDataTable) = c(cols.keep,cols.color)
      
      if(is.null(FlagsExp$EXPcol)){
        EXPcol = ""
        names(EXPcol) = "white"
      }else{
        EXPcol = FlagsExp$EXPcol
      }
      ExpDataTable = datatable(completeExpDataTable,
                               filter = 'none',
                               #server = FALSE,
                               selection = list(mode = 'single', target = 'cell'),
                               rownames= FALSE,
                               options = list(
                                 initComplete = JS(
                                   "function(settings, json) {",
                                   "$(this.api().table().body()).find('td').each(function() {",
                                   "  var bgColor = $(this).css('background-color');",
                                   "  if (bgColor === 'rgb(255, 255, 255)' || bgColor === 'white') {", 
                                   "    $(this).addClass('non-clickable').css({'pointer-events': 'none'});",
                                   "  }",
                                   "});",
                                   "}"),
                                 dom = 't',
                                 pageLength = -1,
                                 info = FALSE,
                                 #scrollX = TRUE,
                                 #lengthChange = FALSE,
                                 columnDefs = list(list(targets = cols.color, visible = FALSE))
                               )) %>%
        formatStyle(cols.keep,
                    cols.color,
                    backgroundColor = styleEqual(names(EXPcol), EXPcol) )
      cell_COLOR <- ExpDataTable.colors
      cell_TIME <- cell_EXP <- cell_REP <- matrix(
        "",
        nrow = length(ExpDataTable$x$data[,1]),
        ncol = length(ExpDataTable$x$data[1,])
      )
      
      if(length(grep(x=names(Result),pattern = "cell_COLOR", value = T))>0)
        Result[[grep(x=names(Result),pattern = "cell_COLOR", value = T)]] <- cell_COLOR
      if(length(grep(x=names(Result),pattern = "cell_EXP", value = T))>0)
        Result[[grep(x=names(Result),pattern = "cell_EXP", value = T)]]<- cell_EXP
      if(length(grep(x=names(Result),pattern = "cell_REP", value = T))>0)
        Result[[grep(x=names(Result),pattern = "cell_REP", value = T)]]<- cell_REP
      if(length(grep(x=names(Result),pattern = "cell_TIME", value = T))>0)
        Result[[grep(x=names(Result),pattern = "cell_TIME", value = T)]]<- cell_TIME
      
      Result$TablePlot = ExpDataTable
  }, "Update" =  {
      ColorsSN = rainbow(n = 50, alpha = 0.5)[sample(50, size = 50, replace = FALSE)]
      
      if(is.null(FlagsExp$EXPcol)) {
        print("No existing color mapping found. Creating new one.")
        EXPcol = setNames(ColorsSN[1:length(FlagsExp$AllExp)], FlagsExp$AllExp)
        EXPcol[names(EXPcol) == ""] <- "white"  
        FlagsExp$EXPcol <- EXPcol
      } else {
        print("Existing color mapping found. Updating if necessary.")
        SNnew = FlagsExp$AllExp[!FlagsExp$AllExp %in% names(FlagsExp$EXPcol)]
        if(length(SNnew) > 0) {
          print(paste("New SNs found:", paste(SNnew, collapse=", ")))
          colNew = ColorsSN[!ColorsSN %in% FlagsExp$EXPcol][1:length(SNnew)]
          names(colNew) = SNnew
          EXPcol = c(FlagsExp$EXPcol, colNew)
          EXPcol[names(EXPcol) == ""] <- "white"
          FlagsExp$EXPcol <- EXPcol
        } else {
          print("No new SNs to update.")
        }
      }
      
      ExpDataTable = Result$TablePlot$x$data
      completeExpDataTable = cbind(Result$Initdata, Result[[grep(x=names(Result), pattern = "cell_COLOR", value = TRUE)]])
      colnames(completeExpDataTable) = colnames(ExpDataTable)
      
      cols.color = grep(x = colnames(ExpDataTable), pattern = "Col", value = TRUE)
      cols.keep = grep(x = colnames(ExpDataTable), pattern = "V", value = TRUE)
      
      Result$TablePlot = datatable(completeExpDataTable,
                               filter = 'none',
                               selection = list(mode = 'single', target = 'cell'),
                               rownames= FALSE,
                               options = list(
                                 initComplete = JS(
                                   "function(settings, json) {",
                                   "$(this.api().table().body()).find('td').each(function() {",
                                   "  var bgColor = $(this).css('background-color');",
                                   "  if (bgColor === 'rgb(255, 255, 255)' || bgColor === 'white') {", 
                                   "    $(this).addClass('non-clickable').css({'pointer-events': 'none'});",
                                   "  }",
                                   "});",
                                   "}"),
                                 dom = 't',
                                 pageLength = -1,
                                 info = FALSE,
                                 columnDefs = list(list(targets = cols.color, visible = FALSE))
                               )) %>%
        formatStyle(cols.keep,
                    cols.color,
                    backgroundColor = styleEqual(names(FlagsExp$EXPcol), FlagsExp$EXPcol))
      
      print("Table and colors updated.")
    }
  )
}

get_formatted_data <- function(colors, color_names, result, singleValue, analysis) {
  if (length(colors) == 0) {
    return(data.frame(Color = character(), Values = character(), ExperimentalCondition = character(), ColorCode = character()))
  }
  
  formatted_data <- vector("list", length(colors))
  column_COLOR <- paste0(analysis, "cell_COLOR")  
  column_TIME <- paste0(analysis, "cell_TIME")  
  
  # Set variable names based on analysis type
  if (analysis == "ELISA" || analysis == "CYTOTOX") {
    value1 = "Sample Name"
    value2 = "Experimental condition"
    column_EXP <- paste0(analysis, "cell_SN") 
  } else if (analysis == "ENDOC") {
    value1 = "Experimental condition"
    value2 = "Time"
    column_EXP <- paste0(analysis, "cell_EXP") 
  }
  
  for (i in seq_along(colors)) {
    matching_indices <- which(result[[column_COLOR]] == color_names[i], arr.ind = TRUE)
    if (nrow(matching_indices) > 0) {
      selected_values <- apply(matching_indices, 1, function(idx) {
        result$Initdata[idx["row"], idx["col"]]
      })
      
      formatted_output <- paste(unlist(selected_values), collapse = " - ")
      
      time_values <- apply(matching_indices, 1, function(idx) {
        if (analysis == "ELISA") {
          val <- result$ELISAcell_EXP[idx["row"], idx["col"]]
        } else if (analysis == "CYTOTOX") {
          val <- result$CYTOTOXcell_EXP[idx["row"], idx["col"]]
        } else {
          val <- result$ENDOCcell_TIME[idx["row"], idx["col"]]
        }
        if (!is.na(val) && !is.null(val) && val != "") val else ""
      })
      
      if (analysis == "CYTOTOX") {
        rep_values <- apply(matching_indices, 1, function(idx) {
          val <- result$CYTOTOXcell_REP[idx["row"], idx["col"]]
          if (!is.na(val) && !is.null(val) && val != "") val else ""
        })
      }
      
      
      time_output <- paste(unlist(time_values), collapse = " - ")
      if (analysis == "CYTOTOX") rep_output <- paste(unlist(rep_values), collapse = " - ")
      
      exp_values <- apply(matching_indices, 1, function(idx) {
        result[[column_EXP]][idx["row"], idx["col"]]
      })
      
      if (length(unique(exp_values)) == 1) {
        exp_condition <- ifelse(exp_values[1] == "" || is.na(exp_values[1]), "-", exp_values[1])
      } else {
        exp_condition <- "No matching between values"
      }
      if (analysis != "CYTOTOX") {
        formatted_data[[i]] <- setNames(
          data.frame(
            ColorCode = color_names[i],
            Color = sprintf("<div style='background-color: %s; padding: 10px; margin-right:20px;'></div>", colors[i]),
            Values = formatted_output,
            exp_condition,
            time_output
          ),
          c("ColorCode", "Color", "Values", value1, value2)
        )
      } else 
        formatted_data[[i]] <- setNames(
          data.frame(
            ColorCode = color_names[i],
            Color = sprintf("<div style='background-color: %s; padding: 10px; margin-right:20px;'></div>", colors[i]),
            Values = formatted_output,
            exp_condition,
            time_output,
            rep_output
          ),
          c("ColorCode", "Color", "Values", value1, value2, "ReplicateNumber")
        )
    } else {
      formatted_data[[i]] <- setNames(
        data.frame(
          ColorCode = color_names[i],
          Color = sprintf("<div style='background-color: %s; padding: 10px; margin-right:20px;'></div>", colors[i]),
          Values = "No matching indices found.",
          "-", "-"
        ),
        c("ColorCode", "Color", "Values", value1, value2)
      )
    }
  }
  return(do.call(rbind, formatted_data))
}

updateTable <- function(position, analysis, info, data, result, flag) {
  req(info) 
  
  selected_row <- info$row
  selected_col <- info$col
  new_value <- info$value
  
  # change the exp_condition column to ENDOC or sample_name to ELISA and CYTOTOX
  if (selected_col == 4) {
    color_code <- data[selected_row, "ColorCode"]
    
    if (!is.na(color_code) && color_code != "" && color_code != "white" && color_code != "#FFFFFF") {
      analysis_lower <- tolower(analysis)
      matching_indices <- which(result[[paste0(analysis, "cell_COLOR")]] == color_code, arr.ind = TRUE)
      
      if (nrow(matching_indices) > 0) {
        current_values <- c()
        
        apply(matching_indices, 1, function(idx) {
          current_values <<- c(current_values, result[[paste0("Initdata")]][idx["row"], idx["col"]])
          old_value_key <- names(flag[[paste0("EXPcol")]])[names(flag[[paste0("EXPcol")]]) == result[[paste0(analysis, "cell_COLOR")]][idx["row"], idx["col"]]]
          
          if (length(old_value_key) > 0) {
            flag[[paste0("EXPcol")]][new_value] <- flag[[paste0("EXPcol")]][old_value_key]
            flag[[paste0("EXPcol")]] <- flag[[paste0("EXPcol")]][!names(flag[[paste0("EXPcol")]]) %in% old_value_key]
            assign(paste0("Flags", analysis), flag, envir = .GlobalEnv)
          }
          
          result[[paste0(analysis, "cell_COLOR")]][idx["row"], idx["col"]] <- new_value
          # if ELISA or CYTOTOX modify SN otherwise modify EXP
          if (analysis == "ELISA" || analysis == "CYTOTOX") {
            result[[paste0(analysis, "cell_SN")]][idx["row"], idx["col"]] <- new_value
          } else if (analysis == "ENDOC") {
            result[[paste0(analysis, "cell_EXP")]][idx["row"], idx["col"]] <- new_value
          }
          
          assign(paste0(analysis_lower, "Result"), result, envir = .GlobalEnv)
        })
        
        if (!new_value %in% flag[[paste0("AllExp")]]) {
          flag[[paste0("AllExp")]] <- unique(c(flag[[paste0("AllExp")]], new_value))
          assign(paste0("Flags", analysis), flag, envir = .GlobalEnv)
        }
      }
    }
  } 
  # change the exp_condition column to ENDOC or sample_name to ELISA and CYTOTOX
  else if (selected_col == 5) {
    color_code <- data[selected_row, "ColorCode"]
    req(color_code != "", color_code != "white", color_code != "#FFFFFF")
    
    analysis_lower <- tolower(analysis)
    matching_indices <- which(result[[paste0(analysis, "cell_COLOR")]] == color_code, arr.ind = TRUE)
    num_matches <- nrow(matching_indices)
    
    # operation to set the unfilled values to NA and save only the modified position
    processed_value <- gsub(" -  - ", " - NA - ", new_value)
    processed_value <- sub("^ - ", "NA - ", processed_value)
    processed_value <- sub(" - $", " - NA", processed_value)
    
    new_values <- strsplit(processed_value, " - ", fixed = TRUE)[[1]]
    new_values[new_values == ""] <- NA  
    
    current_values <- new_values  
    
    if (length(new_values) != num_matches) {
      session$sendCustomMessage(type = "errorNotification", 
                                message = "Number of values does not match the number of matches.")
    } else {
      for (i in seq_along(matching_indices[, "row"])) {
        if (!is.na(new_values[i]) && new_values[i] != "" && new_values[i] != "NA") {
          # if ELISA or CYTOTOX modify SN otherwise modify EXP
          if (analysis == "ELISA" || analysis == "CYTOTOX") {
            result[[paste0(analysis, "cell_EXP")]][matching_indices[i, "row"], matching_indices[i, "col"]] <- new_values[i]
          } else {
            result[[paste0(analysis, "cell_TIME")]][matching_indices[i, "row"], matching_indices[i, "col"]] <- new_values[i]
          }
        } 
      }
      assign(paste0(analysis_lower, "Result"), result, envir = .GlobalEnv)
    }
  # column 6 is used only by cytotoxicity for ReplicateNumber 
  } else if (selected_col == 6) {
    color_code <- data[selected_row, "ColorCode"]
    req(color_code != "", color_code != "white", color_code != "#FFFFFF")
    
    analysis_lower <- tolower(analysis)
    matching_indices <- which(result[[paste0(analysis, "cell_COLOR")]] == color_code, arr.ind = TRUE)
    num_matches <- nrow(matching_indices)
    
    # operation to set the unfilled values to NA and save only the modified position
    processed_value <- gsub(" -  - ", " - NA - ", new_value)
    processed_value <- sub("^ - ", "NA - ", processed_value)
    processed_value <- sub(" - $", " - NA", processed_value)
    
    new_values <- strsplit(processed_value, " - ", fixed = TRUE)[[1]]
    new_values[new_values == ""] <- NA  
    
    current_values <- new_values  
    
    if (length(new_values) != num_matches) {
      session$sendCustomMessage(type = "errorNotification", 
                                message = "Number of values does not match the number of matches.")
    } else {
        for (i in seq_along(matching_indices[, "row"])) {
          if (!is.na(new_values[i]) && new_values[i] != "" && new_values[i] != "NA") {
            if (analysis == "CYTOTOX") {
              result[[paste0(analysis, "cell_REP")]][matching_indices[i, "row"], matching_indices[i, "col"]] <- new_values[i]
          } 
        }
        assign(paste0(analysis_lower, "Result"), result, envir = .GlobalEnv)
      } 
    }
  }
  return(paste("Updated values: ", new_value))
}