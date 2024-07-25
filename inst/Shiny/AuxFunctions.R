updateMultiValues = function(datapaths,result){
  for(dpath in 1:length(datapaths)){
    mess <- readRDS(datapaths[dpath])
    names(mess) -> namesRes
    if("Flags"%in% namesRes) namesRes = namesRes[ namesRes != "Flags"]
    
    if(all(namesRes %in% c("Table","TableTest", "Plot"))){
      result$Stat[[dpath]] <- mess$Table %>% mutate(DataSet = dpath)
    } else if(all(namesRes %in% names(wbquantResult)) || all(namesRes %in% names(DataAnalysisModule))){
      result$WB[[dpath]] <- mess$AdjRelDensity %>% mutate(DataSet = dpath)
    } else if(all(namesRes %in% names(pcrResult)) || all(namesRes %in% names(DataAnalysisModule))){
      result$PCR[[dpath]]  <- mess
    } else if(all(namesRes %in% names(endocResult)) || all(namesRes %in% names(DataAnalysisModule))){
      result$ENDOC[[dpath]]  <- mess
    } else if(all(namesRes %in% names(elisaResult)) || all(namesRes %in% names(DataAnalysisModule))){
      result$ELISA[[dpath]]  <- mess
    } else if(all(namesRes %in% names(cytotoxResult)) || all(namesRes %in% names(DataAnalysisModule))){
      result$CYTOTOX[[dpath]]  <- mess
    } else if(all(namesRes %in% names(ifResult)) || all(namesRes %in% names(DataAnalysisModule))){
      result$IF[[dpath]]  <- mess
    } else if(all(namesRes %in% names(facsResult)) || all(namesRes %in% names(DataAnalysisModule))){
      result$FACS[[dpath]]  <- mess
    }else{
      showAlert("Error", paste(mess[["message"]],"\n The file must be RDs saved through the Data Analysis module."), "error", 5000)
      manageSpinner(FALSE)
      return()
    }
    result$Flag <- TRUE
  }
}

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
         "BCA" = {
           FlagsBCA <- reactiveValues(cellCoo = NULL,
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
           result$BCAcell_EXP <- NULL
           result$BCAcell_SN <- NULL
           result$BCAcell_COLOR <- NULL
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
         "IF" = {
           for(l in names(result))
             result[[l]] = NULL
           for(l in names(flags))
             flags[[l]] = NULL
           
           updateSelectInput(
             inputId = "IF_TTestvariable",
             choices = "",selected = ""
           )
           output$IFtable_stat <- renderDT({data.frame()})
           output$IFtable <- renderDT({data.frame()})
           output$IFsummariseMean <- renderDT({data.frame()})
           output$IFsummarise_plot <- renderPlot({}) 
           
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
      else {LoadImage(filename)}
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
    "ExcelMulti" = {
      result <- list(data = list(), error = NULL)
      filenames <- filename
      
      for(filename in filenames){
        if(!isFileUploaded || !file.exists(filename)) {
          return (list(message = "Please, select an Excel file", call = ""))
        } else if(tolower(tools::file_ext(filename)) != "xls" && tolower(tools::file_ext(filename)) != "xlsx") {
          return(list(message = "Please, upload a file with a .xls or .xlsx extension.", call = ""))
        } 
      }
      
      xmulti = tryCatch(
        {
          do.call(rbind, 
                  lapply(filenames,function(filename){
                    x = readxl::read_excel(filename, col_names = colname)
                    x <- x %>% na.omit()
                    colnames(x) = x[1,]
                    x = x[-1,]
                    return(x)
                  }) 
          )
        }, 
        error = function(e) list(message = "The matrix in the Excel files have different column names", call = "") )
      
      return(xmulti)
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
           insertPlot(wb = wb, sheet = "WBimage and protein bands",
                      fileType = "tiff",
                      units = "in",
                      dpi = 600)
           
           startRow <- 22
           writeDataTable(wb, PanelStructures$data, sheet = "WBimage and protein bands", startRow = startRow, startCol = 1)
           
           addWorksheet(wb, "Plot")
           print(ResultList[["Plots"]])
           insertPlot(wb = wb, sheet="Plot",
                      fileType = "tiff",
                      units = "in",
                      dpi = 600, width = 10,height = 8)
           
           addWorksheet(wb, "Truncated Plot")
           print(ResultList[["TruncatedPlots"]])
           insertPlot(wb = wb, sheet="Truncated Plot",
                      fileType = "tiff",
                      units = "in",
                      dpi = 600, width = 10,height = 8)
           
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
           
           addWorksheet(wb,"AdjRelDensity")
           writeDataTable(wb, sheet = "AdjRelDensity", ResultList[["AdjRelDensity"]])
           
           addWorksheet(wb,"Barplot AdjRelDensity")
           if(!is.null( ResultList[["AdjRelDensity"]])){
             print(
               ResultList[["AdjRelDensity"]]   %>%
                 ggplot() +
                 geom_bar(aes(x = SampleName,
                              y = AdjRelDensity,
                              fill = SampleName ),
                          stat = "identity" ) +
                 theme_bw()
             )
           }
           insertPlot(wb, sheet = "Barplot AdjRelDensity",
                      fileType = "tiff",
                      units = "in",
                      dpi = 600, width = 10,height = 8)
         }, 
         "RT-qPCR" = {
           wb <- createWorkbook("RTqPCR")
           
           addWorksheet(wb,"Table")
           writeDataTable(wb, sheet = "Table", ResultList[["Initdata"]])
           
           
           if(!is.null(ResultList$AllGenesFoldChangePlot)){
             addWorksheet(wb,"AllGenes")
             writeDataTable(wb,ResultList[["AllGenesFoldChangeTable"]], sheet="AllGenes")
             
             print(ResultList[["AllGenesFoldChangePlot"]])
             insertPlot(wb = wb,  sheet="AllGenes",
                        startCol=dim(ResultList[["AllGenesFoldChangeTable"]])[2]+ 2,
                        fileType = "tiff",
                        units = "in",
                        dpi = 600,width = 14,height = 6)
             
           }
           
           for(j in 1:length(names(ResultList[["plotPCR"]])) ){
             i = names(ResultList[["plotPCR"]])[j]
             
             addWorksheet(wb,paste0("Gene ", j) )
             writeDataTable(wb,ResultList[["plotPCR"]][[i]]$table, sheet=paste0("Gene ", j))
             
             print(ResultList[["plotPCR"]][[i]]$plot + labs(title = i))
             insertPlot(wb = wb,  sheet=paste0("Gene ", j),
                        startCol=dim(ResultList[["plotPCR"]][[i]]$table)[2]+ 2,
                        fileType = "tiff",
                        units = "in",
                        dpi = 600,width = 14,height = 6)
             
           }
           
           addWorksheet(wb,"Points Plot")
           plotList = ResultList[["plotPCR"]]
           df = do.call(rbind, lapply(plotList, `[[`, 2)) %>% filter(DDCT == 0)
           
           print(
             ggplot(df,aes(x = Gene, y = -DDCT)) +
               geom_point() +
               facet_wrap(~HousekGene, ncol = 1)+
               theme_bw()+
               labs(x="",y= "Log2(Q)", title = "")+
               theme(axis.text.x=element_blank(), 
                     axis.ticks.x=element_blank())
           )
           insertPlot(wb = wb,  sheet="Points Plot",
                      fileType = "tiff",
                      units = "in",
                      dpi = 600)
           
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
           
           # Se esiste anche dataFinal e vuoi scriverlo su un altro foglio
           if (!is.null(ResultList$dataFinal) && is.data.frame(ResultList$dataFinal)) {
             addWorksheet(wb, "Results Analysis")  # Aggiunge un foglio per l'analisi finale
             writeDataTable(wb, ResultList$dataFinal, sheet = "Results Analysis")  # Scrive i dati di analisi finale
             print("dataFinal scritto nel foglio Excel")
           } else {
             print("dataFinal non disponibile o non è un data.frame")
           }
         },
         "ELISA" = {
           ## Create a new workbook
           wb <- createWorkbook("ELISA")
           
           ## initial data
           addWorksheet(wb,"TablePlot")
           writeDataTable(wb, sheet = "TablePlot", ResultList[["TablePlot"]]$x$data)
           
           ## Linear regression analysis
           addWorksheet(wb,"standard curve")
           standcurve = ResultList[["Tablestandcurve"]]
           lmStancurve = ResultList[["Regression"]]$data
           print(ResultList[["Regression"]]$plot)
           
           writeDataTable(wb,standcurve, sheet="standard curve")
           insertPlot(wb = wb,  sheet="standard curve",
                      startCol=dim(standcurve)[2]+ 2,
                      fileType = "tiff",
                      units = "in",
                      dpi = 600)
           
           ## Analysis
           addWorksheet(wb,"Analysis")
           writeDataTable(wb,ResultList[["dataFinal"]], sheet="Analysis")
           print(ResultList[["dataFinal"]] %>%
                   ggplot( aes(x = Time, y = Quantification,
                               fill= Experiment, group = Experiment ) )+
                   geom_bar(position = "dodge",stat = "identity")+
                   theme_bw()+
                   labs(x = "Time", col = "Experiments",
                        y = "Average quantifications obtained\n from the lm "))
           insertPlot(wb = wb,  sheet="Analysis",
                      startCol=dim(ResultList[["dataFinal"]] )[2]+ 2,
                      fileType = "tiff",
                      units = "in",
                      dpi = 600)
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
           if (!is.null(ResultList$Barplot) ) {
             
             print(ResultList$Barplot)
             insertPlot(wb,  sheet = "Final data",
                        startCol=dim(ResultList[["dataFinal"]] )[2]+ 2,
                        fileType = "tiff",
                        units = "in",
                        dpi = 600,width = 8,height = 6)  
           }
           
         },
         "BCA" = {
           wb <- createWorkbook("BCA")  
           
           if (!is.null(ResultList$Initdata) && is.data.frame(ResultList$Initdata)) {
             addWorksheet(wb, "TablePlot")  # Aggiunge un nuovo foglio al workbook
             writeDataTable(wb, ResultList$Initdata, sheet = "TablePlot")  # Scrive i dati nel foglio
             print("Initdata scritto nel foglio Excel")
           } else {
             print("Errore: Initdata non disponibile o non è un data.frame")
           }
           
           ## Linear regression analysis
           if (!is.null(ResultList$Regression)) {
             addWorksheet(wb,"standard curve")
             standcurve = ResultList[["Tablestandcurve"]]
             lmStancurve = ResultList[["Regression"]]$data
             print(ResultList[["Regression"]]$plot)
             
             writeDataTable(wb,standcurve, sheet="standard curve")
             insertPlot(wb = wb,  sheet="standard curve",
                        startCol=dim(standcurve)[2]+ 2,
                        fileType = "tiff",
                        units = "in",
                        dpi = 600)
           }
           # Se esiste anche dataQuant e vuoi scriverlo su un altro foglio
           if (!is.null(ResultList$dataQuant) && is.data.frame(ResultList$dataQuant)) {
             addWorksheet(wb, "Quantification Analysis")  # Aggiunge un foglio per l'analisi finale
             writeDataTable(wb, ResultList$dataQuant, sheet = "Quantification Analysis")  # Scrive i dati di analisi finale
             print("dataQuant scritto nel foglio Excel")
           } else {
             print("dataQuant non disponibile o non è un data.frame")
           }
           # Se esiste anche dataFinal e vuoi scriverlo su un altro foglio
           if (!is.null(ResultList$dataFinal) && is.data.frame(ResultList$dataFinal)) {
             addWorksheet(wb, "Results Analysis")  # Aggiunge un foglio per l'analisi finale
             writeDataTable(wb, ResultList$dataFinal, sheet = "Results Analysis")  # Scrive i dati di analisi finale
             print("dataFinal scritto nel foglio Excel")
           } else {
             print("dataFinal non disponibile o non è un data.frame")
           }
         },
         "IF" = {
           wb <- createWorkbook("IF")  
           
           if (!is.null(ResultList$Initdata) && is.data.frame(ResultList$Initdata)) {
             addWorksheet(wb, "TablePlot")  # Aggiunge un nuovo foglio al workbook
             writeDataTable(wb, ResultList$Initdata, sheet = "TablePlot")  # Scrive i dati nel foglio
             print("Initdata scritto nel foglio Excel")
           } else {
             print("Errore: Initdata non disponibile o non è un data.frame")
           }
           # Se esiste anche dataFinal e vuoi scriverlo su un altro foglio
           if (!is.null(ResultList$dataFinal) && is.data.frame(ResultList$dataFinal)) {
             addWorksheet(wb, "Results Analysis")  # Aggiunge un foglio per l'analisi finale
             writeDataTable(wb, ResultList$dataFinal, sheet = "Results Analysis")  # Scrive i dati di analisi finale
             print("dataFinal scritto nel foglio Excel")
           } else {
             print("dataFinal non disponibile o non è un data.frame")
           }
           
           # Se esiste anche SubStatData e vuoi scriverlo su un altro foglio
           if (!is.null(ResultList$SubStatData) && is.data.frame(ResultList$SubStatData)) {
             addWorksheet(wb, "Data selected")  # Aggiunge un foglio per l'analisi finale
             writeDataTable(wb, ResultList$SubStatData, sheet = "Data selected")  # Scrive i dati di analisi finale
             print("SubStatData scritto nel foglio Excel")
           } else {
             print("SubStatData non disponibile o non è un data.frame")
           }
           
           # Se esiste anche TTestData e vuoi scriverlo su un altro foglio
           if (!is.null(ResultList$TTestData) && is.data.frame(ResultList$TTestData)) {
             addWorksheet(wb, "T-Test Analysis")  # Aggiunge un foglio per l'analisi finale
             writeDataTable(wb, ResultList$TTestData, sheet = "T-Test Analysis")  # Scrive i dati di analisi finale
             print("TTestData scritto nel foglio Excel")
           } else {
             print("TTestData non disponibile o non è un data.frame")
           }
           
           # Se esiste anche resplot e vuoi scriverlo su un altro foglio
           if (!is.null(ResultList$resplot)) {
             addWorksheet(wb, "Bar plot") 
             print(ResultList$resplot)
             insertPlot(wb = wb,  sheet="Bar plot",
                        fileType = "tiff",
                        units = "in",
                        dpi = 600,width = 8,height = 6)
             print("resplot scritto nel foglio Excel")
           } else {
             print("resplot non disponibile o non è un data.frame")
           }
         }
  )
  
  saveWorkbook(wb, filename)  # Salva il workbook
  return(1)  # Restituisce 1 per indicare il successo
}

tableExcelColored = function(session, output,Result, FlagsExp, type){
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
  if (analysis %in% c("ELISA","BCA","IF") ){
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
        if (analysis == "ELISA")
          val <- result$ELISAcell_EXP[idx["row"], idx["col"]]
        else if (analysis == "IF")
          val <- result$IFcell_EXP[idx["row"], idx["col"]]
        else if (analysis == "BCA")
          val <- result$BCAcell_EXP[idx["row"], idx["col"]]
        else val <- result$ENDOCcell_TIME[idx["row"], idx["col"]]
        
        if (!is.na(val) && !is.null(val) && val != "") val else ""
      })
      
      time_output <- paste(unlist(time_values), collapse = " - ")
      
      exp_values <- apply(matching_indices, 1, function(idx) {
        result[[column_EXP]][idx["row"], idx["col"]]
      })
      
      if (length(unique(exp_values)) == 1) {
        exp_condition <- ifelse(exp_values[1] == "" || is.na(exp_values[1]), "-", exp_values[1])
      } else {
        exp_condition <- "No matching between values"
      }
      
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
  
  # change the exp_condition column to ENDOC or sample_name to ELISA  
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
          # if ELISA, modify SN otherwise modify EXP
          if (analysis %in% c("ELISA","BCA","IF") ) {
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
  # change the time column to ENDOC or exp_condition to ELISA  
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
          # if ELISA, modify EXP otherwise modify TIME
          if (analysis == "ELISA") {
            result[[paste0(analysis, "cell_EXP")]][matching_indices[i, "row"], matching_indices[i, "col"]] <- new_values[i]
          } else {
            result[[paste0(analysis, "cell_TIME")]][matching_indices[i, "row"], matching_indices[i, "col"]] <- new_values[i]
          }
        } 
      }
      assign(paste0(analysis_lower, "Result"), result, envir = .GlobalEnv)
    }
  }
  return(paste("Updated values: ", new_value))
}

updateSelectizeUI <- function(maxDepth) {
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
  return(rowContent)
}

escapeRegex <- function(string) {
  gsub("([\\\\^$.*+?()[{\\]|-])", "\\\\\\1", string)
}

loadDrop <- function(facsResult, FlagsFACS, session) {
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

UploadRDs = function(Flag, session, output,
                     DataAnalysisModule,
                     Result, FlagsExp,PanelStructures=NULL){
  
  if(Flag == "WB"){
    
    for(nameList in names(Result)) 
      Result[[nameList]] <- DataAnalysisModule$wbResult[[nameList]]
    
    for(nameList in names(DataAnalysisModule$wbResult$Flags)) 
      FlagsExp[[nameList]] <- DataAnalysisModule$wbResult$Flags[[nameList]]
    
    # update image WB
    if(!is.null(Result$Im)){
      im = Result$Im$WB
      output$TifPlot2 <- renderPlot({
        plot(c(1,dim(im)[2]),c(1,dim(im)[1]), type='n',ann=FALSE)
        rasterImage(im,1,1,dim(im)[2],dim(im)[1])
        if (nrow(Result$Planes) > 0) {
          r <- Result$Planes
          rect(r$xmin, r$ymin, r$xmax, r$ymax, border = "red")
        }
      })
      
      output$PlanesStructureTable <- renderDT(
        Result$Planes,
        server = FALSE,
        editable = list(target = "cell", 
                        disable = list(columns = 1:4)),
        options = list(lengthChange = FALSE, autoWidth = TRUE),
        rownames= FALSE
      )
      
      PanelStructures$data = Result$Planes
      
    }
    
    # update lanes
    
    if(!is.null(Result$PanelsValue)){
      updateSelectInput(session = session, "LaneChoice",
                        choices = unique(Result$PanelsValue$ID),
                        selected = ""
      )
    }
    
    if(!is.null(Result$TruncatedPlots)){
      output$DataPlot <- renderPlot({Result$TruncatedPlots})
      FlagsExp$LanesCut = T
    }else if(!is.null(Result$Plots)){
      output$DataPlot <- renderPlot({Result$Plots})
      FlagsExp$LanesCut = T
    }
    
    output$AUC <- renderDT({
      Result$AUCdf %>% 
        dplyr::select(SampleName,Truncation, AUC)
    },
    selection = 'none',
    #session = session,
    rownames= FALSE
    )
    
    # output$AUC <- renderTable({
    #   wbResult$AUCdf  %>% dplyr::select(Lane,Truncation, AUC)
    # },width = "100%")
    
    # change pannel
    updateTabsetPanel(session = session, "SideTabs",
                      selected = "grey")
    
  }else if(Flag == "WBquant"){
    
    for(nameList in names(DataAnalysisModule$wbquantResult)) 
      Result[[nameList]] <- DataAnalysisModule$wbquantResult[[nameList]]
    
    for(nameList in names(DataAnalysisModule$wbquantResult$Flags)) 
      FlagsExp[[nameList]] <- DataAnalysisModule$wbquantResult$Flags[[nameList]]
    
    # change pannel
    updateTabsetPanel(session = session, "SideTabs",
                      selected = "quantification")
    
  }else if(Flag == "PRCC"){
    
    for(nameList in names(DataAnalysisModule$pcrResult)) 
      Result[[nameList]] <- DataAnalysisModule$pcrResult[[nameList]]
    
    for(nameList in names(DataAnalysisModule$pcrResult$Flags)) 
      FlagsExp[[nameList]] <- DataAnalysisModule$pcrResult$Flags[[nameList]]
    
    choices = ""
    selected = rep("",3)
    
    if(!is.null(Result$Initdata)){
      choices = c("",colnames(Result$Initdata))
    }
    
    if(!is.null(Result$selectPCRcolumns)){
      selected = Result$selectPCRcolumns
    }
    
    updateSelectInput(session = session,"PCR_gene",
                      choices = choices,
                      selected = selected[1]
    )
    updateSelectInput(session = session,"PCR_sample",
                      choices = choices,
                      selected = selected[2]
    )
    updateSelectInput(session = session,"PCR_value",
                      choices = choices,
                      selected = selected[3]
    )
    
    
    if(!is.null(Result$PCRnorm))
      FlagsExp$norm = T
    
    if(!is.null(Result$PCRbaseline))
      FlagsExp$baseline = T
    
    # change pannel
    updateTabsetPanel(session = session, "SideTabs",
                      selected = "uploadPCR")
    
    
  }
  else if(Flag == "ENDOC"){
    
    for(nameList in names(DataAnalysisModule$endocResult)) 
      Result[[nameList]] <- DataAnalysisModule$endocResult[[nameList]]
    
    for(nameList in names(DataAnalysisModule$endocResult$Flags)) 
      FlagsExp[[nameList]] <- DataAnalysisModule$endocResult$Flags[[nameList]]
    
    if(!is.null(Result$TablePlot)){
      output$ENDOCmatrix <- renderDT(
        Result$TablePlot,
        server = FALSE
      )
    }
    
    if(!is.null(Result$ENDOCcell_TIME)){
      
      updateSelectizeInput(inputId = "ENDOCcell_TIME",session = session,
                           choices = unique(c(Result$ENDOCcell_TIME))
      )
      
      updateSelectizeInput(inputId = "ENDOCcell_SN",session = session,
                           choices = unique(c(Result$ENDOCcell_SN))
      )
      
      FlagsExp$AllExp = unique(c(Result$ENDOCcell_SN))
    }
    
    # change pannel
    updateTabsetPanel(session = session, "SideTabs",
                      selected = "uploadENDOC")
    
  } 
  else if(Flag == "BCA"){
    
    for(nameList in names(DataAnalysisModule$bcaResult)) 
      Result[[nameList]] <- DataAnalysisModule$bcaResult[[nameList]]
    
    for(nameList in names(DataAnalysisModule$bcaResult$Flags)) 
      FlagsExp[[nameList]] <- DataAnalysisModule$bcaResult$Flags[[nameList]]
    
    if(!is.null(Result$TablePlot)){
      output$BCAmatrix <- renderDT(
        Result$TablePlot,
        server = FALSE
      )
    }
    
    if(!is.null(Result$BCAcell_EXP)){
      
      updateSelectizeInput(inputId = "BCAcell_EXP",
                           session = session,
                           choices = c("",unique(c(Result$BCAcell_EXP)) ),
                           selected = ""
      )
      
      updateSelectizeInput(inputId = "BCAcell_SN",
                           session =session,
                           choices = c("",unique(c(Result$BCAcell_SN)) ),
                           selected = ""
      )
      
      FlagsExp$AllExp = unique(c(Result$BCAcell_SN))
    }
    if(!is.null(Result$MapBaseline)){
      FlagsExp$BASEselected = unique(Result$MapBaseline$Baseline)
    }
    if(!is.null(Result$MapBlanche)){
      FlagsExp$BLANCHEselected = unique(Result$MapBlanche$Blanche)
    }
    if(!is.null(Result$Tablestandcurve)){
      output$BCA_Table_stdcurve <- DT::renderDataTable({
        DT::datatable( Result$Tablestandcurve,
                       selection = 'none',
                       editable = list(target = "cell",
                                       disable = list(columns = 0:2) ),
                       options = list(lengthChange = FALSE, autoWidth = TRUE),
                       rownames= FALSE
        )
      })
      
      if(length(FlagsExp$AllExp) > 1){
        exp = FlagsExp$AllExp
        exp = exp[exp != ""]
        
        bool.tmp = exp %in% unique(c(FlagsExp$BLANCHEselected,FlagsExp$BASEselected))
        if( length(bool.tmp) > 0  )
          exp = exp[!bool.tmp]
        
        updateSelectizeInput(session = session,
                             inputId = "BCA_standcurve",
                             choices = exp,
                             selected = unique(Result$Tablestandcurve$exp) )
        
        FlagsExp$STDCselected= unique(Result$Tablestandcurve$exp)
      }
    }
    
    ### updating check boxes
    if(!is.null(FlagsExp$BLANCHEselected)){
      exp = FlagsExp$AllExp
      exp = exp[exp != ""]
      
      bool.tmp = exp %in% unique(c(FlagsExp$STDselected,FlagsExp$BASEselected))
      if( length(bool.tmp) > 0  )
        exp = exp[!bool.tmp]
      
      updateCheckboxGroupInput(session = session,
                               inputId = "BCA_blanches",
                               choices = exp,
                               selected = unique(FlagsExp$BLANCHEselected) )
    }
    if(!is.null(FlagsExp$BASEselected)){
      exp = FlagsExp$AllExp
      exp = exp[exp != ""]
      
      bool.tmp = exp %in% unique(c(FlagsExp$STDselected,FlagsExp$BLANCHEselected))
      if( length(bool.tmp) > 0  )
        exp = exp[!bool.tmp]
      
      updateCheckboxGroupInput(session = session,
                               inputId = "BCA_baselines",
                               choices = exp,
                               selected = unique(FlagsExp$BASEselected) )
    }
    ###
    if(!is.null(Result$Regression)){
      Result$Regression$data -> lmStancurve
      Result$Tablestandcurve -> standcurve
      
      output$BCAregression <- renderPlot(
        Result$Regression$plot
      )
    }
    
    # change pannel
    updateTabsetPanel(session = session, "SideTabs",
                      selected = "uploadBCA")
    
  }
  else if(Flag == "ELISA"){
    
    for(nameList in names(DataAnalysisModule$elisaResult)) 
      Result[[nameList]] <- DataAnalysisModule$elisaResult[[nameList]]
    
    for(nameList in names(DataAnalysisModule$elisaResult$Flags)) 
      FlagsExp[[nameList]] <- DataAnalysisModule$elisaResult$Flags[[nameList]]
    
    if(!is.null(Result$TablePlot)){
      output$ELISAmatrix <- renderDT(
        Result$TablePlot,
        server = FALSE
      )
    }
    
    if(!is.null(Result$ELISAcell_EXP)){
      
      updateSelectizeInput(inputId = "ELISAcell_EXP",
                           session = session,
                           choices = unique(c(Result$ELISAcell_EXP))
      )
      
      updateSelectizeInput(inputId = "ELISAcell_SN",
                           session =session,
                           choices = unique(c(Result$ELISAcell_SN))
      )
      
      FlagsExp$AllExp = unique(c(Result$ELISAcell_SN))
    }
    if(!is.null(Result$MapBaseline)){
      FlagsExp$BASEselected = unique(Result$MapBaseline$Baseline)
    }
    if(!is.null(Result$MapBlanche)){
      FlagsExp$BLANCHEselected = unique(Result$MapBlanche$Blanche)
    }
    if(!is.null(Result$Tablestandcurve)){
      output$ELISA_Table_stdcurve <- DT::renderDataTable({
        DT::datatable( Result$Tablestandcurve,
                       selection = 'none',
                       editable = list(target = "cell",
                                       disable = list(columns = 0:2) ),
                       options = list(lengthChange = FALSE, autoWidth = TRUE),
                       rownames= FALSE
        )
      })
      
      if(length(FlagsExp$AllExp) > 1){
        exp = FlagsExp$AllExp
        exp = exp[exp != ""]
        
        bool.tmp = exp %in% unique(c(FlagsExp$BLANCHEselected,FlagsExp$BASEselected))
        if( length(bool.tmp) > 0  )
          exp = exp[!bool.tmp]
        
        updateSelectizeInput(session = session,
                             inputId = "ELISA_standcurve",
                             choices = exp,
                             selected = unique(Result$Tablestandcurve$exp) )
        
        FlagsExp$STDCselected= unique(Result$Tablestandcurve$exp)
      }
    }
    
    ### updating check boxes
    if(!is.null(FlagsExp$BLANCHEselected)){
      exp = FlagsExp$AllExp
      exp = exp[exp != ""]
      
      bool.tmp = exp %in% unique(c(FlagsExp$STDselected,FlagsExp$BASEselected))
      if( length(bool.tmp) > 0  )
        exp = exp[!bool.tmp]
      
      updateCheckboxGroupInput(session = session,
                               inputId = "ELISA_blanches",
                               choices = exp,
                               selected = unique(FlagsExp$BLANCHEselected) )
    }
    if(!is.null(FlagsExp$BASEselected)){
      exp = FlagsExp$AllExp
      exp = exp[exp != ""]
      
      bool.tmp = exp %in% unique(c(FlagsExp$STDselected,FlagsExp$BLANCHEselected))
      if( length(bool.tmp) > 0  )
        exp = exp[!bool.tmp]
      
      updateCheckboxGroupInput(session = session,
                               inputId = "ELISA_baselines",
                               choices = exp,
                               selected = unique(FlagsExp$BASEselected) )
    }
    ###
    if(!is.null(Result$Regression)){
      Result$Regression$data -> lmStancurve
      Result$Tablestandcurve -> standcurve
      
      output$ELISAregression <- renderPlot(
        Result$Regression$plot
      )
    }
    
    # change pannel
    updateTabsetPanel(session = session, "SideTabs",
                      selected = "uploadELISA")
    
  }
  else if(Flag == "CYTOTX"){
    
    for(nameList in names(DataAnalysisModule$cytotoxResult)) 
      Result[[nameList]] <- DataAnalysisModule$cytotoxResult[[nameList]]
    
    for(nameList in names(DataAnalysisModule$cytotoxResult$Flags)) 
      FlagsExp[[nameList]] <- DataAnalysisModule$cytotoxResult$Flags[[nameList]]
    
    if(!is.null(Result$TablePlot)){
      output$CYTOTOXmatrix <- renderDT(
        Result$TablePlot,
        server = FALSE
      )
    }
    
    if(!is.null(Result$CYTOTOXcell_EXP)){
      
      updateSelectizeInput(inputId = "CYTOTOXcell_EXP",
                           session = session,
                           choices = unique(c(Result$CYTOTOXcell_EXP))
      )
      
      updateSelectizeInput(inputId = "CYTOTOXcell_SN",
                           session =session,
                           choices = unique(c(Result$CYTOTOXcell_SN))
      )
      updateSelectizeInput(inputId = "CYTOTOXcell_REP",
                           session =session,
                           choices = unique(c(Result$CYTOTOXcell_REP))
      )
      # FlagsExp$AllExp = unique(c(Result$CYTOTOXcell_SN))
    }
    
    # change pannel
    updateTabsetPanel(session = session, "SideTabs",
                      selected = "uploadCYTOTOX")
    
  }
  else if(Flag == "FACS"){
    
    for(nameList in names(DataAnalysisModule$facsResult)) 
      Result[[nameList]] <- DataAnalysisModule$facsResult[[nameList]]
    
    for(nameList in names(DataAnalysisModule$facsResult$Flags)) 
      FlagsExp[[nameList]] <- DataAnalysisModule$facsResult$Flags[[nameList]]
    
    # facsResult = reactiveValues(
    #   Initdata= NULL,
    #   data = NULL,
    #   dataFinal = NULL,
    #   depth = NULL,
    #   depthCount = NULL,
    #   originalName = NULL,
    #   name = NULL,
    #   statistics = NULL,
    #   cells = NULL,
    #   ExpConditionDF = NULL,
    #   Barplot = NULL,
    #   StatDF = NULL
    # )
    
    
    # change pannel
    updateTabsetPanel(session, "SideTabs", selected = "tablesFACS")
    
  }
  else if(Flag == "IF"){
    
    for(nameList in names(DataAnalysisModule$ifResult)) 
      Result[[nameList]] <- DataAnalysisModule$ifResult[[nameList]]
    
    for(nameList in names(DataAnalysisModule$ifResult$Flags)) 
      FlagsExp[[nameList]] <- DataAnalysisModule$ifResult$Flags[[nameList]]
    
    
    # change pannel
    
    updateTabsetPanel(session, "SideTabs", selected = "tablesIF")
    
  }
}

# testStat.function = function(data, var = NULL) {
#   steps <- c()  
#   
#   print("Contenuto della colonna prima della conversione:")
#   print(data[, 2])
#   
#   # Pulizia dei dati: rimuove spazi bianchi e verifica di valori non numerici
#   data[, 2] <- gsub("^\\s+|\\s+$", "", data[, 2])  # Rimuove spazi bianchi
#   data[, 2] <- as.numeric(data[, 2])  # Converte a numerico
#   
#   # Shapiro-Wilk normality test
#   shapiro_test = shapiro.test(data[, 2])
#   steps <- c(steps, paste("Shapiro-Wilk test performed, p-value:", shapiro_test$p.value))
#   
#   if (shapiro_test$p.value > 0.05) {
#     steps <- c(steps, "Data is normally distributed")
#     
#     vars = data[, 1] %>% distinct() %>% pull()
#     if (length(vars) == 2) {
#       steps <- c(steps, "Two groups, performing t-test")
#       
#       combo = combn(vars, 2)
#       combo = data.frame(Var1 = combo[1, ], Var2 = combo[2, ])
#       combo = combo[combo$Var1 != combo$Var2, ]
#       
#       resTTest = do.call(rbind,
#                          lapply(1:dim(combo)[1], function(x) {
#                            sn = combo[x, ]
#                            ttest = t.test(data[data[, 1] == sn$Var1, 2],
#                                           data[data[, 1] == sn$Var2, 2])
#                            data.frame(Test = "t-test",
#                                       Condition = paste(sn$Var1, " vs ", sn$Var2),
#                                       pValue = ttest$p.value,
#                                       conf.int = paste(ttest$conf.int, collapse = ";")
#                            )
#                          })
#       )
#       
#       steps <- c(steps, "t-test completed")
#     } else if (length(vars) > 2) {
#       steps <- c(steps, "More than two groups, performing ANOVA")
#       
#       colnames(data) = c("SampleName", "Value")
#       data$SampleName <- as.factor(data$SampleName)
#       # Perform ANOVA
#       anova_model = aov(Value ~ SampleName, data = data)
#       summary(anova_model) -> a
#       print(a)
#       # Calculate group means and standard errors
#       group_stats = data %>%
#         group_by(SampleName) %>%
#         summarize(
#           mean = mean(Value),
#           sd = sd(Value),
#           n = n()
#         )
#       
#       # Calculate standard error
#       group_stats = group_stats %>%
#         mutate(se = sd / sqrt(n))
#       
#       # Calculate critical value for 95% confidence interval
#       alpha = 0.05
#       t_critical = qt(1 - alpha/2, df = df.residual(anova_model))
#       
#       # Calculate confidence intervals
#       group_stats = group_stats %>%
#         mutate(
#           ci_lower = mean - t_critical * se,
#           ci_upper = mean + t_critical * se
#         )
#       
#       # Display ANOVA results
#       resTTest = data.frame(Test = "Anova",
#                             Condition = paste(anova_model$call)[2],
#                             pValue = a[[1]]$`Pr(>F)`[1],
#                             conf.int = paste("-", collapse = ";"))
#       
#       steps <- c(steps, "ANOVA completed")
#     } else {
#       steps <- c(steps, "Invalid number of groups for t-test or ANOVA")
#       resTTest = data.frame(Test = "Error", Condition = "Invalid number of groups", pValue = NA, conf.int = NA)
#     }
#     
#     if (!is.null(var))
#       resTTest$Var = var
#     
#     return(list(resTTest = resTTest, steps = steps))
#   } else {
#     steps <- c(steps, "Data is not normally distributed")
#     
#     vars = data[, 1] %>% distinct() %>% pull()
#     if (length(vars) > 2) {
#       steps <- c(steps, "More than two groups, performing Kruskal-Wallis test")
#       
#       # Perform Kruskal-Wallis test
#       kruskal_test = kruskal.test(data[, 2] ~ data[, 1])
#       resKruskal = data.frame(Test = "Kruskal-Wallis",
#                               Condition = paste(kruskal_test$data.name),
#                               pValue = kruskal_test$p.value,
#                               conf.int = paste("-", collapse = ";"))
#       steps <- c(steps, "Kruskal-Wallis test completed")
#       return(list(resTTest = resKruskal, steps = steps))
#     } else if (length(vars) == 2) {
#       steps <- c(steps, "Two groups, performing Wilcoxon test")
#       
#       # Perform Wilcoxon test
#       wilcox_test = wilcox.test(data[data[, 1] == vars[1], 2],
#                                 data[data[, 1] == vars[2], 2])
#       resWilcox = data.frame(Test = "Wilcoxon",
#                              Condition = paste(vars[1], " vs ", vars[2]),
#                              pValue = wilcox_test$p.value,
#                              conf.int = paste("-", collapse = ";"))
#       steps <- c(steps, "Wilcoxon test completed")
#       return(list(resTTest = resWilcox, steps = steps))
#     } else {
#       steps <- c(steps, "Invalid number of groups")
#       return(list(resTTest = data.frame(Test = "Error", Condition = "Invalid number of groups", pValue = NA, conf.int = NA), steps = steps))
#     }
#   }
# }

testStat.function <- function(data) {
  steps <- ""
  step_counter <- 1
  resTTest <- NULL
  resANOVA <- NULL
  resPairwise <- NULL
  
  data$Value <- as.numeric(gsub("^\\s+|\\s+$", "", as.character(data[[2]])))
  
  shapiro_results <- data %>%
    group_by(data[[1]]) %>%
    filter(n() >= 3) %>%
    summarize(
      p.value = ifelse(
        length(unique(Value)) > 1, 
        shapiro.test(Value)$p.value, 
        NA
      ),
      message = ifelse(
        length(unique(Value)) > 1, 
        "", 
        "All values are identical or less than 3 samples"
      )
    )
  
  steps <- c(steps, paste("Step", step_counter, ". Shapiro-Wilk test performed for each group:\n"))
  step_counter <- step_counter + 1
  for (i in 1:nrow(shapiro_results)) {
    if (is.na(shapiro_results$p.value[i])) {
      steps <- c(steps, paste("  Group", shapiro_results[[1]][i], ": p-value = NA (", shapiro_results$message[i], ")\n"))
    } else {
      steps <- c(steps, paste("  Group", shapiro_results[[1]][i], ": p-value =", shapiro_results$p.value[i], "\n"))
    }
  }
  
  group_counts <- data %>%
    group_by(data[[1]]) %>%
    summarize(count = n())
  
  if (all(group_counts$count < 30) || all(shapiro_results$p.value > 0.05, na.rm = TRUE)) {
    steps <- c(steps, paste("Step ", step_counter, ". the data is normalized", "\n"))
    step_counter <- step_counter + 1 
    
    vars <- data[,1] %>% distinct() %>% pull() 
    
    if (length(vars) == 2) {
      steps <- c(steps, paste("Step ", step_counter, ". there are 2 groups, I will use t-test for analysis", "\n"))
      step_counter <- step_counter + 1  
      combo <- combn(vars, 2)
      combo <- data.frame(Var1 = combo[1,], Var2 = combo[2,])
      
      combo <- combo[combo$Var1 != combo$Var2, ]
      resTTest <- do.call(rbind,
                          lapply(1:dim(combo)[1], function(x){
                            sn <- combo[x,]
                            ttest <- t.test(data[data[,1] == sn$Var1, "Value"],
                                            data[data[,1] == sn$Var2, "Value"]) 
                            data.frame(Test = "t-test",
                                       Condition = paste(sn$Var1, " vs ", sn$Var2), 
                                       pValue = ttest$p.value,
                                       conf.int = paste(ttest$conf.int, collapse = ";")
                            )
                          })
      )
    } else if(length(vars) > 2){
      steps <- c(steps, paste("Step ", step_counter, ". groups are more than 2, I will use ANOVA for analysis", "\n"))
      step_counter <- step_counter + 1  
      colnames(data) <- c("SampleName", "Value")
      data$SampleName <- as.factor(data$SampleName)
      
      anova_model <- aov(Value ~ SampleName, data = data)
      summary(anova_model) -> a
      
      resANOVA <- data.frame(Test = "Anova",
                             Condition = paste(anova_model$call)[2], 
                             pValue = a[[1]]$`Pr(>F)`[1],
                             conf.int = paste("-", collapse = ";"))
      
      if (resANOVA$pValue < 0.05) {
        steps <- c(steps, paste("Step ", step_counter, ". ANOVA p-value <", resANOVA$pValue, ", performing pairwise t-tests", "\n"))
        step_counter <- step_counter + 1  
        
        combo <- combn(vars, 2)
        combo <- data.frame(Var1 = combo[1,], Var2 = combo[2,])
        
        combo <- combo[combo$Var1 != combo$Var2, ]
        resPairwise <- do.call(rbind,
                               lapply(1:dim(combo)[1], function(x){
                                 sn <- combo[x,]
                                 ttest <- t.test(data[data[,1] == sn$Var1, "Value"],
                                                 data[data[,1] == sn$Var2, "Value"]) 
                                 data.frame(Test = "t-test",
                                            Condition = paste(sn$Var1, " vs ", sn$Var2), 
                                            pValue = ttest$p.value,
                                            conf.int = paste(ttest$conf.int, collapse = ";")
                                 )
                               })
        )
      } else {
        steps <- c(steps, paste("Step ", step_counter, ". ANOVA p-value >=", resANOVA$pValue, ", no pairwise t-tests performed", "\n"))
      }
    }
    
    return(list(resTTest = resTTest, test = resANOVA, pairwise = resPairwise, steps = steps))
  } else {
    steps <- c(steps, paste("Step ", step_counter, ". the data is not normalized", "\n"))
    step_counter <- step_counter + 1 
    
    vars <- data[,1] %>% distinct() %>% pull()
    
    if (length(vars) == 2) {
      steps <- c(steps, paste("Step ", step_counter, ". there are 2 groups, I will use Wilcoxon test for analysis", "\n"))
      step_counter <- step_counter + 1  
      combo <- combn(vars, 2)
      combo <- data.frame(Var1 = combo[1,], Var2 = combo[2,])
      
      combo <- combo[combo$Var1 != combo$Var2, ]
      resTTest <- do.call(rbind,
                          lapply(1:dim(combo)[1], function(x){
                            sn <- combo[x,]
                            wilcox_test <- wilcox.test(data[data[,1] == sn$Var1, "Value"],
                                                       data[data[,1] == sn$Var2, "Value"]) 
                            data.frame(Test = "Wilcoxon",
                                       Condition = paste(sn$Var1, " vs ", sn$Var2), 
                                       pValue = wilcox_test$p.value,
                                       conf.int = paste(wilcox_test$conf.int, collapse = ";")
                            )
                          })
      )
    } else if(length(vars) > 2){
      steps <- c(steps, paste("Step ", step_counter, ". groups are more than 2, I will use Kruskal-Wallis for analysis", "\n"))
      step_counter <- step_counter + 1  
      colnames(data) <- c("SampleName", "Value")
      data$SampleName <- as.factor(data$SampleName)
      
      kruskal_test <- kruskal.test(Value ~ SampleName, data = data)
      
      resKRUSKAL <- data.frame(Test = "Kruskal-Wallis",
                             Condition = paste(kruskal_test$call)[2], 
                             pValue = kruskal_test$p.value,
                             conf.int = paste("-", collapse = ";"))
      
      if (resKRUSKAL$pValue < 0.05) {
        steps <- c(steps, paste("Step ", step_counter, ". Kruskal-Wallis p-value <", resKRUSKAL$pValue, ", performing pairwise Wilcoxon tests", "\n"))
        step_counter <- step_counter + 1  
        
        combo <- combn(vars, 2)
        combo <- data.frame(Var1 = combo[1,], Var2 = combo[2,])
        
        combo <- combo[combo$Var1 != combo$Var2, ]
        resPairwise <- do.call(rbind,
                               lapply(1:dim(combo)[1], function(x){
                                 sn <- combo[x,]
                                 wilcox_test <- wilcox.test(data[data[,1] == sn$Var1, "Value"],
                                                            data[data[,1] == sn$Var2, "Value"]) 
                                 data.frame(Test = "Wilcoxon",
                                            Condition = paste(sn$Var1, " vs ", sn$Var2), 
                                            pValue = wilcox_test$p.value,
                                            conf.int = paste(wilcox_test$conf.int, collapse = ";")
                                 )
                               })
        )
      } else {
        steps <- c(steps, paste("Step ", step_counter, ". Kruskal-Wallis p-value >=", resKRUSKAL$pValue, ", no pairwise Wilcoxon tests performed", "\n"))
      }
    }
    
    return(list(resTTest = resTTest, test = resKRUSKAL, pairwise = resPairwise, steps = steps))
  }
}