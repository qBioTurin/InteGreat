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

resetPanel <- function(type, flags, panelStructures, numberOfPlanes, planeSelected, result, output, panelData, pcrResult = NULL) {
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
           pcrResult$Initdata <- NULL
           pcrResult$selectPCRcolumns <- NULL
           pcrResult$data <- NULL
           pcrResult$PCRnorm <- NULL
           pcrResult$BaselineExp <- NULL
           pcrResult$plotPRC <- NULL
           pcrResult$NewPCR <- NULL
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
        return (list(message = "Please, select an Excel file",call = ""))
      } else if(tolower(tools::file_ext(filename)) != "xls" && tolower(tools::file_ext(filename)) != "xlsx") {
        return(list(message = "Please, upload a file with a .xls or .xlsx extension.", call = ""))
      } 
      
      x = readxl::read_excel(filename, col_names = colname)
      if (allDouble) {
        xstr = which(sapply(1:dim(x)[2], function(i) !is.double(x[[i]])))
        if (length(xstr) > 0) {
          for (i in xstr) {
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
            df = data.frame(row = x$rows, col = x$cols,
                            fill = ifelse(length(x$style$fill$fillFg) > 0, 
                                          color,
                                          "white"))
          }
        })
        l = do.call(rbind, l[lengths(l) > 0])
        SN = table(l$fill)
        l$SN = paste0("Color ", match(l$fill, names(SN)))
        
        tb.SN = matrix("", nrow = max(l$row), ncol = max(l$col))
        
        for (j in 1:ncol(tb.SN)) {
          fill_col = l %>% filter(col == j)
          tb.SN[fill_col$row, j] = fill_col$SN
        }
        
        col = l %>% select(fill, SN) %>% distinct()
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
loadImage = function(pathImage){
  im <- OpenImageR::readImage(pathImage,as.is = T,convert=TRUE)
  
  if( length(dim(im)) != 2  ) 
    im = rgb_2gray(im)
  
  JpegImage = tempfile(fileext = "jpeg")
  
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
  
  if(analysis =="WB"){
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
  } else if(analysis =="WB comparison"){
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
  }
  else if(analysis =="RT-qPCR"){
    wb <- createWorkbook("RTqPCR")
    
    addWorksheet(wb,"Table")
    writeDataTable(wb, sheet = "Table", ResultList[["Initdata"]])
    
    addWorksheet(wb,"Norm PRC")
    writeDataTable(wb,ResultList[["NewPCR"]], sheet="Norm PRC")
    
    print(ResultList[["plotPRC"]])
    insertPlot(wb = wb,  sheet="Norm PRC",
               startCol=dim(ResultList[["NewPCR"]])[2]+ 2)
    
  }
  
  saveWorkbook(wb, filename)
  return(1)
}



