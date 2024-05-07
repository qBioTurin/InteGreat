library(jsonlite)
library(shiny)
library(shinybusy)
options(shiny.maxRequestSize= 3*1024^4) # more Size for videos

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
                                       videoResult = NULL,
                                       images=NULL,
                                       list_b_n=NULL
                                       )
  
  DataIntegrationModule <- reactiveValues(dataLoaded = NULL,
                                          data = NULL,
                                          wbTabs = NULL, 
                                          pcrTabs = NULL,
                                          cytotoxTabs= NULL,
                                          endocTabs=NULL,
                                          otherTabs = NULL,
                                          otherTabsMean = NULL)
  
  MapAnalysisNames =c("WB", "WB comparison", "Endocytosis", "ELISA", "RT-qPCR", "Cytotoxicity", "Border") 
  names(MapAnalysisNames) =c("wbResult", "wbquantResult", "endocResult", "elisaResult", "pcrResult", "cytotoxResult","borderResult") 
  
  
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
  
  
  # Define the reactive expression to call the Docker container
  callDocker_start <- function(video_path) {
    cartella_start= tempdir()
    file.copy(video_path,paste0(cartella_start,"/movies.lif")) #rinomino guardo plotoutput in aux function e render plot raster plot 
    #file.copy(video_path,paste0("/Users/lorenzopace/Documents/Universita/Stage/ORCA/inst/Shiny/prova","/movies.lif"))
    share_command=paste0(cartella_start,":/home")
    docker_command <- paste("docker run -d -v", share_command, "lorenzo/border python3 script_docker.py first_frame ../home/movies.lif")
    #docker_command<-paste("docker run -d -v /Users/lorenzopace/Documents/Universita/Stage/ORCA/inst/Shiny/prova:/home lorenzo/border python3 script_docker.py first_frame ../home/movies.lif")
    
    system(docker_command, intern = TRUE)
    
    dockerFinished <- function() {
      docker_ps <- system("docker ps -q", intern = TRUE)
      if (length(docker_ps) == 0) {
        return(TRUE)  # No running containers
      } else {
        return(FALSE)  # There are running containers
      }
    }

    # Wait until either result.json is created or Docker process is finished
    while (!dockerFinished()) {
      Sys.sleep(1)  # Wait for 1 second before checking again
    }
    
    # Read the result from the file
    result_file <- paste0(cartella_start,"/resultFirstFrames.json")
    result <- read_json(result_file)
    unlink(cartella_start)
    #result <- read_json("/Users/lorenzopace/Documents/Universita/Stage/ORCA/inst/Shiny/prova/resultFirstFrames.json")
    return(result)
  }
  
  
    # Observer for uploading video and calling Docker
    observeEvent(input$uploadVideo_button, {
      show_spinner()
      video_path <- input$uploaded_video$datapath
      images <- callDocker_start(video_path)
      num_images <- length(images)
      # Convert from list to matrix for plotting
      for (i in 1:num_images) {
        images[[i]] <- matrix(unlist(images[[i]]), nrow = length(images[[i]]),ncol=length(images[[i]]), byrow = TRUE)  
        # Scale pixel values to range [0, 1]
        images[[i]] <- images[[i]] / 255
      }
      
      # Initialize a reactiveValues object to store checkbox states
      checkbox_states <- reactiveValues()
      for (i in 1:num_images) {
        checkbox_states[[paste0("image_", i)]] <- FALSE
      }
      # Current image index (starts at 1)
      current_image <- reactiveVal(1)
      
    
      output$imageGallery <- renderUI({
        # Only render the current image plot
        image_ui <- column(
          12,
          align = "center",
          
          # Add horizontal spacing using spacer tags and margins
          fluidRow(
            column(6,
                   plotOutput(paste0("plotImage_", current_image())),
                   # Add margin-right for the plot
                   style = "margin-right: 20px;"  # Adjust as needed
            ),
            # Spacer tag for horizontal space
            #tags$span(style = "flex: 1;"),  # Adjust flex value for desired space
            column(6,
                   # Add margin-left for the checkbox
                   style = "margin-top: 120px;",
                   checkboxInput(inputId = paste0("checkbox_", current_image()), label = paste("Select Image", current_image()), value = FALSE)
            )
          ),
          actionButton("nextButton", "Next Video"),
          actionButton("start", "Start Analysis"),
        )
        image_ui
        #hide_spinner() 
      })
      
      # Define a function to render each plot
      renderImagePlot <- function(i) {
        renderPlot({
          print(paste("Image", i, "dimensions:", dim(images[[i]])))
          plot(c(1, dim(images[[i]])[2]), c(1, dim(images[[i]])[1]), type = 'n', ann = FALSE)
          rasterImage(images[[i]], 1, 1, dim(images[[i]])[2], dim(images[[i]])[1])
        }, width = dim(images[[i]])[2], height = dim(images[[i]])[1])
      }
      
      
        observeEvent(input[[paste0("checkbox_", current_image())]], {
          checkbox_states[[paste0("image_", current_image())]] <- input[[paste0("checkbox_", current_image())]]
        })
  
      
      output$plotImage_1 <- renderImagePlot(1)  # Initially render the first image
      hide_spinner()
      
      observeEvent(input$nextButton, {
        if(current_image()<num_images){
        current_image(current_image() + 1)
        if(current_image()==num_images){
          shinyjs::hide("nextButton")
        }
        # Update rendered plot with the new image
        output[[paste0("plotImage_", current_image())]] <- renderImagePlot(current_image())
        }
      })
      
      observeEvent(input$start, {
        show_spinner()
          border_detection(checkbox_states,video_path)
      })
  })
    


 callDocker_border <- function(checkbox_states,video_path) {
   cartella_border= tempdir()
   
   file.copy(video_path,paste0(cartella_border,"/movies.lif")) #rinomino guardo plotoutput in aux function e render plot raster plot 
   
   checkbox_states <- reactiveValuesToList(checkbox_states)  
   # Convert boolean values to strings (false -> "false", true -> "true")
   checkbox_states_str <- sapply(checkbox_states, function(x) toString(x))
   checkbox_states_str <- paste(checkbox_states_str, collapse = " ")
   checkbox_states <- unlist(strsplit(checkbox_states_str, " "))
   writeLines(checkbox_states, paste0(cartella_border,"/checkbox_states.txt"))
   #writeLines(checkbox_states, paste0("/Users/lorenzopace/Documents/Universita/Stage/ORCA/inst/Shiny/prova","/checkbox_states.txt"))
   share_command=paste0(cartella_border,":/home")
   docker_command <- paste("docker run -d -v", share_command, "lorenzo/border python3 script_docker.py border_detection ../home/movies.lif ../home/checkbox_states.txt")
   #docker_command <- paste("docker run -d -v /Users/lorenzopace/Documents/Universita/Stage/ORCA/inst/Shiny/prova:/home lorenzo/border python3 script_docker.py border_detection ../home/movies.lif ../home/checkbox_states.txt")
   # Run the Docker command and capture its output

   # Run the Docker container
   system(docker_command,intern = TRUE)
   
   dockerFinished <- function() {
     docker_ps <- system("docker ps -q", intern = TRUE)
     if (length(docker_ps) == 0) {
       return(TRUE)  # No running containers
     } else {
       return(FALSE)  # There are running containers
     }
   }

   while (!dockerFinished()) {
     Sys.sleep(2)  # Wait for 2 seconds before checking again
   }
   
    result_border <- paste0(cartella_border,"/resultBorder.json")
    result_nucleus <-paste0(cartella_border,"/resultNucleus.json")
   # 
   # # # Read the result from the file
    result_b <- read_json(result_border)
    result_n <- read_json(result_nucleus)
   #result_b <- read_json("/Users/lorenzopace/Documents/Universita/Stage/ORCA/inst/Shiny/prova/resultBorder.json")
   #result_n <- read_json("/Users/lorenzopace/Documents/Universita/Stage/ORCA/inst/Shiny/prova/resultNucleus.json")
   
   unlink(cartella_border)
   return(list(border=result_b,nucleus=result_n))
 }

 border_detection<- function(checkbox_states,video_path){
   #print(checkbox_states[[paste0("image_", 1)]])
  list_b_n <- callDocker_border(checkbox_states,video_path)
  num_images=length(list_b_n$border)
   images <- vector("list", length = num_images)
   for (i in 1:num_images) {
     images[[i]] <- matrix(unlist(list_b_n$border[[i]]) + unlist(list_b_n$nucleus[[i]]) ,nrow = length(list_b_n$border[[i]]),ncol=length(list_b_n$border[[i]]), byrow = TRUE)  
   }
   
   for (i in 1:num_images) {
     list_b_n$border[[i]] <- matrix(unlist(list_b_n$border[[i]]),nrow = length(list_b_n$border[[i]]),ncol=length(list_b_n$border[[i]]), byrow = TRUE)
   }
   for (i in 1:num_images) {
     list_b_n$nucleus[[i]] <- matrix(unlist(list_b_n$nucleus[[i]]),nrow = length(list_b_n$nucleus[[i]]),ncol=length(list_b_n$nucleus[[i]]), byrow = TRUE)
   }
   DataAnalysisModule$images<-images
   DataAnalysisModule$list_b_n<-list_b_n
   # # Convert from list to matrix for plotting
   # for (i in 1:num_images) {
   #   images[[i]] <- matrix(unlist(images[[i]]), nrow = length(images[[i]]),ncol=length(images[[i]]), byrow = TRUE)  
   #   # Scale pixel values to range [0, 1]
   #   #images[[i]] <- images[[i]] / 255
   # }
   
   # Current image index (starts at 1)
   current_image <- reactiveVal(1)
   
   # Reactive values to store click coordinates over the plot
   click_coords <- reactiveValues(click1 = NULL, click2 = NULL)
   
   output$imageGallery <- renderUI({
     # Only render the current image plot
     image_ui <- column(
       12,
       align = "center",
       
       fluidRow(
         column(6,
                plotOutput(paste0("plotImage_", current_image()),click = "plot_click"),
                style = "margin-right: 20px;" ,
                tags$script(HTML(sprintf("
                 $(document).ready(function() {
                   var plotID = '%s'; // Get the plot ID dynamically
                   
                   // Function to handle click events
                   function handleClick(e) {
                     var parentOffset = $('#' + plotID).offset(); 
                     var relX = e.pageX - parentOffset.left;
                     var relY = e.pageY - parentOffset.top;
                     Shiny.setInputValue('click_coords', {x: relX, y: relY});
                     
                   }
                   
                   // Bind click event to the plot
                   $('#' + plotID).click(handleClick);
                   
                   // Listen for custom message captureClick
                   Shiny.addCustomMessageHandler('captureClick', function(message) {
                     // Execute the handleClick function when the custom message is received
                     $('#' + plotID).click(handleClick);
                   });
                 });
               ", paste0("plotImage_", current_image()))))
         ),
         column(6,
                # Add margin-left for the checkbox
                style = "margin-top: 120px;",
         )
       ),
       
       actionButton("prevAnalysis", "Previous Video"),
       actionButton("nextAnalysis", "Next Video"),
       actionButton("startTracking","Start tracking"),
       use_busy_spinner(spin = "fading-circle")
     )
     shinyjs::hide("prevAnalysis")
     image_ui
   })
   
   observe({
     session$sendCustomMessage("captureClick", NULL)
   })
   
   observeEvent(input$plot_click, {
     click <- input$plot_click
     if (is.null(click_coords$click1)) {
       click_coords$click1 <- click
     } else if (is.null(click_coords$click2)) {
       click_coords$click2 <- click
       addLine(click_coords$click1$x,click_coords$click1$y,click_coords$click2$x,click_coords$click2$y)
       click_coords$click1 =NULL
       click_coords$click2 =NULL
     }
   })
   
   addLine<-function(x1,y1,x2,y2){
     i=current_image()
    images<-DataAnalysisModule$images
     y_length=dim(images[[i]])[1] 
     x1 <- round(x1)
     x2 <- round(x2)
     y1 <- y_length-round(y1) # y coordinates by click are opposite than matrix
     y2 <- y_length-round(y2)
     #matrix<-images[[i]]
     #border<-list_b_n$border[[i]]
     if (x1 == x2) {
       # Vertical line
       for (y in floor(min(y1, y2)):floor(max(y1, y2))) {
         images[[i]][y, x1] <- 1  # Set the value at each coordinate to 1
         list_b_n$border[[i]][y, x1] <- 1
       }
     } else if (y1 == y2) {
       # Horizontal line
       for (x in floor(min(x1, x2)):floor(max(x1, x2))) {
         images[[i]][y1, x] <- 1  # Set the value at each coordinate to 1 
         list_b_n$border[[i]][y1, x] <- 1
         
       }
     } else {
       # Diagonal line (using Bresenham's line algorithm)
       dx <- abs(x2 - x1)
       dy <- abs(y2 - y1)
       sx <- ifelse(x1 < x2, 1, -1)
       sy <- ifelse(y1 < y2, 1, -1)
       err <- dx - dy
       
       while (x1 != x2 || y1 != y2) {
         images[[i]][y1, x1] <- 1
         list_b_n$border[[i]][y1, x1] <- 1
         e2 <- 2 * err
         if (e2 > -dy) {
           err <- err - dy
           x1 <- x1 + sx
         }
         if (e2 < dx) {
           err <- err + dx
           y1 <- y1 + sy
         }
       }
     }
     DataAnalysisModule$images[[i]]<-images[[i]]
     DataAnalysisModule$list_b_n$border<-list_b_n
     output[[paste0("plotImage_", i)]] <- renderImagePlot(images[[i]])
   }
   
   # renderImagePlot <- function(i) {
   #   renderPlot({
   #     print(paste("Image", i, "dimensions:", dim(images[[i]])))
   #     plot(c(1, dim(images[[i]])[2]), c(1, dim(images[[i]])[1]), type = 'n', ann = FALSE)
   #     rasterImage(images[[i]], 1, 1, dim(images[[i]])[2], dim(images[[i]])[1])
   #   }, width = dim(images[[i]])[2], height = dim(images[[i]])[1])
   # 
   # }
   
   renderImagePlot <- function(image) {
     renderPlot({
       print(paste("Image", i, "dimensions:", dim(image)))
       plot(c(1, dim(image)[2]), c(1, dim(image)[1]), type = 'n', ann = FALSE)
       rasterImage(image, 1, 1, dim(image)[2], dim(image)[1])
     }, width = dim(image)[2], height = dim(image)[1])

   }

   
   output$plotImage_1 <- renderImagePlot(images[[1]])  # Initially render the first image
   hide_spinner()
   
   observeEvent(input$nextAnalysis, {
     if(current_image()<num_images){
      current_image(current_image() + 1)
      if(current_image()==num_images){
        shinyjs::hide("nextAnalysis")
      }
      if(current_image()==2)
        shinyjs::show("prevAnalysis")
     
      # Update rendered plot with the new image
      output[[paste0("plotImage_", current_image())]] <- renderImagePlot(images[[current_image()]])
     }
   })
   
   observeEvent(input$prevAnalysis, {
     if(current_image()>1){
       current_image(current_image() - 1)
      if(current_image()==1)
        shinyjs::hide("prevAnalysis")
     if(current_image()==num_images-1)
       shinyjs::show("nextAnalysis")
       # Update rendered plot with the new image
       output[[paste0("plotImage_", current_image())]] <- renderImagePlot(images[[current_image()]])
     }
   })
   
   observeEvent(input$startTracking, {
     show_spinner()
     DataAnalysisModule$list_b_n<-list_b_n
     tracking(checkbox_states,video_path,list_b_n$border,list_b_n$nucleus)
   })
  
 }
 
 
}
 