# 
# Appui <- system.file("Shiny","ui.R", package = "ORCA")
# Appserver <- system.file("Shiny","server.R", package = "ORCA")
# 
source("inst/Shiny/server.R")
source("inst/Shiny/ui.R")

shinyApp(ui, server,
         options =  options(shiny.maxRequestSize=1000*1024^2,
                            shiny.launch.browser = .rs.invokeShinyWindowExternal)
         )


