library(shinydashboard)
library(shinyjs)
library(shinydisconnect)

ui <- dashboardPage(
  dashboardHeader(title = "Tally"),
  
  dashboardSidebar(),
  dashboardBody(
    useShinyjs(),
    
    fluidRow(
      box(
        uiOutput("fastaInput"),
        textOutput("test"),
        
        uiOutput("fileInput"),
        uiOutput("intensityInput"),
        uiOutput("lfqInput"),
        uiOutput("cutoffInput"),
        uiOutput("fileTypeInput"),
        uiOutput("uniprotInput"),
        uiOutput("ctrInput"),
        uiOutput("cohortInput"),
        uiOutput("parametricInput"),
        uiOutput("yInput"),
        uiOutput("pInput"),
        
        disabled(actionButton("confirmFasta", "Confirm FASTA Input")),
        
        uiOutput("confirmFile"),
        disabled(downloadButton('downloadPlot', "Download the table")),
        disabled(actionButton("clickThrough", "Next"))
      )
    )
    ,
    fluidRow(uiOutput("plot")),
    disconnectMessage(
      text = "An error occurred. For parametric analysis please ensure that there is at least one protein with at least two samples across all cohorts. Also ensure that the cohorts listed match the input file names.",
      refresh = "Refresh",
      background = "#FFFFFF",
      colour = "#444444",
      refreshColour = "#337AB7",
      overlayColour = "#000000",
      overlayOpacity = 0.6,
      width = 450,
      top = 50,
      size = 22,
      css = ""
    )
    
  )
)
