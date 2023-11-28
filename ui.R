library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(title = "Tally"),
  
  dashboardSidebar(),
  dashboardBody(useShinyjs(),
                
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
                    
                    disabled(actionButton("confirmFasta", "Confirm FASTA Input")),
                    
                    uiOutput("confirmFile"),
                    disabled(downloadButton('downloadPlot', "Download the table")),
                    disabled(actionButton("clickThrough", "Next"))
                  )
                )
                , fluidRow(uiOutput("plot")))
)
