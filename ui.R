library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(title = "Tally"),
  
  dashboardSidebar(),
  dashboardBody(useShinyjs(),
                
                fluidRow(box(
                  uiOutput("fastaInput"),
                  textOutput("test"),
                  
                  uiOutput("fileInput"),
                  uiOutput("intensityInput"),
                  uiOutput("lfqInput"),
                  uiOutput("cutoffInput"),
                  uiOutput("fileTypeInput"),
                  disabled(actionButton("confirmFasta", "Confirm FASTA Input")),
                 
                  uiOutput("uniprotInput"),
                  uiOutput("confirmFile"),
                  disabled(downloadButton('download',"Download the table"))
                  
                ))
                , fluidRow(dataTableOutput("contents")))
)
