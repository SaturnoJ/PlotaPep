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
                  disabled(actionButton("confirmFasta", "Confirm FASTA Input")),
                 
                  uiOutput("uniprot"),
                  uiOutput("confirmFile"),
                  disabled(downloadButton('download',"Download the table"))
                  
                ))
                , fluidRow(dataTableOutput("contents")))
)
