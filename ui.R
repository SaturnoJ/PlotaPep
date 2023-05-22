library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(title = "Tally"),
  
  dashboardSidebar(),
  dashboardBody(useShinyjs(),
                
                fluidRow(box(
                  uiOutput("fastaInput"),
                  
                  
                  uiOutput("fileInput"),
                  disabled(actionButton("confirmFasta", "Confirm FASTA Input")),
                 
                  uiOutput("uniprot"),
                  uiOutput("confirmFile"),
                  
                ))
                , fluidRow(dataTableOutput("contents")))
)
