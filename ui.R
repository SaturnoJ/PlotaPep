library(shinydashboard)
library(shinyjs)
library(shinydisconnect)
library(shinyBS)
library(tippy)

ui <- dashboardPage(
  dashboardHeader(title = "PlotaPep"),
  
  dashboardSidebar(),
  dashboardBody(
    useShinyjs(),
    
    fluidRow(
      box(
        uiOutput("fastaInput"),
        textOutput("test"),
        
        uiOutput("fileInput"),
        uiOutput("file_typeInput"),
        uiOutput("intensityInput"),
        uiOutput("lfqInput"),
        uiOutput("cutoffInput"),
        uiOutput("uniprotInput"),
        uiOutput("ctrInput"),
        uiOutput("cohortInput"),
        uiOutput("parametricInput"),
        uiOutput("comparativeInput"),
        # uiOutput("yInput"),
        uiOutput("pInput"),
        uiOutput("colorsInput"),
        uiOutput("stdInput"),
        uiOutput("svgInput"),
        
        disabled(actionButton("confirmFasta", "Confirm FASTA Input")),
        
        uiOutput("confirmFile"),
        disabled(downloadButton('downloadPlot', "Download the table")),
        disabled(actionButton("clickThrough", "Next"))
      ),
      bsTooltip(
        id = "cohortInput",
        title = "Input cohort identifier from your file input. These should the renamed combined_peptide.tsv from msfragger.",
        placement = "right",
        trigger = "hover"
      ),
      bsTooltip(
        id = "fastaInput",
        title = "This is an input",
        placement = "right",
        trigger = "hover"
      ),
      bsTooltip(
        id = "colorsInput",
        title = "These can be common names for colors such as red or blue but for special color values use the respective hex value such as #FF0000 or #0000FF.",
      placement = "right",
      trigger = "hover"
    ),
    bsTooltip(
      id = "fileInput",
      title = "One or multiple combined_peptide.tsv that have been renamed with the cohort that is to be displayed on the plots.",
      placement = "right",
      trigger = "hover"
    ),
    bsTooltip(
      id = "intensityInput",
      title = "This is an input",
      placement = "right",
      trigger = "hover"
    ),
    bsTooltip(
      id = "lfqInput",
      title = "This is an input",
      placement = "right",
      trigger = "hover"
    ),
    bsTooltip(
      id = "cutoffInput",
      title = "This is an input",
      placement = "right",
      trigger = "hover"
    ),
    bsTooltip(
      id = "uniprotInput",
      title = "Input uniprot ids to plot only selected proteins.",
      placement = "right",
      trigger = "hover"
    ),
    bsTooltip(
      id = "ctrInput",
      title = "Input control identifier from your file input. These should the sample file names from the acquisition.",
      placement = "right",
      trigger = "hover"
    ),
    bsTooltip(
      id = "parametricInput",
      title = "Choose if your statistical test is parametric or nonparametric",
      placement = "right",
      trigger = "hover"
    ),
    bsTooltip(
      id = "comparativeInput",
      title = "Takes two or more cohorts performs independent ttest on each and plot them to the same plot for direct comparison.",
      placement = "right",
      trigger = "hover"
    ),
    # bsTooltip(
    #   id = "yInput",
    #   title = "Height of the y-axis for the plots",
    #   placement = "right",
    #   trigger = "hover"
    # ),
    bsTooltip(
      id = "pInput",
      title = "Choose max p-value for statisical significance",
      placement = "right",
      trigger = "hover"
    ),
    bsTooltip(
      id = "stdInput",
      title = "Choose how many standard deviations from the mean are used for outlier removal. Max: 6 Min: 2",
      placement = "right",
      trigger = "hover"
    ),
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
