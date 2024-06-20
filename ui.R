library(shinydashboard)
library(shinyjs)
library(shinydisconnect)
library(shinyBS)
library(tippy)
library(shinythemes)

ui <- navbarPage(
  "PlotaPep",
  theme = shinytheme("flatly"),
  useShinyjs(),
  
  tabPanel(
    "Fasta Input",
    mainPanel(fluidRow(box(
      fileInput(
        "fastaInput",
        "Choose FASTA file : ",
        multiple = FALSE,
        accept = c("fasta", ".FASTA", "FASTA", ".fasta"),
        # close c
        buttonLabel = "Browse...",
        placeholder = "No FASTA Selected"
      )
      
    )),
    fluidRow(
      bsTooltip(
        id = "fastaInput",
        title = "This is an input",
        placement = "right",
        trigger = "manual"
      ),
    )),
    fluidRow(column(4, wellPanel(
      p("How to use:"),
      p(
        " 1. Choose FASTA File. Can be any FASTA file but if it doesn't have the proteins you're looking for it will not be plotted"
      ),
      p(
        " 2. Input data set(s) to be plotted can be one or multiple. If multiple dataset are being plotted then be sure to choose the comparative option"
      ),
      p(
        " 3. Select file output type. Currently supported MSFragger combined peptide with and without modifications and DIANN pr.matrix."
      ),
      p(
        " 4. Set cohorts to be plotted, the cohort set for the control is what fold change will be based off of."
      ),
      p(" 5. Set other options to preference then run."),
      p(" 6. Use the text input above the plots to search for specific plots.")
    )))
  ),
  
  
  
  tabPanel("File Input", (mainPanel(
    fluidRow(
      box(
        fileInput(
          "fileInput",
          "Choose file : ",
          multiple = TRUE,
          accept = c("csv",
                     "comma-separated-values",
                     ".csv",
                     ".tsv"),
          buttonLabel = "Browse...",
          placeholder = "No Files Selected"
        ),
        radioButtons(
          "file_typeInput",
          "Select Search Output Type : ",
          c(
            "MSFragger No modifications " = 0,
            "MSFragger Modified/Ion" = 1,
            "DIA-NN" = 2,
            "Generic Dataframe" = 3
          )
        ),
        radioButtons(
          "intensityInput",
          "Select Intensity : ",
          c(
            "Intensity" = 0,
            "Total" = 1,
            "Unique" = 2
          )
        ),
        radioButtons("lfqInput", "Select LFQ Type : ", c("None" = 1, "Max" = 0)),
        numericInput(
          "cutoffInput",
          "Select Cutoff : ",
          0,
          min = 0,
          max = 100
        ),
        textInput(
          "ctrInput",
          "Input control identifier : ",
          placeholder = "Example CTR, AD, DLB etc."
        ),
        textInput(
          "cohortInput",
          "Input cohort identifiers : ",
          
          placeholder = "Example AD, DLB etc."
        ),
        radioButtons(
          "parametricInput",
          "Select Parametric or Nonparametric statistical test : ",
          c("Parametric" = 1, "Nonparametric" = 0)
        ),
        radioButtons(
          "comparativeInput",
          "Select Comparative TTest Results (plots multiple dataset ttest results on single plot) : ",
          c("Noncomparative" = 0, "Comparative" = 1)
        ),
        numericInput(
          "pInput",
          "Select p value for significance : ",
          0.05,
          min = 0,
          max = 100
        ),
        textInput(
          "colorsInput",
          "Input hex value color codes :",
          placeholder = "Example CTR, AD, DLB etc.",
          value = "#E495A5,#86B875,#7DB0DD,#FFA500"
        ),
        numericInput(
          "stdInput",
          "Select standard deviations : ",
          2,
          min = 2,
          max = 6
        ),
        actionButton("confirmFile", "Run Plotter"),
      )
      
    ),
    fluidRow(
      bsTooltip(
        id = "cohortInput",
        title = "Input cohort identifier from your file input. These should the renamed combined_peptide.tsv from msfragger.",
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
        id = "intensityInput",
        title = "This is an Fragger option only",
        placement = "right",
        trigger = "hover"
      ),
      bsTooltip(
        id = "lfqInput",
        title = "This is an Fragger option only",
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
        id = "ctrInput",
        title = "Input control identifier from your file input. This will determine how the fold change is plotted",
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
      )
    )
  )),
  fluidRow(column(
    4, wellPanel(
      p("How to use:"),
      p(
        " 1. Choose FASTA File. Can be any FASTA file but if it doesn't have the proteins you're looking for it will not be plotted"
      ),
      p(
        " 2. Input data set(s) to be plotted can be one or multiple. If multiple dataset are being plotted then be sure to choose the comparative option"
      ),
      p(
        " 3. Select file output type. Currently supported MSFragger combined peptide with and without modifications and DIANN pr.matrix."
      ),
      p(
        " 4. Set cohorts to be plotted, the cohort set for the control is what fold change will be based off of."
      ),
      p(" 5. Set other options to preference then run."),
      p(" 6. Use the text input above the plots to search for specific plots.")
    )
  ))),
  
  tabPanel(
    "Plot Output",
    
    mainPanel(fluidRow(
      box(
        disabled(
          downloadButton('downloadPlot', "Download Plots and Final Data Table")
        ),
        disabled(actionButton("previousPlot", "Previous")),
        disabled(actionButton("nextPlot", "Next")),
        radioButtons(
          "svgInput",
          "Select .SVG or .PNG file formats : ",
          c("PNG" = 0, "SVG" = 1)
        ),
        tagList(
          textInput(
            "uniprotInput",
            "Input Uniprot IDs separated by commas (,) : ",
            placeholder = "Example P05067, P02649, etc."
          ),
          value = NULL
        ),
        
      )
    )),
    fluidRow(
      bsTooltip(
        id = "uniprotInput",
        title = "Input uniprot ids to plot only selected proteins.",
        placement = "right",
        trigger = "hover"
      ),
    ),
    fluidRow(uiOutput("plot"))
  ),
)
