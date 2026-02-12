library(shinydashboard)
library(shinyjs)
library(shinydisconnect)
library(shinyBS)
library(tippy)
library(shinythemes)
library(DT)

ui <- navbarPage(
  "PlotaPep",
  theme = shinytheme("flatly"),
  useShinyjs(),
########
  tabPanel(
    "Key Input",
    mainPanel(fluidRow(box(
      fileInput(
        "originalInput",
        "Input output file from database search tools to be renamed : ",
        multiple = FALSE,
        accept =c("csv",
                  "comma-separated-values",
                  ".csv",
                  ".tsv",".txt"),
        # close c
        buttonLabel = "Browse...",
        placeholder = "No Original Selected"
      ),

      downloadButton("confirmKey", "Run Key"),
      DTOutput('tableInput')
    ))),
    fluidRow(column(4, wellPanel(
      p("How to use:"),
      p(
        " 1. Choose file that needs its sample names unobscured. This can be any output file such as report.pr_matrix, combined_peptide, combined_modified, etc."
      ),
      p(
        " 2. Put new names into table below. Press CTRL-Enter to save the values."
      ),
      p(
        " 3. Run then download"
      )
    )))
  ),
#########

#######
  tabPanel("Woods Input", (mainPanel(
    fluidRow(
     shinydashboard::box(
        fileInput(
          "fastaInput",
          "Choose FASTA file : ",
          multiple = FALSE,
          accept = c(".fasta",".FASTA"),
          buttonLabel = "Browse...",
          placeholder = "No FASTA selected"
        )
        ,
        fileInput(
          "fileInput",
          "Choose file : ",
          multiple = TRUE,
          accept = c("csv",
                     "comma-separated-values",
                     ".csv",
                     ".tsv",".txt"),
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
          "ctr_cutoffInput",
          "Select Control Cutoff : ",
          3,
          min = 2,
          max = 100
        ),
        numericInput(
          "case_cutoffInput",
          "Select Case Cutoff : ",
          3,
          min = 2,
          max = 100
        ),
        numericInput(
          "peptide_minInput",
          "Select Min Peptide Amount for Plots : ",
          3,
          min = 1,
          max = 100
        ),
        textInput("ctrInput",
                  "Input Control Identifier : ",
                  placeholder = "Example CTR, AD, DLB etc."),
        textInput("caseInput",
                  "Input Case Identifiers : ",

                  placeholder = "Example AD, DLB etc."),
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
          4,
          min = 2,
          max = 6
        ),
        actionButton("confirmFile", "Run Plotter"),
      )

    ),
    fluidRow(
      bsTooltip(
        id = "caseInput",
        title = "Input case identifier from your file input. These are sample names from the file. Be specific i.e. Old_AD or Young_AD if multiple cases have similar names.",
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
        id = "file_typeInput",
        title = "This can be any output file such as report.pr_matrix, combined_peptide, combined_modified, etc.",
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
        id = "ctr_cutoffInput",
        title = "Amount of control values needed for a protein to be used for statistical test. Min = 2",
        placement = "right",
        trigger = "hover"
      ),  bsTooltip(
        id = "case_cutoffInput",
        title = "Amount of case values needed for a protein to be used for statistical test. Min = 2",
        placement = "right",
        trigger = "hover"
      ),
      bsTooltip(
        id = "peptide_minInput",
        title = "Set min value for peptides to be plotted. If the amount of peptides is less than specified number they will not be plotted.",
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
        title = "Takes two or more cases performs independent ttest on each and plot them to the same plot for direct comparison.",
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
        " 4. Set cases to be plotted, the case set for the control is what fold change will be based off of."
      ),
      p(" 5. Set other options to preference then run."),
      p(" 6. Use the text input above the plots to search for specific plots.")
    )
  ))),
############
  tabPanel(
    "Plot Output",

    mainPanel(fluidRow(
     shinydashboard::box(
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

  tabPanel(
    "PCA Output",

    mainPanel(

    fluidRow(
     shinydashboard::box(

        fileInput(
          "fileInput",
          "Choose file : ",
          multiple = TRUE,
          accept = c("csv",
                     "comma-separated-values",
                     ".csv",
                     ".tsv",".txt"),
          buttonLabel = "Browse...",
          placeholder = "No Files Selected"
        ),



        textInput("conditionAInput",
                  "Input Conditions : ",
                  placeholder = "Example CTR, AD, DLB, F_12, etc."),


        textInput("labelAInput",
                  "Input Labels : ",
                  placeholder = "Example CTR, AD, DLB, F_12, etc."),




        actionButton("confirmPCA", "Run Plotter"),

      ),


    fluidRow(uiOutput("pcaPlot"))
  ))),


tabPanel(
  "Volcano Output",

  mainPanel(

    fluidRow(
     shinydashboard::box(

        fileInput(
          "fileInput",
          "Choose file : ",
          multiple = TRUE,
          accept = c("csv",
                     "comma-separated-values",
                     ".csv",
                     ".tsv",".txt"),
          buttonLabel = "Browse...",
          placeholder = "No Files Selected"
        ),
        # radioButtons(
        #   "file_typeInput",
        #   "Select Search Output Type : ",
        #   c(
        #     "MSFragger No modifications " = 0,
        #     "MSFragger Modified/Ion" = 1,
        #     "DIA-NN" = 2,
        #     "Generic Dataframe" = 3
        #   )
        # ),
        # radioButtons(
        #   "intensityInput",
        #   "Select Intensity : ",
        #   c(
        #     "Intensity" = 0,
        #     "Total" = 1,
        #     "Unique" = 2
        #   )
        # ),


        textInput("conditionAInput",
                  "Input A Conditions : ",
                  placeholder = "Example CTR, AD, DLB, F_12, etc."),
        textInput("conditionBInput",
                  "Input B Condition : ",

                  placeholder = "Example AD, DLB, F_9, etc."),
        textInput("labelAInput",
                  "Input A Label : ",
                  placeholder = "Example CTR, AD, DLB, F_12, etc."),
        textInput("labelBInput",
                  "Input B Label : ",

                  placeholder = "Example AD, DLB, F_9, etc."),




        actionButton("confirmVol", "Run Plotter"),

      ),


      fluidRow(uiOutput("volPlot"))
    )))


)
