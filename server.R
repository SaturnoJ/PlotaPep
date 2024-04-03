library(tidyr)
library(shiny)
library(shinyjs)
library(shinydashboard)
library("tidyverse")
library("org.Hs.eg.db") #if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("org.Hs.eg.db")
library("Biostrings")
library("arsenal")
library("shiny")
library("qdapRegex")
library("dplyr")
library("stringr")
library("seqinr")
library("rstatix")
library(data.table)
library(magrittr)
library(shinydisconnect)
library("colorspace")

source("functions.R", local = T)


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024 ^ 3)
  #Reactive Variables
  #######
  
  fasta_reactive <- reactiveVal()
  files <- reactiveVal()
  uniprot <- reactiveVal()
  intensity <- reactiveVal()
  lfq <- reactiveVal()
  cutoff <- reactiveVal()
  downloadData <- reactiveVal()
  df <- data.frame()
  protein_ids <- reactiveVal()
  j <- reactiveVal(0)
  located_peptides <- reactiveVal()
  plotting_protein <- reactiveVal()
  filtered_results <- reactiveVal()
  protein_plotr <- reactiveVal()
  filtered_results <- reactiveVal()
  plot_fasta <- reactiveVal()
  cohort_name <- reactiveVal()
  parametric <- reactiveVal()
  ctr_name <- reactiveVal()
  final_dataframe <- reactiveVal()
  removed_proteins <- reactiveVal()
  y_axis <- reactiveVal()
  p_cutoff <- reactiveVal()
  comparative <- reactiveVal()
  pal <- reactiveVal()
  colors_reactive <- reactiveVal()
  std <- reactiveVal()
  svg <- reactiveVal()
  less_than <- reactiveVal()
  file_type <- reactiveVal()
  #######
  #Fasta Input
  #######
  
  output$fastaInput <- renderUI({
    "txt"
    tagList(
      fileInput(
        "fastaInput",
        "Choose FASTA file : ",
        multiple = FALSE,
        accept = c("fasta", ".FASTA", "FASTA", ".fasta"),
        # close c
        buttonLabel = "Browse...",
        placeholder = "No FASTA Selected"
      )
    )
  })
  
  
  
  
  output$uniprotInput <- renderUI({
    "txt"
    
    tagList(
      textInput(
        "uniprotInput",
        "Input Uniprot IDs separated by commas (,) : ",
        placeholder = "Example P05067, P02649, etc."
      ),
      value = NULL
    )
  })
  
  
  observeEvent(input$fastaInput, {
    fasta_reactive(readAAStringSet(input$fastaInput$datapath))
    observeEvent(input$uniprotInput, {
      uniprotTemp <- input$uniprotInput
      uniprotTemp <- strsplit(uniprotTemp, ", |,| , ")
      uniprot(uniprotTemp)
      
    })
    shinyjs::enable("confirmFasta")
    
  })
  
  observeEvent(input$confirmFasta, {
    shinyjs::hide("fastaInput")
    shinyjs::hide("confirmFasta")
    shinyjs::hide("uniprotInput")
    
  })
  
  
  
  #######
  
  #File Input
  ######
  
  output$fileInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(
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
      )
    )
  })
  
  
  output$file_typeInput <- renderUI({
    "txt"
    
    req(input$confirmFasta)
    tagList(radioButtons(
      "file_typeInput",
      "Select MSFragger Output Type : ",
      c("No modifications " = 0, "Modified Ion " = 1)
    ))
  })
  
  output$lfqInput <- renderUI({
    "txt"
    
    req(input$confirmFasta)
    tagList(radioButtons("lfqInput", "Select LFQ Type : ", c("None" = 1, "Max" = 0)))
  })
  
  output$intensityInput <- renderUI({
    "txt"
    
    req(input$confirmFasta)
    tagList(radioButtons(
      "intensityInput",
      "Select Intensity : ",
      c(
        "Intensity" = 0,
        "Total" = 1,
        "Unique" = 2
      )
    ))
  })
  
  output$parametricInput <- renderUI({
    "txt"
    
    req(input$confirmFasta)
    tagList(
      radioButtons(
        "parametricInput",
        "Select Parametric or Nonparametric statistical test : ",
        c("Parametric" = 1, "Nonparametric" = 0)
      )
    )
  })
  output$svgInput <- renderUI({
    "txt"
    
    req(input$confirmFasta)
    tagList(radioButtons(
      "svgInput",
      "Select .SVG or .PNG file formats : ",
      c("PNG" = 0, "SVG" = 1)
    ))
    
  })
  output$comparativeInput <- renderUI({
    "txt"
    
    req(input$confirmFasta)
    tagList(
      radioButtons(
        "comparativeInput",
        "Select Comparative TTest Results (Compares multiple cohort's ttest results) : ",
        c("Noncomparative" = 0, "Comparative" = 1)
      )
    )
  })
  
  
  output$cutoffInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(numericInput(
      "cutoffInput",
      "Select Cutoff : ",
      0,
      min = 0,
      max = 100
    ))
  })
  
  output$pInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(numericInput(
      "pInput",
      "Select p value for significance : ",
      0.05,
      min = 0,
      max = 100
    ))
  })
  
  # output$yInput <- renderUI({
  #   "txt"
  #   req(input$confirmFasta)
  #   tagList(
  #     numericInput(
  #       "yInput",
  #       "Select Y-Axis (Will be used for both upregulation and downregulation) : ",
  #       2.5,
  #       min = 0,
  #       max = 100
  #     )
  #   )
  # })
  
  
  output$confirmFile <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(actionButton("confirmFile", "Confirm File(s) Input"))
  })
  
  output$ctrInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    
    tagList(
      textInput(
        "ctrInput",
        "Input control identifiers separated by commas : ",
        placeholder = "Example CTR, AD, DLB etc."
      ),
      value = NULL
      
    )
  })
  
  output$cohortInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    
    
    tagList(
      textInput(
        "cohortInput",
        "Input cohort identifiers separated by commas (,). [These should be the file names of the combined_petide.tsv] : ",
        
        placeholder = "Example AD, DLB etc."
      ),
      value = NULL
    )
  })
  
  
  output$colorsInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    
    tagList(
      textInput(
        "colorsInput",
        "Input hex value color codes :",
        placeholder = "Example CTR, AD, DLB etc.",
        value = "#E495A5,#86B875,#7DB0DD,#FFA500"
      ),
      value =  NULL
      
    )
  })
  output$stdInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(numericInput(
      "stdInput",
      "Select standard deviations : ",
      2,
      min = 2,
      max = 6
    ))
  })
  #######
  
  #File Observe
  ########
  observeEvent(input$fileInput, {
    files <- input$fileInput
    df <- list()
    for (i in 1:nrow(files)) {
      df[[i]] <-
        as.data.frame(read_delim(files$datapath[i], show_col_types = FALSE))
      
      df[[i]]$File <- str_sub(files$name[i], end = -5)
      
    }
    files(df)
    
  })
  
  observeEvent(input$file_typeInput, {
    file_type(input$file_typeInput)
    
  })
  
  observeEvent(input$intensityInput, {
    intensity(input$intensityInput)
    
  })
  observeEvent(input$lfqInput, {
    lfq(input$lfqInput)
    
    
  })
  observeEvent(input$parametricInput, {
    parametric(input$parametricInput)
    
  })
  
  
  observeEvent(input$cutoffInput, {
    cutoff(input$cutoffInput)
    
  })
  
  
  # observeEvent(input$yInput, {
  #   y_axis(input$yInput)
  #
  # })
  
  observeEvent(input$pInput, {
    p_cutoff(input$pInput)
    
  })
  observeEvent(input$stdInput, {
    std(input$stdInput)
    
  })
  observeEvent(input$svgInput, {
    svg(input$svgInput)
    
  })
  
  observeEvent(input$comparativeInput, {
    comparative(input$comparativeInput)
    
  })
  
  observeEvent(input$ctrInput, {
    name_temp <- input$ctrInput
    # name_temp <- strsplit(name_temp, ", |,| , ")
    ctr_name(unlist(name_temp))
    
  })
  
  observeEvent(input$cohortInput, {
    name_temp <- input$cohortInput
    name_temp <- strsplit(name_temp, ", |,| , ")
    cohort_name(unlist(name_temp))
    
  })
  
  observeEvent(input$colorsInput, {
    colors_temp <- input$colorsInput
    colors_temp <- strsplit(colors_temp, ", |,| , ")
    colors_reactive(unlist(colors_temp))
    
    
    
    
  })
  
  
  
  ########
  
  
  
  
  
  
  observeEvent(input$confirmFile, {
    #UI LOGIC
    #######
    shinyjs::hide("fileInput")
    shinyjs::hide("file_typeInput")
    shinyjs::hide("confirmFile")
    shinyjs::hide("lfqInput")
    shinyjs::hide("cutoffInput")
    shinyjs::hide("intensityInput")
    shinyjs::hide("cohortInput")
    shinyjs::hide("ctrInput")
    shinyjs::hide("parametricInput")
    shinyjs::hide("pInput")
    # shinyjs::hide("yInput")
    shinyjs::hide("comparativeInput")
    shinyjs::hide("colorsInput")
    shinyjs::hide("stdInput")
    shinyjs::hide("svgInput")
    
    #######
    
    #Progress bar
    #######
    progress <- shiny::Progress$new()
    progress$set(message = "Preparing Data", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 10
      }
      progress$set(value = value, detail = detail)
    }
    #######
    #Variable Declaration
    #######
    combined_cutoff_df <- data.frame()
    comparative_combined <- data.frame()
    less_than_two <- data.frame()
    uniprot_ids <- as.data.frame(uniprot())
    protein_plot <- list()
    colnames(uniprot_ids)  <- c("Protein.ID")
    cohorts <- toupper(cohort_name())
    parametric_type <- parametric()
    ctr <- toupper(ctr_name())
    protein_id_loop <- list
    colors_df <- as.list(colors_reactive())
    fasta <- as.data.frame(fasta_reactive())
    fasta <- tibble::rownames_to_column(fasta, "accession")
    suppressWarnings(
      fasta %<>%
        separate(accession, sep = '\\|', c("sp", "Protein.ID", "info")) %>%
        separate(info, sep = 'GN=', c("info", "Gene")) %>%
        separate(Gene, sep = " ", c("Gene", "other"))
    )
    df <- files()
    #######
    #Dateframe creation
    #######
    if (comparative() == 0) {
      for (i in df) {
        for (j in 1:length(cohorts)) {
          temp <-
            cleanData(
              i,
              as.integer(intensity()),
              as.integer(lfq()),
              as.integer(cutoff()),
              as.integer(file_type()),
              updateProgress
            )
          
          temp <-
            cohortSplit(temp,
                        ctr,
                        cohorts[j],
                        as.integer(std()),
                        updateProgress)
          
          combined_cutoff_df <- rbind(combined_cutoff_df, temp)
        }
        
        possible <-
          combined_cutoff_df %>%  spread(Cohort, Intensity) %>% group_by(Protein) %>% summarise(ctr = sum(!is.na(get("CTR"))),
                                                                                                cohort1 = sum(!is.na(get(cohorts[j])))) %>%
          mutate(possible = ifelse(ctr < 2 |
                                     cohort1 < 2, FALSE, TRUE)) %>%
          filter(possible)
        
        cohort_df <-
          combined_cutoff_df %>% filter(Protein %in% possible$Protein)
        
        # cohort_df <-
        #   combined_cutoff_df[duplicated(combined_cutoff_df$Protein), ]
        less_than_two <-
          anti_join(cohort_df, combined_cutoff_df, by = "Protein")
        
        
        if (nrow(cohort_df) < 1) {
          session$close()
          
        }
        ########
        
        cohort_df$Intensity <- log2(cohort_df$Intensity)
        
        #Nonparametric
        #####
        if (parametric_type == 0 & length(cohorts) == 1) {
          statisical_test <-
            getMannWhit(cohort_df, p_cutoff(), updateProgress)
          
          
        }
        
        
        else if (parametric_type == 0) {
          statisical_test <- getKruskal(cohort_df, updateProgress())
          
        }
        
        #####
        
        #Parametric
        #######
        
        
        else{
          if (length(cohorts) > 1) {
            statisical_test <- getANOVA(cohort_df, updateProgress)
            
            
          }
          else{
            statisical_test <- getTtest(cohort_df, p_cutoff(), updateProgress)
            
            
          }
          
        }
        
        ########
        
      }
    }
    
    else{
      for (i in df) {
        #Comparative
        #########
        for (j in 1:length(cohorts)) {
          temp <-
            cleanData(
              i,
              as.integer(intensity()),
              as.integer(lfq()),
              as.integer(cutoff()),
              as.integer(file_type()),
              updateProgress
            )
          temp[, 1] <- toupper(temp[, 1])
          if (!isEmpty(temp[temp$Sample %like% toupper(cohorts[j]),])) {
            cohort_df <-
              cohortSplit(temp,
                          ctr,
                          cohorts[j],
                          as.integer(std()),
                          updateProgress)
            
          }
          else{
            next
          }
          
          possible <-
            cohort_df %>%  spread(Cohort, Intensity) %>% group_by(Protein) %>% summarise(ctr = sum(!is.na(get("CTR"))),
                                                                                         cohort1 = sum(!is.na(get(cohorts[j])))) %>%
            mutate(possible = ifelse(ctr < 2 |
                                       cohort1 < 2, FALSE, TRUE)) %>%
            filter(possible)
          
          cohort_df <-
            cohort_df %>% filter(Protein %in% possible$Protein)
          # cohort_df[cohort_df$Protein %in% cohort_df$Protein[duplicated(cohort_df$Protein)], ]
          cohort_df$Intensity <- log2(cohort_df$Intensity)
          
          
          if (parametric_type == 0) {
            statisical_test <-
              getMannWhit(cohort_df, p_cutoff(), updateProgress)
            
            
            
          }
          else{
            statisical_test <- getTtest(cohort_df, p_cutoff(), updateProgress)
            
          }
          temp <-  anti_join(cohort_df, temp, by = "Protein")
          
          less_than_two <- rbind(less_than_two, temp)
          comparative_combined <-
            rbind(comparative_combined, statisical_test)
        }
      }
      
      if (is.function(updateProgress)) {
        text <- "Removing Outliers"
        updateProgress(detail = text)
      }
    }
    
    #Plotting loop
    ######
    if (comparative() == 0) {
      leftovers <-
        anti_join(statisical_test, fasta, by = "Protein.ID")
      
      statisical_test %<>% inner_join(fasta, by = "Protein.ID", multiple = "all") %>% arrange(Protein.ID)
      
      located_peptides(locatePeptides(statisical_test, updateProgress))
      #######
      if (nrow(uniprot_ids) > 0) {
        protein_id_loop <- uniprot_ids
      }
      else{
        protein_id_loop <-  unique(statisical_test$Protein.ID)
      }
      for (i in protein_id_loop) {
        plotting_protein(filter(located_peptides(), Protein.ID == i))
        filtered_results(filter(located_peptides(), Protein.ID == i))
        filtered_results(inner_join(
          fasta,
          plotting_protein(),
          by  = c('Protein.ID', 'x', 'Gene', 'info'),
          multiple = "all"
        ))
        protein_length <- unique(nchar(filtered_results()$x))
        
        protein_plot[[i]] = plotting_protein() %>%
          
          #Lets make a column based on significance
          mutate(isSignificant = ifelse(filtered_results()$p.adj < p_cutoff(), "yes", "no")) %>%
          
          #Start plotting
          ggplot() +
          
          #We take here the 'mean' but this is of course X-times the same value
          scale_y_continuous(limits = c(NA, NA)) +
          xlim(0, protein_length + 1) +
          
          #Create some lines to help visualise the start and end of the protein.
          geom_vline(xintercept = 0,
                     lwd = 2,
                     alpha = .5) +
          geom_vline(xintercept = protein_length + 1,
                     lwd = 2,
                     alpha = 0.5) +
          
          #set a horizontal line to inform about FC = 0
          geom_hline(yintercept = 0, color = "black") +
          
          #Here we build the peptide "blocks". Do note that all column-info goes INSIDE the aes() part.
          geom_rect(
            inherit.aes = FALSE,
            aes(
              xmin = plotting_protein()$start_seq,
              xmax = (
                plotting_protein()$start_seq + (
                  plotting_protein()$end_seq - plotting_protein()$start_seq
                )
              ),
              ymin = filtered_results()$FC - 0.05,
              ymax = filtered_results()$FC + 0.05,
              fill = isSignificant
             
            ),
            
            #Here I specify some stuff that will be universal for all blocks irrespective of column info.
            col = "black",
            alpha = 0.75,
            linewidth = .15,
          ) +
          
          #Set the Ggplot theme, limit the y-axis for FC.
          theme_bw() +
          
          #Specify the colours I want to use for the isSignificant column
          scale_color_manual(values = c("yes" = "black", "no" = "red")) +
          theme(legend.position = "bottom") +
          
          #x and yaxis titles
          xlab("Protein Sequence") +
          ylab("FC") +
          labs(
            subtitle = paste0(
              as.character(i),
              " p == ",
              p_cutoff(),
              " Standard Deviations ==",
              std()
            ),
            title = as.character(filtered_results()$Gene)
          )
        
        
        
        
      }
    }
    ######
    else{
      print("hello")
      leftovers <-
        anti_join(comparative_combined, fasta, by = "Protein.ID")
      comparative_combined %<>% inner_join(fasta, by = "Protein.ID", multiple = "all") %>% arrange(Protein.ID) %>% mutate(
        colors = case_when(
          Cohort == as.character(cohorts[1]) ~  as.character(colors_df[1]),
          Cohort == as.character(cohorts[2]) ~ as.character(colors_df[2]),
          Cohort == as.character(cohorts[3]) ~ as.character(colors_df[3]),
          Cohort == as.character(cohorts[4]) ~ as.character(colors_df[4])
        )
      )
      
      
      located_peptides(locatePeptides(comparative_combined, updateProgress))
      colors <- distinct(comparative_combined, Cohort, colors)
      pal_temp <- colors$color
      names(pal_temp) <- colors$Cohort
      pal(pal_temp)
      
      #Plotting loop
      #######
      if (nrow(uniprot_ids) > 0) {
        protein_id_loop <- uniprot_ids
        
      }
      else{
        protein_id_loop <-  unique(comparative_combined$Protein.ID)
      }
      # browser()
      for (i in protein_id_loop)  {
        plotting_protein(filter(located_peptides(), Protein.ID == i))
        filtered_results(filter(located_peptides(), Protein.ID == i))
        filtered_results(inner_join(
          fasta,
          plotting_protein(),
          by  = c('Protein.ID', 'x', 'Gene', 'info'),
          multiple = "all"
        ))
        protein_length <- unique(nchar(filtered_results()$x))
        plotting_protein(
          plotting_protein() %>%
            
            #Lets make a column based on significance
            mutate(
              isSignificant = ifelse(filtered_results()$p.adj < p_cutoff(), "yes", "no")
            )
          %>% mutate(color = ifelse(
            isSignificant == "yes",
            darken(plotting_protein()$color, 0.4),
            lighten(plotting_protein()$color, 0.2)
          ))
          
          %>% mutate(
            Cohort = ifelse(
              isSignificant == "yes",
              paste0(plotting_protein()$Cohort, "_Significant"),
              plotting_protein()$Cohort
            )
          )
          # %>% 
          #   column_to_rownames(var = "Cohort")
        )
        
        # names(pal_temp) <- colors$Cohort
        # pal(pal_temp)
        protein_plot[[i]] = plotting_protein() %>%
          
          #Start plotting
          ggplot() +
          
          #We take here the 'mean' but this is of course X-times the same value
          scale_y_continuous(limits = c(NA, NA)) +
          xlim(0, protein_length + 1) +
          
          #Create some lines to help visualise the start and end of the protein.
          geom_vline(xintercept = 0,
                     lwd = 2,
                     alpha = .5) +
          geom_vline(xintercept = protein_length + 1,
                     lwd = 2,
                     alpha = 0.5) +
          
          #set a horizontal line to inform about FC = 0
          geom_hline(yintercept = 0, color = "black") +
          
          #Here we build the peptide "blocks". Do note that all column-info goes INSIDE the aes() part.
          geom_rect(
            #Here I specify some stuff that will be universal for all blocks irrespective of column info.
            
            aes(
              xmin = plotting_protein()$start_seq,
              xmax = (
                plotting_protein()$start_seq + (
                  plotting_protein()$end_seq - plotting_protein()$start_seq
                )
              ),
              ymin = filtered_results()$FC - 0.05,
              ymax = filtered_results()$FC + 0.05,
              #fill = isSignificant,
              fill = color,
            ),
            alpha = 0.75,
            linewidth = .25,
            
          ) +
          
          #Set the Ggplot theme, limit the y-axis for FC.
          theme_bw() +
          
          #Specify the colours I want to use for the isSignificant column
          #scale_fill_manual(values = c("yes" = darken(plotting_protein()$color,0.4), "no" = lighten(plotting_protein()$color,0.2))) +
          scale_fill_identity("Cohort",  labels = setNames(unique(plotting_protein()$Cohort), unique(plotting_protein()$color)), guide = "legend") +
          theme(legend.position = "bottom",
                legend.key = element_rect(colour = "white")) +
          guides(fill = guide_legend(order = 1),
                 color = guide_legend(order = 2)) +
          #x and yaxis titles
          xlab("Protein Sequence") +
          ylab("FC") +
          
          labs(
            subtitle = paste0(
              as.character(i),
              " p == ",
              p_cutoff(),
              " standard deviations ==",
              std()
            ),
            title = as.character(filtered_results()$Gene)
          )
        
        
        
      }
      
      
      
    }
    #Reactive Variables
    ######
    protein_plotr(protein_plot)
    plot_fasta(fasta)
    if (nrow(uniprot_ids) == 0) {
      # browser()
      if (comparative() == 0) {
        temp_id <- as.data.frame(unique(statisical_test$Protein.ID))
        final_dataframe(located_peptides())
        
      }
      else{
        temp_id <- as.data.frame(unique(comparative_combined$Protein.ID))
        final_dataframe(located_peptides())
        
      }
      colnames(temp_id) <- c("Protein.ID")
      protein_ids(temp_id)
    }
    else{
      protein_ids(uniprot_ids)
    }
    removed_proteins(leftovers)
    less_than(less_than_two)
    final_dataframe(located_peptides())
    ######
    
    
    shinyjs::enable("clickThrough")
    shinyjs::enable("downloadPlot")
    
    
    
    
    
    
    
    
    
    
  })
  
  output$downloadPlot <- downloadHandler(
    filename = 'plots.zip',
    content = function(file) {
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      saved <- data.frame(Protein.ID = character())
      
      plots <- as.list(protein_plotr())
      if (as.integer(svg()) == 0) {
        for (i in unique(protein_ids()$Protein.ID)) {
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          if (nrow(plotting_protein()) > 0) {
            filtered_results(filter(located_peptides(), Protein.ID == i))
            filtered_results(inner_join(
              plot_fasta(),
              plotting_protein(),
              by  = c('Protein.ID', 'x', 'Gene', 'info'),
              multiple = "all"
            ))
            
            protein_length <- unique(nchar(filtered_results()$x))
            
            
            ggsave(
              paste(
                as.character(filtered_results()$Gene[1]),
                "_",
                as.character(filtered_results()$Protein.ID[1]),
                ".png",
                sep = ""
              ),
              plot = plots[[i]],
              device = "png"
            )
            
          }
        }
      }
      
      else{
        for (i in unique(protein_ids()$Protein.ID)) {
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          if (nrow(plotting_protein()) > 0) {
            filtered_results(filter(located_peptides(), Protein.ID == i))
            filtered_results(inner_join(
              plot_fasta(),
              plotting_protein(),
              by  = c('Protein.ID', 'x', 'Gene', 'info'),
              multiple = "all"
            ))
            
            protein_length <- unique(nchar(filtered_results()$x))
            
            
            ggsave(
              paste(
                as.character(filtered_results()$Gene[1]),
                "_",
                as.character(filtered_results()$Protein.ID[1]),
                ".svg",
                sep = ""
              ),
              plot = plots[[i]],
              device = "svg"
            )
          }
        }
        
      }
      
      write.csv(final_dataframe()
                , file = "final_dataframe.csv"
                , row.names = F)
      write.csv(removed_proteins()
                ,
                file = "removed_proteins_fasta.csv"
                ,
                row.names = F)
      write.csv(less_than()
                ,
                file = "removed_proteins_less_two.csv"
                ,
                row.names = F)
      
      zip_files <-
        list.files(path = getwd(), pattern = ".svg$|.csv$|.png$")
      
      zip::zip(file, files = zip_files)
      files_to_delete <-
        dir(path = getwd() , pattern = "*.png$|*.svg$|*.csv$")
      file.remove(file.path(getwd(), files_to_delete))
      
    }
    
  )
  
  
  
  
  observeEvent(input$clickThrough, {
    plots <- as.list(protein_plotr())
    ids <- as.data.frame(protein_ids())
    if (j() < nrow(ids)) {
      output$plot <- renderUI({
        i <-  unique(protein_ids()[j(), 'Protein.ID'])
        
        plotting_protein(filter(located_peptides(), Protein.ID == i))
        filtered_results(filter(located_peptides(), Protein.ID == i))
        filtered_results(inner_join(
          plot_fasta(),
          plotting_protein(),
          by  = c('Protein.ID', 'x', 'Gene', 'info'),
          multiple = "all"
        ))
        
        
        protein_length <- unique(nchar(filtered_results()$x))
        
        
        renderPlot({
          plots[[i]]
          
        })
        
      })
      j(j() + 1)
      
    }
    else {
      output$plot <- renderUI({
        i <-  unique(ids[j(), 'Protein.ID'])
        
        
        
        plotting_protein(filter(located_peptides(), Protein.ID == i))
        filtered_results(filter(located_peptides(), Protein.ID == i))
        filtered_results(inner_join(
          plot_fasta(),
          plotting_protein(),
          by  = c('Protein.ID', 'x', 'Gene', 'info'),
          multiple = "all"
        ))
        
        protein_length <- unique(nchar(filtered_results()$x))
        
        
        renderPlot({
          plots[[i]]
        })
        
      })
      j(1)
      
    }
    
    
  })
  
  
  
  
  
}