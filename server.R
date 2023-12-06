library(tidyr)
library(shiny)
library(shinyjs)
library(shinydashboard)
library("tidyverse")
library("org.Hs.eg.db")
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
library(cowplot)



# clean the data
#################
cleanData <-
  function(df,
           intensity,
           lfq,
           cutoff,
           updateProgress = NULL) {
    if (is.function(updateProgress)) {
      text <- "Cleaning Data"
      updateProgress(detail = text)
    }
    
    file <- df$file[1]
    
    cleaned_df <-
      fraggerCleaner(df, intensity, lfq, cutoff)
    
    if (is.function(updateProgress)) {
      text <- "Removing Outliers"
      updateProgress(detail = text)
    }
    
    outlier_removed_df <- removeOutlier(cleaned_df)
    
    
    
    combined_cutoff_df <- outlier_removed_df %>%
      mutate(uniqueSamplesInDF = length(unique(Sample))) %>%
      group_by(Protein) %>%
      mutate(countMissingValues = sum(!is.na(Intensity))) %>%
      mutate(percentageMissing = countMissingValues / uniqueSamplesInDF) %>%
      filter(percentageMissing > 0.95) %>% # Apply cutoff of 95%
      dplyr::select(Sample, Protein, Intensity) #clean up by selecting the original columns
    
    combined_cutoff_df$File <- file
    
    return(as.data.frame(combined_cutoff_df))
  }


fraggerCleaner <-
  function(df, intensity, lfq, cutoff) {

        # change Protein name
    if (grepl('Peptide Sequence', names(df))) {
      names(df)[names(df) == 'Peptide Sequence'] <-
        'Sequence'
    }
    
    #Select columns of interest
    df <-
      unite(df,
            Protein,
            c(Protein, Sequence),
            sep = "|",
            remove = FALSE)
    # Select the specified Intensity.
    # We also start creating the substring we will remove from sample names.
    if (intensity == 0) {
      df <-
        dplyr::select(df,
                      Protein,
                      contains("Intensity"),
                      -contains("Total"),
                      -contains("Unique"))
      removeSubString <- paste0(" ", "Intensity")
    } else if (intensity == 1) {
      df <-
        dplyr::select(df, Protein, contains("Total Intensity"))
      removeSubString <- paste0(" ", "Total", " Intensity")
    } else {
      df <-
        dplyr::select(df, Protein, contains("Unique Intensity"))
      removeSubString <- paste0(" ", "Unique", " Intensity")
    }
    
    #select or not select MaxLFQ
    if (lfq == 0) {
      df <- dplyr::select(df, Protein, contains("MaxLFQ"))
      removeSubString <- paste0(" MaxLFQ", removeSubString)
    } else {
      df <-
        dplyr::select(df,
                      Protein,
                      contains("Intensity"), -contains("MaxLFQ"))
    }
    
    #col to rowname
    df <- data.frame(df[,-1], row.names = df[, 1])
    
    
    #Get rid of substring part we do not need
    
    colnames(df) <-
      stringr::str_remove(colnames(df), removeSubString)
    
    
    #transpose, remove zeros, log2
    df_transposed <-  transpose(df)
    
    rownames(df_transposed) <- colnames(df)
    colnames(df_transposed) <- rownames(df)
    df <- df_transposed
    
    df <- removeZero(df)
    df <- log2(df)
    
    
    
    #apply cutoff on the protein level
    
    df <-
      purrr::discard(df, ~ sum(is.na(.x)) / length(.x) * 100 >= (100 - cutoff))
    
    
    #Long Format it
    df <- setDT(df, keep.rownames = "Sample")
    df <- tidyr::gather(df,
                        colnames(df)[2:ncol(df)],
                        key = "Protein",
                        value = "Intensity")
    
    
    return(df)
  }


removeOutlier <- function(df) {
  #Determine correlations between sample and theoretical
  SampleCorrelations <- df %>%
    group_by(Protein) %>%
    dplyr::mutate(medianProtein = median(Intensity, na.rm = T)) %>%
    ungroup() %>%
    group_by(Sample) %>%
    dplyr::mutate(
      SampleCor = stats::cor(Intensity, medianProtein, method = "pearson", use =
                               "pairwise.complete.obs")
    ) %>%
    ungroup() %>%
    distinct(Sample, .keep_all = T)
  
  SampleCorrelations <- SampleCorrelations %>%
    dplyr::mutate(Outlier = SampleCor < (1 - 1 * sd(
      unique(SampleCorrelations$SampleCor), na.rm = T
    )))
  
  #Make the DF that we will return (i.e. without outliers)
  returndf <- df %>%
    dplyr::filter(Sample %in% dplyr::filter(SampleCorrelations, Outlier == F)$Sample) %>%
    dplyr::select(Protein, Sample, Intensity)
  
  #Half minimum value (per protein) imputation
  PCA_DF <- df %>%
    group_by(Protein) %>%
    dplyr::mutate(IntensImputed = replace_na(Intensity, mean(Intensity, na.rm = T) /
                                               2)) %>%
    dplyr::select(Protein, Sample, IntensImputed) %>%
    spread(key = "Protein", value = "IntensImputed") %>%
    column_to_rownames("Sample")
  
  #scale and PCA
  PCA_DF_Results <- prcomp(scale(as.matrix(PCA_DF)))
  eigs <- PCA_DF_Results$sdev ^ 2
  variance_percentage <- (eigs / sum(eigs)) * 100
  pc1var <- round(variance_percentage[1], digits = 0)
  pc2var <- round(variance_percentage[2], digits = 0)
  
  return(returndf)
}

removeZero <- function(df) {
  for (j in seq_len(ncol(df))) {
    set(df, which((df[[j]]) == 0), j, NA)
  }
  return(df)
}

#################

# Cohort Split
########
cohortSplit <-
  function(combined_cutoff_df,
           cohort,
           updateProgress = NULL) {
    if (is.function(updateProgress)) {
      text <- "Creating Cohorts"
      updateProgress(detail = text)
    }
    browser()
    cohort_df <- data.frame()
    
    for (i in cohort) {
      temp_df <- splitCtr(combined_cutoff_df, i)
      cohort_df <- cbind(temp_df, i)
      
    }
    
    
    
    if (is.function(updateProgress)) {
      text <- "Cohort Created"
      updateProgress(detail = text)
    }
    
    return(as.data.frame(cohort_df))
  }


##extracts the alzheimers df from the merged df
splitDisease <- function(df, ctr) {
  df <- subset(df,!grepl("CTR", df$Sample))
  if (nrow(df) < 0) {
    return(NULL)
  }
  return(df)
  
}

##extracts the control df from the merged df
splitCtr <- function(combined_cutoff_df, cohort) {
  df <-
    subset(combined_cutoff_df,
           grepl(cohort, combined_cutoff_df$Sample))
  if (nrow(df) < 0) {
    return(NULL)
  }
  return(df)
}


###########


#TTest
########

getTtest <- function(cohort_df, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    text <- "Performing Ttest (This may take a while)"
    updateProgress(detail = text)
  }
  
  ttest_results <-  cohort_df %>%
    
    
    #Run the T-test and adjustments
    group_by(Protein) %>%
    t_test(Intensity ~ DX, detailed = T) %>%
    adjust_pvalue(method = "BH") %>%
    
    #Split the Protein name in Uniprot and Gene
    separate(Protein, sep = '\\|', c("SP", "Protein.ID", "Gene", "Peptide"))  %>%
    
    #Determine Fold change. Since we work with log-transformed values we can just substract
    mutate(FC = estimate1 - estimate2) %>%
    
    #Create log10 p-vals
    mutate(log10adjustP = -1 * log10(p.adj)) %>%
    
    #Determine if up or down regulated
    mutate(Direction = ifelse(p.adj > 0.05, "NotSignificant", ifelse(FC < 0, "Down", "Up")))
  
  ##Locate the peptides in the protein sequence and add the start and end point to the dataframe
  # located_peptides() <- ttest_results[(grepl(paste(dfIds),ttest_results$UniprotID))]
  ttest_results <-
    cbind(ttest_results,
          start_seq = NA,
          end_seq = NA)
  
  if (is.function(updateProgress)) {
    text <- "Ttest Complete"
    updateProgress(detail = text)
  }
  
  return(as.data.frame(ttest_results))
}


#########







locatePeptides <- function(df,
                           updateProgress = NULL) {
  if (is.function(updateProgress)) {
    text <- "Matching peptides to proteins"
    updateProgress(detail = text)
  }
  
  for (i in 1:nrow(df)) {
    postion <-
      as.data.frame(str_locate(df$x[i], df$Peptide[i]))
    if (nrow(postion) != 0) {
      df$start_seq[i] <- postion$start
      df$end_seq[i] <- postion$end
    }
    else{
      df$start_seq[i] <- NA
      df$end_seq[i] <- NA
    }
    
  }
  
  if (is.function(updateProgress)) {
    text <- "Peptides matched"
    updateProgress(detail = text)
  }
  
  return(df)
  
}




server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024 ^ 3)
  
  fastaR <- reactiveVal()
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
  
  #Fasta Input
  #######
  
  output$fastaInput <- renderUI({
    "txt"
    tagList(
      fileInput(
        "fastaInput",
        "Choose FASTA file",
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
        "Input Uniprot IDs ated by commas (,)",
        placeholder = "Example P05067, P02649, etc."
      ),
      value = NULL
    )
  })
  
  
  observeEvent(input$fastaInput, {
    fastaR(readAAStringSet(input$fastaInput$datapath))
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
        "Choose file",
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
  
  
  output$lfqInput <- renderUI({
    "txt"
    
    req(input$confirmFasta)
    tagList(radioButtons("lfqInput", "Select LFQ Type: ", c("None" = 1,"Max" = 0)))
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
  
  output$cutoffInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(numericInput(
      "cutoffInput",
      "Select Cutoff : ",
      99,
      min = 0,
      max = 100
    ))
  })
  
  
  output$confirmFile <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(disabled(actionButton(
      "confirmFile", "Confirm File(s) Input"
    )))
  })
  
  
  output$cohortInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    
    tagList(
      textInput(
        "cohortInput",
        "Input cohort identifiers seperated by commas (,). [These should be in the filename for the Fragger search]",
        placeholder = "Example CTR, AD, DLB etc."
      ),
      value = NULL
    )
  })
  
  
  ######
  
  observeEvent(input$fileInput, {
    files <- input$fileInput
    df <- list()
    for (i in 1:nrow(files)) {
      df[[i]] <-
        as.data.frame(read_delim(files$datapath[i], show_col_types = FALSE))
      
      df[[i]]$fileName <- str_sub(files$name[i], end = -5)
      
      
      
      
    }
    files(df)
    
  })
  
  
  observeEvent(input$intensityInput, {
    intensity(input$intensityInput)
    
  })
  observeEvent(input$lfqInput, {
    lfq(input$lfqInput)
    
    
  })
  observeEvent(input$cutoffInput, {
    cutoff(input$cutoffInput)
    
  })
  
  
  observeEvent(input$cohortInput, {
    nameTemp <- input$cohortInput
    nameTemp <- strsplit(nameTemp, ", |,| , ")
    cohort_name(unlist(nameTemp))
    
    if (!is.null(input$fileInput))
      shinyjs::enable("confirmFile")
    
    
  })
  
  
  
  
  
  
  
  ###### end file screen ######
  
  
  observeEvent(input$confirmFile, {
    #UI LOGIC
    #######
    shinyjs::hide("fileInput")
    shinyjs::hide("confirmFile")
    shinyjs::hide("lfqInput")
    shinyjs::hide("cutoffInput")
    shinyjs::hide("intensityInput")
    shinyjs::hide("cohortInput")
    shinyjs::hide("cohortNInput")
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
    uniprot_ids <- as.data.frame(uniprot())
    protein_plot <- list()
    colnames(uniprot_ids)  <- c("Protein.ID")
    cohorts <- cohort_name()
    
    fasta <- as.data.frame(fastaR())
    fasta <- tibble::rownames_to_column(fasta, "accession")
    fasta %<>%
      separate(accession, sep = '\\|', c("sp", "Protein.ID", "info"))
    
    df <- files()
    #######
    
    for (i in df) {
      
      temp <-
        cleanData(
          i,
          as.integer(intensity()),
          as.integer(lfq()),
          as.integer(cutoff()),
          updateProgress
        )
      
      combined_cutoff_df <- rbind(combined_cutoff_df, temp)
      
    }
    
    browser()
     
    cohort_df <-
      cohortSplit(combined_cutoff_df, cohorts, updateProgress)
    
    ttest_results <- getTtest(cohort_df, updateProgress)
    
    leftovers <- anti_join(ttest_results, fasta, by = "Protein.ID")
    
    ttest_results %<>% inner_join(fasta, by = "Protein.ID", multiple = "all") %>% arrange(Protein.ID)
    
    located_peptides(locatePeptides(ttest_results, updateProgress))
    
    
    #Plotting loop
    #######
    for (i in unique(uniprot_ids$Protein.ID)) {
      plotting_protein(filter(located_peptides(), Protein.ID == i))
      filtered_results(filter(located_peptides(), Protein.ID == i))
      filtered_results(inner_join(
        fasta,
        plotting_protein(),
        by  = c('Protein.ID', 'x'),
        multiple = "all"
      ))
      
      
      protein_length <- unique(nchar(filtered_results()$x))
      
      
      protein_plot[[i]] = plotting_protein() %>%
        
        #Lets make a column based on significance
        mutate(isSignificant = ifelse(filtered_results()$p.adj < 0.05, "yes", "no")) %>%
        
        #Start plotting
        ggplot() +
        
        #We take here the 'mean' but this is of course X-times the same value
        #ylim(-2.5, 2.5)
        xlim(1, protein_length) +
        
        #Create some lines to help visualise the start and end of the protein.
        geom_vline(xintercept = 1,
                   lwd = 2,
                   alpha = .5) +
        geom_vline(xintercept = protein_length,
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
        ) +
        
        #Set the Ggplot theme, limit the y-axis for FC.
        theme_bw() +
        ylim(-2, 2) +
        
        #Specify the colours I want to use for the isSignificant column
        scale_fill_manual(values = c("yes" = "red", "no" = "grey")) +
        theme(legend.position = "bottom") +
        
        #x and yaxis titles
        xlab("Protein Sequence") +
        ylab("FC") +
        
        labs(title = as.character(i))
      
      
      
      
      
    }
    ######
    
    
    #Reactive Variables
    ######
    protein_plotr(protein_plot)
    plot_fasta(fasta)
    protein_ids(uniprot_ids)
    ######
    
    
    shinyjs::enable("clickThrough")
    shinyjs::enable("downloadPlot")
    
    
    
    
    
    
    
    
    
    
  })
  
  output$downloadPlot <- downloadHandler(
    filename = 'plots.zip',
    content = function(file) {
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      
      
      
      browser()
      
      plots <- as.list(protein_plotr())
      
      for (i in unique(protein_ids()$Protein.ID)) {
        plotting_protein(filter(located_peptides(), Protein.ID == i))
        if (nrow(plotting_protein()) > 0) {
          filtered_results(filter(located_peptides(), Protein.ID == i))
          filtered_results(inner_join(
            plot_fasta(),
            plotting_protein(),
            by  = c('Protein.ID', 'x'),
            multiple = "all"
          ))
          
          protein_length <- unique(nchar(filtered_results()$x))
          
          
          ggsave(paste(as.character(i), ".png", sep = ""),
                 plot = plots[[i]],
                 device = "png")
        }
      }
      
      zip_files <- list.files(path = getwd(), pattern = ".png$")
      
      zip::zip(file, files = zip_files)
      
      
    }
    
  )
  
  
  
  
  observeEvent(input$clickThrough, {
    plots <- as.list(protein_plotr())
    ids <- as.data.frame(protein_ids())
    if (j() < nrow(as.data.frame(uniprot()))) {
      output$plot <- renderUI({
        i <-  unique(ids[j(), 'Protein.ID'])
        
        plotting_protein(filter(located_peptides(), Protein.ID == i))
        filtered_results(filter(located_peptides(), Protein.ID == i))
        filtered_results(inner_join(
          plot_fasta(),
          plotting_protein(),
          by  = c('Protein.ID', 'x'),
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
          by  = c('Protein.ID', 'x'),
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