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
library(shinydisconnect)

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
    
    file <- df$File[1]
    
    cleaned_df <-
      fraggerCleaner(df, intensity, lfq, cutoff)
    
    if (is.function(updateProgress)) {
      text <- "Removing Outliers"
      updateProgress(detail = text)
    }
    
    outlier_removed_df <- removeOutlier(cleaned_df)
    
    
    # combined_cutoff_df <- outlier_removed_df %>%
    #   mutate(uniqueSamplesInDF = length(unique(Sample))) %>%
    #   group_by(Protein) %>%
    #   mutate(countMissingValues = sum(!is.na(Intensity))) %>%
    #   mutate(percentageMissing = countMissingValues / uniqueSamplesInDF) %>%
    #   filter(percentageMissing > 0.95) %>% # Apply cutoff of 95%
    #   dplyr::select(Sample, Protein, Intensity) #clean up by selecting the original columns
    #
    #
    # combined_cutoff_df$File <- file
    #
    #
    
    
    outlier_removed_df$File <- file
    
    return(as.data.frame(outlier_removed_df))
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
                      contains("Intensity"),-contains("Total"),-contains("Unique"))
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
                      contains("Intensity"),-contains("MaxLFQ"))
    }
    
    #col to rowname
    df <- data.frame(df[, -1], row.names = df[, 1])
    
    
    #Get rid of substring part we do not need
    
    colnames(df) <-
      stringr::str_remove(colnames(df), removeSubString)
    
    
    #transpose, remove zeros, log2
    df_transposed <-  data.frame(t(df))
    
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
  sample_correlations <- df %>%
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
  
  sample_correlations <- sample_correlations %>%
    dplyr::mutate(Outlier = SampleCor < (1 - 1 * sd(
      unique(sample_correlations$SampleCor), na.rm = T
    )))
  
  #Make the DF that we will return (i.e. without outliers)
  returndf <- df %>%
    dplyr::filter(Sample %in% dplyr::filter(sample_correlations, Outlier == F)$Sample) %>%
    dplyr::select(Protein, Sample, Intensity)
  
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
           ctr,
           updateProgress = NULL) {
    if (is.function(updateProgress)) {
      text <- "Creating Cohorts"
      updateProgress(detail = text)
    }
    br
    
    cohort_df <- data.frame()
    
    temp_ctr <- splitCtr(combined_cutoff_df, ctr)
    temp_dx <- splitDisease(combined_cutoff_df, temp_ctr)
    
    cohort_df <- rbind(temp_ctr, temp_dx)
    
    
    
    
    
    if (is.function(updateProgress)) {
      text <- "Cohort Created"
      updateProgress(detail = text)
    }
    
    return(as.data.frame(cohort_df))
  }


##extracts the dx df from the merged df
splitDisease <- function(df, ctr) {
  df <- anti_join(df, ctr, by = "Sample")
  df$Cohort <- df$File[1]
  
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
  df$Cohort <- "CTR"
  
  if (nrow(df) < 0) {
    return(NULL)
  }
  return(df)
}


######



#Kruskal
######

getKruskal <- function(cohort_df, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    text <- "Performing Kruskal(This may take a while)"
    updateProgress(detail = text)
  }
  
  
  
  kruskal_results <-    cohort_df %>%
    group_by(Protein) %>%
    kruskal_test(Intensity ~ Cohort) %>%
    adjust_pvalue(method = "BH") %>%
    separate(Protein, sep = '\\|', c("sp", "Protein.ID", "Accession", "Peptide"))
  
  kruskal_results <-
    cbind(kruskal_results,
          start_seq = NA,
          end_seq = NA)
  
  if (is.function(updateProgress)) {
    text <- "Kruskal Complete"
    updateProgress(detail = text)
  }
  
  return(as.data.frame(kruskal_results))
}

######


#ANOVA
########
getANOVA <- function(cohort_df, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    text <- "Performing ANOVA (This may take a while)"
    updateProgress(detail = text)
  }
  
  
  anova_results <- cohort_df %>%
    drop_na(Intensity) %>%
    group_by(Protein) %>%
    rstatix::anova_test(Intensity ~ Cohort, detailed = T) %>%
    adjust_pvalue(method = "BH") %>%
    as_tibble() %>%
    separate(Protein, sep = '\\|', c("sp", "Protein.ID", "Accession", "Peptide")) %>%
    magrittr::set_class(c("anova_test", "rstatix_test", "data.frame"))
  
  
  
  if (is.function(updateProgress)) {
    text <- "ANOVA Complete"
    updateProgress(detail = text)
  }
  
  return(anova_results)
  
  
}



########

#TTest
########

getTtest <- function(cohort_df, p, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    text <- "Performing Ttest (This may take a while)"
    updateProgress(detail = text)
  }
  
  
  ttest_results <-  cohort_df %>%
    
    
    #Run the T-test and adjustments
    group_by(Protein) %>%
    t_test(Intensity ~ Cohort, detailed = T) %>%
    adjust_pvalue(method = "BH") %>%
    
    #Split the Protein name in Uniprot and Accession
    separate(Protein, sep = '\\|', c("sp", "Protein.ID", "Accession", "Peptide"))  %>%
    
    #Determine Fold change. Since we work with log-transformed values we can just substract
    mutate(FC = estimate1 - estimate2) %>%
    
    #Create log10 p-vals
    mutate(log10adjustP = -1 * log10(p.adj)) %>%
    
    #Determine if up or down regulated
    mutate(Direction = ifelse(p.adj > p, "NotSignificant", ifelse(FC < 0, "Down", "Up")))
  
  ##Locate the peptides in the protein sequence and add the start and end point to the dataframe
  # located_peptides() <- ttest_results[(grepl(paste(dfIds),ttest_results$UniprotID))]
  ttest_results <-
    cbind(ttest_results,
          start_seq = NA,
          end_seq = NA)
  
  if (ttest_results$group1[1] == "CTR") {
    ttest_results %<>% rename(group1 = "CTR", group2 = "Cohort")
    
  }
  
  else{
    ttest_results %<>% rename(group2 = "CTR", group1 = "Cohort")
    
  }
  
  if (is.function(updateProgress)) {
    text <- "Ttest Complete"
    updateProgress(detail = text)
  }
  
  return(as.data.frame(ttest_results))
}


#########

#MannWhitney
######


getMannWhit <- function(cohort_df, p, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    text <- "Performing MannWhitney (This may take a while)"
    updateProgress(detail = text)
  }
  
  mannwhit_results <-  cohort_df %>%
    
    
    #Run the T-test and adjustments
    group_by(Protein) %>%
    wilcox_test(Intensity ~ Cohort, detailed = T) %>%
    adjust_pvalue(method = "BH") %>%
    
    #Split the Protein name in Uniprot and Accession
    separate(Protein, sep = '\\|', c("sp", "Protein.ID", "Accession", "Peptide"))  %>%
    
    #Determine Fold change. Since we work with log-transformed values we can just substract
    mutate(FC = estimate) %>%
    
    #Create log10 p-vals
    mutate(log10adjustP = -1 * log10(p.adj)) %>%
    
    #Determine if up or down regulated
    mutate(Direction = ifelse(p.adj > p, "NotSignificant", ifelse(FC < 0, "Down", "Up")))
  
  
  
  mannwhit_results <-
    cbind(mannwhit_results,
          start_seq = NA,
          end_seq = NA)
  
  if (is.function(updateProgress)) {
    text <- "MannWhit Complete"
    updateProgress(detail = text)
  }
  
  return(as.data.frame(mannwhit_results))
}





######


#Locate Peptides
#######
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
#######



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
  
  #######
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
    tagList(radioButtons("lfqInput", "Select LFQ Type: ", c("None" = 1, "Max" = 0)))
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
    tagList(radioButtons(
      "parametricInput",
      "Select Parametric or Nonparametric: ",
      c("Parametric" = 1, "Nonparametric" = 0)
    ))
  })
  
  output$comparativeInput <- renderUI({
    "txt"
    
    req(input$confirmFasta)
    tagList(
      radioButtons(
        "comparativeInput",
        "Select Comparative TTest Results (Plotting multiple Ttest results to one plot): ",
        c("Comparative" = 1, "Noncomparative" = 0)
      )
    )
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
  
  output$pInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(numericInput(
      "pInput",
      "Select p value for signifcance : ",
      0.005,
      min = 0,
      max = 100
    ))
  })
  
  output$yInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(
      numericInput(
        "yInput",
        "Select Y-Axis (Will be used for both upregulation and downregulation) : ",
        2.5,
        min = 0,
        max = 100
      )
    )
  })
  
  
  output$confirmFile <- renderUI({
    "txt"
    req(input$confirmFasta)
    tagList(disabled(actionButton(
      "confirmFile", "Confirm File(s) Input"
    )))
  })
  
  output$ctrInput <- renderUI({
    "txt"
    req(input$confirmFasta)
    
    tagList(
      textInput(
        "ctrInput",
        "Input control identifiers seperated by commas",
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
        #"Input cohort identifiers seperated by commas (,). Comparative Ttests use parentheses () to separate different groupings. I.E. (cohort1,ctr1),(cohort2,ctr2) [These should be the file names of the combined_petide.tsv]",
        "Input cohort identifiers seperated by commas (,). [These should be the file names of the combined_petide.tsv]",
        
        placeholder = "Example AD, DLB etc."
      ),
      value = NULL
    )
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
  
  
  observeEvent(input$yInput, {
    y_axis(input$yInput)
    
  })
  
  observeEvent(input$pInput, {
    p_cutoff(input$pInput)
    
  })
  
  observeEvent(input$comparativeInput, {
    comparative(input$comparativeInput)
    
  })
  
  observeEvent(input$ctrInput, {
    nameTemp <- input$ctrInput
    # nameTemp <- strsplit(nameTemp, ", |,| , ")
    ctr_name(unlist(nameTemp))
    
    # if (!is.null(input$fileInput))
    #   shinyjs::enable("confirmFile")
    
    
  })
  
  observeEvent(input$cohortInput, {
    nameTemp <- input$cohortInput
    nameTemp <- strsplit(nameTemp, ", |,| , ")
    cohort_name(unlist(nameTemp))
    
    if (!is.null(input$fileInput))
      shinyjs::enable("confirmFile")
    
    
  })
  
  
  
  ########
  
  
  
  
  
  
  observeEvent(input$confirmFile, {
    #UI LOGIC
    #######
    shinyjs::hide("fileInput")
    shinyjs::hide("confirmFile")
    shinyjs::hide("lfqInput")
    shinyjs::hide("cutoffInput")
    shinyjs::hide("intensityInput")
    shinyjs::hide("cohortInput")
    shinyjs::hide("ctrInput")
    shinyjs::hide("parametricInput")
    shinyjs::hide("pInput")
    shinyjs::hide("yInput")
    shinyjs::hide("comparativeInput")
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
    uniprot_ids <- as.data.frame(uniprot())
    protein_plot <- list()
    colnames(uniprot_ids)  <- c("Protein.ID")
    cohorts <- cohort_name()
    parametric_type <- parametric()
    ctr <- ctr_name()
    protein_id_loop <- list
    
    
    fasta <- as.data.frame(fasta_reactive())
    fasta <- tibble::rownames_to_column(fasta, "accession")
    fasta %<>%
      separate(accession, sep = '\\|', c("sp", "Protein.ID", "info")) %>%
      separate(info, sep = 'OX=[A-Za-z0-9]+', c("info", "Gene"))
    
    df <- files()
    #######
    
    #Dateframe creation
    #######
    if (comparative() == 0) {
      for (i in df) {
        temp <-
          cleanData(
            i,
            as.integer(intensity()),
            as.integer(lfq()),
            as.integer(cutoff()),
            updateProgress
          )
        
        temp <- cohortSplit(temp, ctr, updateProgress)
        
        combined_cutoff_df <- rbind(combined_cutoff_df, temp)
      }
      
      if (length(cohorts) < 2) {
        combined_cutoff_df_possible <-
          combined_cutoff_df %>%  spread(Cohort, Intensity) %>% group_by(Protein) %>% summarise(ctr = sum(!is.na(get("CTR"))),
                                                                                                cohort1 = sum(!is.na(get(
                                                                                                  .env$cohorts[1]
                                                                                                )))) %>%
          mutate(possible = ifelse(ctr < 2 |
                                     cohort1 < 2, FALSE, TRUE)) %>%
          filter(possible)
      }
      else{
        combined_cutoff_df_possible <-
          combined_cutoff_df %>%  spread(Cohort, Intensity) %>% group_by(Protein) %>% summarise(
            ctr = sum(!is.na(get(.env$ctr))),
            cohort1 = sum(!is.na(get(
              .env$cohorts[1]
            ))),
            cohort2 = sum(!is.na(get(
              .env$cohorts[2]
            )))
          ) %>%
          mutate(possible = ifelse(ctr < 2 |
                                     cohort1 < 2 |
                                     cohort2 < 2, FALSE, TRUE)) %>%
          filter(possible)
      }
      
      cohort_df <-
        combined_cutoff_df %>% filter(Protein %in% combined_cutoff_df_possible$Protein)
      
      
      if (nrow(cohort_df) < 1) {
        session$close()
        
      }
      ########
      
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
    
    #Comparative
    #########
    browser()
    for (i in df) {
      temp <-
        cleanData(
          i,
          as.integer(intensity()),
          as.integer(lfq()),
          as.integer(cutoff()),
          updateProgress
        )
      
      temp <- cohortSplit(temp, ctr, updateProgress)
      
      
      for (i in cohorts) {
        if (i == temp$File[1]) {
          combined_cutoff_df_possible <-
            temp %>%  spread(Cohort, Intensity) %>% group_by(Protein) %>% summarise(ctr = sum(!is.na(get("CTR"))),
                                                                                    cohort1 = sum(!is.na(get(.env$i)))) %>%
            mutate(possible = ifelse(ctr < 2 |
                                       cohort1 < 2, FALSE, TRUE)) %>%
            filter(possible)
        }
      }
      
      
      cohort_df <-
        temp %>% filter(Protein %in% combined_cutoff_df_possible$Protein)
      
      
      
      
      if (parametric_type == 0) {
        statisical_test <-
          getMannWhit(cohort_df, p_cutoff(), updateProgress)
        
        
        
      }
      else{
        statisical_test <- getTtest(cohort_df, p_cutoff(), updateProgress)
        
      }
      
      comparative_combined <-
        rbind(comparative_combined, statisical_test)
    }
    
    
    
    
    
    ######
    
    
    if (comparative() == 0) {
      leftovers <-
        anti_join(statisical_test, fasta, by = "Protein.ID")
      
      statisical_test %<>% inner_join(fasta, by = "Protein.ID", multiple = "all") %>% arrange(Protein.ID)
      
      located_peptides(locatePeptides(statisical_test, updateProgress))
      # browser()
      
      #Plotting loop
      #######
      if (nrow(uniprot_ids) > 0) {
        protein_id_loop <- uniprot_ids
        
      }
      else{
        protein_id_loop <-  unique(comparative_combined$Protein.ID)
      }
      for (i in protein_id_loop) {
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
          ylim(-y_axis(), y_axis()) +
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
          
          #Specify the colours I want to use for the isSignificant column
          scale_fill_manual(values = c("yes" = "red", "no" = "grey")) +
          theme(legend.position = "bottom") +
          
          #x and yaxis titles
          xlab("Protein Sequence") +
          ylab("FC") +
          labs(subtitle = as.character(i),
               title = as.character(filtered_results()$Gene),) +
          
          annotate(
            "text",
            x = 15 ,
            y = y_axis(),
            label = paste0("p == ", p_cutoff()),
            parse = TRUE
          )
        
        
        
        
      }
    }
    ######
    else{
      leftovers <-
        anti_join(comparative_combined, fasta, by = "Protein.ID")
      
      comparative_combined %<>% inner_join(fasta, by = "Protein.ID", multiple = "all") %>% arrange(Protein.ID) %>% mutate(
        color = case_when(
          Cohort == as.character(cohorts[1]) ~  "red",
          Cohort == as.character(cohorts[2]) ~ "orange",
          Cohort == as.character(cohorts[3]) ~ "green"
        )
      )
      
      # cols <- c(as.character(cohorts[1]) = "red", as.character(cohorts[2]) = "blue", as.character(cohorts[1]) = "green")
      
      located_peptides(locatePeptides(comparative_combined, updateProgress))
      # browser()
      colors <- distinct(comparative_combined, Cohort, color)
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
        
        
        protein_plot[[i]] = plotting_protein() %>%
          
          #Lets make a column based on significance
          mutate(isSignificant = ifelse(filtered_results()$p.adj < 0.05, "yes", "no")) %>%
          
          #Start plotting
          ggplot() +
          
          #We take here the 'mean' but this is of course X-times the same value
          ylim(-y_axis(), y_axis()) +
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
              fill = Cohort,
            ),
            col = "black",
            alpha = 0.75,
          ) +
          
          #Set the Ggplot theme, limit the y-axis for FC.
          theme_bw() +
          
          #Specify the colours I want to use for the isSignificant column
          scale_fill_manual(values = pal()) +
          
          theme(legend.position = "bottom") +
          
          #x and yaxis titles
          xlab("Protein Sequence") +
          ylab("FC") +
          
          labs(subtitle = as.character(i),
               title = as.character(filtered_results()$Gene),) +
          
          annotate(
            "text",
            x = 10 ,
            y = y_axis(),
            label = paste0("p == ", p_cutoff()),
            parse = TRUE
          )
        
        
        
        
      }
      
      
      
    }
    
    #Reactive Variables
    ######
    protein_plotr(protein_plot)
    plot_fasta(fasta)
    if (nrow(uniprot_ids) == 0) {
      if (comparative() == 0) {
        temp_id <- as.data.frame(unique(statisical_test$Protein.ID))
      }
      else{
        temp_id <- as.data.frame(unique(comparative_combined$Protein.ID))
        
      }
      colnames(temp_id) <- c("Protein.ID")
      protein_ids(temp_id)
    }
    else{
      protein_ids(uniprot_ids)
    }
    final_dataframe(statisical_test)
    removed_proteins(leftovers)
    ######
    
    
    shinyjs::enable("clickThrough")
    shinyjs::enable("downloadPlot")
    
    
    
    
    
    
    
    
    
    
  })
  
  output$downloadPlot <- downloadHandler(
    filename = 'plots.zip',
    content = function(file) {
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
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
      
      write.csv(final_dataframe()
                , file = "final_dataframe.csv"
                , row.names = F)
      write.csv(removed_proteins()
                ,
                file = "removed_proteins.csv"
                ,
                row.names = F)
      
      zip_files <-
        list.files(path = getwd(), pattern = ".svg$|.csv$|.png$")
      
      zip::zip(file, files = zip_files)
      
      
    }
    
  )
  
  
  
  
  observeEvent(input$clickThrough, {
    plots <- as.list(protein_plotr())
    ids <- as.data.frame(protein_ids())
    if (j() < nrow(ids)) {
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