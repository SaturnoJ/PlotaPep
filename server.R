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

# clean the data
clean_combined_peptide_fragger <-
  function(dataset, Intensity, MaxLFQ, cutoff) {
    #change Protein name
    names(dataset)[names(dataset) == 'Protein ID'] <- 'ProteinID'
    names(dataset)[names(dataset) == 'Peptide Sequence'] <-
      'Sequence'
    
    #Select columns of interest
    dataset <-
      dplyr::mutate(dataset, Protein = paste0(ProteinID, "_", Gene, "_", Sequence))
    
    # Select the specified Intensity.
    # We also start creating the substring we will remove from sample names.
    if (Intensity == 0) {
      dataset <-
        dplyr::select(
          dataset,
          Protein,
          contains("Intensity"),-contains("Total"),-contains("Unique")
        )
      removeSubString <- paste0(" ", "Intensity")
    } else if (Intensity == 1) {
      dataset <-
        dplyr::select(dataset, Protein, contains("Total Intensity"))
      removeSubString <- paste0(" ", "Total", " Intensity")
    } else {
      dataset <-
        dplyr::select(dataset, Protein, contains("Unique Intensity"))
      removeSubString <- paste0(" ", "Unique", " Intensity")
    }
    
    #select or not select MaxLFQ
    if (MaxLFQ == 0) {
      dataset <- dplyr::select(dataset, Protein, contains("MaxLFQ"))
      removeSubString <- paste0(" MaxLFQ", removeSubString)
    } else {
      dataset <-
        dplyr::select(dataset,
                      Protein,
                      contains("Intensity"),-contains("MaxLFQ"))
    }
    
    #col to rowname
    dataset <- tibble::column_to_rownames(dataset, "Protein")
    
    #Get rid of substring part we do not need
    colnames(dataset) <-
      stringr::str_remove(colnames(dataset), removeSubString)
    
    #transpose, remove zeros, log2
    dataset <- data.frame(t(dataset))
    dataset[dataset == 0] <- NA
    dataset <- log2(dataset)
    
    #apply cutoff on the protein level
    
    dataset <-
      purrr::discard(dataset, ~ sum(is.na(.x)) / length(.x) * 100 >= (100 - cutoff))
    
    
    #Long Format it
    dataset <- tibble::rownames_to_column(dataset, "Sample")
    dataset <-
      tidyr::gather(dataset,
                    colnames(dataset)[2:ncol(dataset)],
                    key = "Protein",
                    value = "Intensity")
    
    return(dataset)
  }

##Remove the outliers
OutlierRemover <- function(dataset) {
  #Determine correlations between sample and theoretical
  SampleCorrelations <- dataset %>%
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
  returnDataset <- dataset %>%
    dplyr::filter(Sample %in% dplyr::filter(SampleCorrelations, Outlier == F)$Sample) %>%
    dplyr::select(Protein, Sample, Intensity)
  
  #Half minimum value (per protein) imputation
  PCA_DF <- dataset %>%
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
  
  return(returnDataset)
}


##extracts the alzheimers dataset from the merged dataset
split_disease <- function(df, ctr) {
  df <- subset(df, !grepl("CTR", df$Sample))
  if (isEmpty(df)) {
    return(NULL)
  }
  return(df)
  
}

##extracts the control dataset from the merged dataset
split_ctr <- function(combined_cutoff_df) {
  df <-
    subset(combined_cutoff_df,
           grepl("CTR", combined_cutoff_df$Sample))
  if (isEmpty(df)) {
    return(NULL)
  }
  return(df)
}
# Loads already found semi-trypic values into aggregated dataset.
load_semi <- function(cleavages, df) {
  for (i in 1:nrow(cleavages)) {
    for (j in 1:nrow(df)) {
      if (cleavages$Sequence[i] == df$Sequence[j]) {
        cleavages$N_terminus[i] <- df$N_terminus[j]
        cleavages$C_terminus[i] <- df$C_terminus[j]
        break
      }
    }
  }
  return(cleavages)
}
# Loads protein symbols into aggregate dataset
load_protein <- function(cleavages, df) {
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(cleavages)) {
      if (df$Sequence[i] == cleavages$Sequence[j] &&
          !is.null(cleavages$protein[j])) {
        cleavages$protein[j] <- df$Accession[i]
      }
    }
  }
  return(cleavages)
}
# Loads where missed cleavage location occurred into aggregate dataset
cut_at <- function(cleavages) {
  for (i in 1:nrow(cleavages)) {
    if (cleavages$N_terminus[i] == "False") {
      cleavages$cleavage_loc[i] <- cleavages$start_seq[i]
    } else if (cleavages$C_terminus[i] == "False") {
      cleavages$cleavage_loc[i] <- cleavages$end_seq[i]
    }
  }
  
  return(cleavages)
}
# Filters out proteins above cutoff and not being searched for
filterCutSym <- function(df, cutoff, symbols) {
  df <- filter(df, cutoff > df$Expect)
  x <- NULL
  for (i in 1:nrow(symbols)) {
    y <- filter(df,
                grepl(symbols$UNIPROT[i], df$Accession))
    x <- rbind(x, y)
  }
  df <- x
  
  return(df)
}
# Adds symbols to Accession column in dataframe
add_symbols <- function(df, symbols) {
  for (i in 1:nrow(symbols)) {
    for (j in 1:nrow(df)) {
      if (grepl(symbols$UNIPROT[i], df$Accession[j])) {
        df$Accession[j] <- symbols$SYMBOL[i]
      }
    }
  }
  return(df)
}

locate <- function(cleavages, fasta) {
  for (i in 1:nrow(fasta)) {
    for (j in 1:nrow(cleavages)) {
      if (grepl(cleavages$protein[j], fasta$accession[i])) {
        locate <-
          as.data.frame(str_locate(fasta$seq[i], cleavages$Sequence[j]))
        if (!is.na(locate)) {
          cleavages$start_seq[j] <- locate$start
          cleavages$end_seq[j] <- locate$end
        }
      }
    }
  }
  return(cleavages)
}

remove_false <- function(df) {
  return(subset(df, !grepl("XXX", df$Accession)))
}
quant <- function(proteins, x) {
  if (missing(x)) {
    return(subset(proteins, Freq > quantile(proteins$Freq, probs = 0.95)))
  } else {
    return(subset(proteins, Freq > quantile(proteins$Freq, probs = x)))
  }
}

top_n <- function(df, x) {
  df <- df[order(-df$Freq), ]
  if (missing(x)) {
    return(df <- head(df, -(nrow(df) - 100)))
  } else {
    return(df <- head(df, -(nrow(df) - x)))
  }
}


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024 ^ 3)
  
  fastaR <- reactiveVal()
  files <- reactiveVal()
  uniprot <- reactiveVal()
  downloadData <- reactiveVal()
  df <- data.frame()
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
    tagList(radioButtons("lfqInput", "Select LFQ Type: ", c(
      "Max: " = 0, "None: " = 1
    )))
  })
  
  output$intensityInput <- renderUI({
    "txt"
    
    req(input$confirmFasta)
    tagList(radioButtons(
      "intensityInput",
      "Select Intensity : ",
      c(
        "Intensity: " = 0,
        "Total: " = 1,
        "Unique: " = 2
      )
    ))
  })
  
  output$cutoffInput <- renderUI({
    "txt"
    req(input$confirmFasta)    
    tagList(numericInput(
      "cutoffInput",
      "Select Cutoff : ",
      95,
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
  
  output$uniprotInput <- renderUI({
    "txt"
    req(input$fastaInput)
    
    tagList(
      textInput(
        "uniprotInput",
        "Input Uniprot IDs seperated by commas (,)",
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
      if (!is.null(uniprot)) {
        shinyjs::enable("confirmFile")
        
      }
      
    })
    shinyjs::enable("confirmFasta")
    
  })
  observeEvent(input$confirmFasta, {
    shinyjs::hide("fastaInput")
    shinyjs::hide("confirmFasta")
    
  })
  
  
  
  observeEvent(input$fileInput, {
    files <- input$fileInput
    for (i in 1:nrow(files)) {
      if (i == 1) {
        df <- read_delim(files$datapath[i], show_col_types = FALSE)
        
        df$fileName <- files$name[i]
        
      }
      else{
        df1 <- read_delim(files$datapath[i],
                          show_col_types = FALSE)
        df1$fileName <- files$name[i]
        
        df <- rbind(df, df1)
        
        
      }
      
    }
    files(df)
    
    
  })
  
  observeEvent(input$confirmFile, {
    shinyjs::hide("fileInput")
    shinyjs::hide("confirmFile")
    shinyjs::hide("uniprotInput")
    fasta <- as.data.frame(fastaR())
    df <- as.data.frame(files())
    
    cleaned_df <-
      clean_combined_peptide_fragger(df, 0, 1, 95)
    
    outlier_removed_df <- OutlierRemover(cleaned_df)
    
    
    combined_cutoff_df <- outlier_removed_df %>%
      mutate(uniqueSamplesInDF = length(unique(Sample))) %>%
      group_by(Protein) %>%
      mutate(countMissingValues = sum(!is.na(Intensity))) %>%
      mutate(percentageMissing = countMissingValues / uniqueSamplesInDF) %>%
      filter(percentageMissing > 0.95) %>% # Apply cutoff of 95%
      dplyr::select(Sample, Protein, Intensity) #clean up by selecting the original columns
    combined_cutoff_df <- as.data.frame(combined_cutoff_df)
    #######
    
    control_df <- split_ctr(combined_cutoff_df)
    control_df <- cbind(control_df, "CTR")
    colnames(control_df) <-
      c('Sample', 'Protein', 'Intensity', 'DX')
    alzheimers_df <-
      split_disease(combined_cutoff_df, control_df)
    alzheimers_df <- cbind(alzheimers_df, "AD")
    colnames(alzheimers_df) <-
      c('Sample', 'Protein', 'Intensity', 'DX')
    cohort_df <- rbind(control_df, alzheimers_df)
    cohort_df$Protein <- gsub('\\.', 'A', cohort_df$Protein)
    
    ## change this shit
    for (i in 1:nrow(cohort_df)) {
      fn <- "_"
      rp <- "A"
      n <- 2
      a <- str_count(cohort_df$Protein[i], "_")
      if (a == 3) {
        regmatches(cohort_df$Protein[i], gregexpr(fn, cohort_df$Protein[i])) <-
          list(c(rep(fn, n - 1), rp))
      }
      i <- i - 1
    }
    
    
    ##Apply a peptide level cutoff of 95%
    ttest_results <- cohort_df %>%
      
      
      #Run the T-test and adjustments
      group_by(Protein) %>%
      t_test(Intensity ~ DX, detailed = T) %>%
      adjust_pvalue(method = "BH") %>%
      
      #Split the Protein name in Uniprot and Gene
      separate(Protein, c("UniprotID", "Gene", "Peptide")) %>%
      
      #Determine Fold change. Since we work with log-transformed values we can just substract
      mutate(FC = estimate1 - estimate2) %>%
      
      #Create log10 p-vals
      mutate(log10adjustP = -1 * log10(p.adj)) %>%
      
      #Determine if up or down regulated
      mutate(Direction = ifelse(p.adj > 0.05, "NotSignificant", ifelse(FC < 0, "Down", "Up")))
    
    ##Locate the peptides in the protein sequence and add the start and end point to the dataframe
    located_peptides <- ttest_results %>%
      filter(UniprotID == protein)
    located_peptides <-
      cbind(located_peptides,
            start_seq = NA,
            end_seq = NA)
    locate_peptides <- function(cleavages, fasta) {
      for (i in 1:nrow(cleavages)) {
        locate <-
          as.data.frame(str_locate(toString(fasta$seq), toString(cleavages$Peptide[i])))
        cleavages$start_seq[i] <- locate$start
        cleavages$end_seq[i] <- locate$end
      }
      return(cleavages)
    }
    located_peptides <-
      locate_peptides(located_peptides, fasta)
    
    
    ##Plot the peptides within the protein
    #We need to select the protein we want first
    proteinWeWantToPlot <- ttest_results %>%
      filter(UniprotID == protein)
    filtered_results <- ttest_results %>%
      filter(UniprotID == protein)
    
    protein_length <- nchar(fasta$seq)
    
    #Now we can make the plot
    proteinWeWantToPlot %>%
      
      #Lets make a column based on significance
      mutate(isSignificant = ifelse(filtered_results$p.adj < 0.05, "yes", "no")) %>%
      
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
        aes(
          xmin = located_peptides$start_seq,
          xmax = (
            located_peptides$start_seq + (located_peptides$end_seq - located_peptides$start_seq)
          ),
          ymin = filtered_results$FC - 0.05,
          ymax = filtered_results$FC + 0.05,
          fill = isSignificant
        ),
        
        #Here I specify some stuff that will be universal for all blocks irrespective of column info.
        col = "black",
        alpha = 0.75
      ) +
      
      #Set the Ggplot theme, limit the y-axis for FC.
      theme_bw() +
      ylim(-2, 2) +
      
      #Specify the colours I want to use for the isSignificant column
      scale_fill_manual(values = c("yes" = "red", "no" = "grey")) +
      theme(legend.position = "bottom") +
      
      #x and yaxis titles
      xlab("Protein Sequence") +
      ylab("FC")
    
    shinyjs::enable("download")
    
  })
  
  
  output$download <- downloadHandler(
    filename = function() {
      "processed_data.csv"
    },
    content = function(fname) {
      write.csv(downloadData(), fname)
    }
  )
  
  
  
  
  
  
  
}
