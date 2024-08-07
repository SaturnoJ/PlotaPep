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


cleanData <-
  function(df,
           intensity,
           lfq,
           file_type,
           updateProgress = NULL) {
    if (is.function(updateProgress)) {
      text <- "Cleaning Data"
      updateProgress(detail = text)
    }
    
    file <- df$File[1]
    df <-  subset(df, select=-c(File))
    if(file_type == 2){
      cleaned_df <-
        diannCleaner(df)
      
      }
    else if(file_type == 1){
      cleaned_df <-
        fraggerCleanerIon(df, intensity, lfq, )
    }
    else if(file_type == 3){
      cleaned_df <-
        genericCleaner(df)
    }
    else{
      cleaned_df <-
        fraggerCleaner(df, intensity, lfq)
    }
    if (is.function(updateProgress)) {
      text <- "Data Cleaned"
      updateProgress(detail = text)
    }
    
    
    cleaned_df$File <- file
    
    return(as.data.frame(cleaned_df))
  }

diannCleaner <-function(df ){
  
  #combine the ProteinGroup and Genes
 df <-
      unite(df,
            Protein,
            c(Protein.Ids,Genes, Stripped.Sequence, Modified.Sequence, Precursor.Charge),
            sep = "|",
            remove = FALSE)  
  #Proteins to rownames
  df <- tibble::column_to_rownames(df, "Protein")
  
  #deselect the not-user stuff
  df <- dplyr::select(df, - Protein.Group, - Protein.Names, - Proteotypic, - Stripped.Sequence, - Precursor.Charge, - Modified.Sequence, - First.Protein.Description, - Protein.Group, - Genes, - Precursor.Id, - Protein.Ids)
  
  #Clean the Sample Names
  colnames(df) <- sapply(strsplit(colnames(df), "\\", fixed=TRUE), tail, 1)
  colnames(df) <- stringr::str_remove(colnames(df), ".d")
  
  #transpose, remove zeros, log2
  df_transposed <-  data.frame(t(df))
  rownames(df_transposed) <- colnames(df)
  colnames(df_transposed) <- rownames(df)
  df <- df_transposed
   
  df <- removeZero(df)
  
  
  
  #apply cutoff on the protein level

  
  #Long Format it
  df <- setDT(df, keep.rownames = "Sample")
  df <- tidyr::gather(df,
                      colnames(df)[2:ncol(df)],
                      key = "Protein",
                      value = "Intensity")
  df$Intensity <- as.numeric(df$Intensity)
  df <- removeZero(df)
  
    return(df)
  
  
}
  

fraggerCleanerIon <-
  function(df, intensity, lfq) {
    # change Protein name
    if (grep('Peptide Sequence', names(df))) {
      names(df)[names(df) == 'Peptide Sequence'] <-
        'Sequence'
    }
    if (grep('Modified Sequence', names(df))) {
      names(df)[names(df) == 'Modified Sequence'] <-
        'Modified'
    }
    if (grep('Protein ID', names(df))) {
      names(df)[names(df) == 'Protein ID'] <-
        'Protein.Ids'
    }


    #Select columns of interest
    df <-
      unite(df,
            Protein,
            c(Protein.Ids,Gene, Sequence, Modified, Charge),
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
    df_transposed <-  data.frame(t(df))
    rownames(df_transposed) <- colnames(df)
    colnames(df_transposed) <- rownames(df)
    df <- df_transposed
    
    df <- removeZero(df)
    
    
    
    #apply cutoff on the protein level
   
    
    #Long Format it
    df <- setDT(df, keep.rownames = "Sample")
    df <- tidyr::gather(df,
                        colnames(df)[2:ncol(df)],
                        key = "Protein",
                        value = "Intensity")
    
    return(df)
  }

fraggerCleaner <-
  function(df, intensity, lfq) {
    # change Protein name
    if (grep('Peptide Sequence', names(df))) {
      names(df)[names(df) == 'Peptide Sequence'] <-
        'Sequence'
    }
    if (grep('Protein ID', names(df))) {
      names(df)[names(df) == 'Protein ID'] <-
        'Protein.Ids'
    }
    
    #Select columns of interest
    df <-
      unite(df,
            Protein,
            c(Protein.Ids,Gene, Sequence),
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
    df_transposed <-  data.frame(t(df))
    rownames(df_transposed) <- colnames(df)
    colnames(df_transposed) <- rownames(df)
    df <- df_transposed
    
    df <- removeZero(df)
    
    
    
    
    
    #Long Format it
    df <- setDT(df, keep.rownames = "Sample")
    df <- tidyr::gather(df,
                        colnames(df)[2:ncol(df)],
                        key = "Protein",
                        value = "Intensity")
    
    return(df)
  }


genericCleaner <-function(df){

  #combine the ProteinGroup and Genes
  df <-
    unite(df,
          Protein,
          c(Protein.Ids,Genes, Sequence, Unique),
          sep = "|",
          remove = TRUE)  
  #Proteins to rownames
  df <- tibble::column_to_rownames(df, "Protein")
  
  #transpose, remove zeros, log2
  df_transposed <-  data.frame(t(df))
  rownames(df_transposed) <- colnames(df)
  colnames(df_transposed) <- rownames(df)
  df <- df_transposed
  
  df <- removeZero(df)
  
  
  
  #apply cutoff on the protein level
 
  
  #Long Format it
  df <- setDT(df, keep.rownames = "Sample")
  df <- tidyr::gather(df,
                      colnames(df)[2:ncol(df)],
                      key = "Protein",
                      value = "Intensity")
  df$Protein <- as.character(df$Protein)
  df$Intensity <- as.numeric(df$Intensity)
  df <- removeZero(df)
  
  return(df)
  
  
}

removeZero <- function(df) {
  for (j in seq_len(ncol(df))) {
    set(df, which((df[[j]]) == 0), j, NA)
  }
  return(df)
}


cohortSplit <-
  function(combined_cutoff_df,
           ctr,
           cohort,
           std,
           updateProgress = NULL) {
    if (is.function(updateProgress)) {
      text <- "Creating Cohorts"
      updateProgress(detail = text)
    }
    cohort_df <- data.frame()
    temp_ctr <- splitCtr(combined_cutoff_df, ctr)
    temp_ctr <- removeOutlier(temp_ctr, std,updateProgress)
    temp_dx <- splitDisease(combined_cutoff_df, cohort)
    temp_dx <- removeOutlier(temp_dx, std,updateProgress)
    
    cohort_df <- rbind(temp_ctr, temp_dx)
    
    
    
    
    
    if (is.function(updateProgress)) {
      text <- "Cohort Created"
      updateProgress(detail = text)
    }
    
    return(as.data.frame(cohort_df))
  }
removeOutlier <- function(df, std,updateProgress = NULL) {
  if (is.function(updateProgress)) {
    text <- "Removing Outliers"
    updateProgress(detail = text)
  }
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
  
  sample_correlations_outlier <- sample_correlations %>%
    dplyr::mutate(Outlier = SampleCor < (1 - std * sd(
      unique(sample_correlations$SampleCor), na.rm = T
    )))
  
  #Make the DF that we will return (i.e. without outliers)
  returndf <- df %>%
    dplyr::filter(Sample %in% dplyr::filter(sample_correlations_outlier, Outlier == F)$Sample) %>%
    dplyr::select(Protein, Sample, Intensity, File, Cohort)
  if (is.function(updateProgress)) {
    text <- "Removed Outliers"
    updateProgress(detail = text)
  }
  return(returndf)
}

splitDisease <- function(df, cohort) {
  df <-
    subset(df,
           grepl(paste(cohort, collapse = "|"), df$Sample))
  
  df$Cohort <- cohort
  
  if (nrow(df) < 0) {
    return(NULL)
  }
  return(df)
  
}

splitCtr <- function(df, cohort) {
  df <-
    subset(df,
           grepl(paste(cohort, collapse = "|"), df$Sample))
  
  df$Cohort <- "CTR"
  
  if (nrow(df) < 0) {
    return(NULL)
  }
  return(df)
}






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
    separate(Protein, sep = '\\|', c("Protein.ID", "Gene", "Peptide"))   %>%
    
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
  ttest_results$File <- unique(cohort_df$File)

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
    separate(Protein, sep = '\\|', c("Protein.ID", "Gene", "Peptide"))  %>%
    
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
  
  
  if (mannwhit_results$group1[1] == "CTR") {
    mannwhit_results %<>% rename(group1 = "CTR", group2 = "Cohort")
    
  }
  
  else{
    mannwhit_results %<>% rename(group2 = "CTR", group1 = "Cohort")
    
  }
  
  if (is.function(updateProgress)) {
    text <- "MannWhit Complete"
    updateProgress(detail = text)
  }
  
  return(as.data.frame(mannwhit_results))
}






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


