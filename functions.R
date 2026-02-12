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
library(ComplexHeatmap)


#Top level function of cleaning data. Determines what file type was inputted and
#then chooses the proper cleaner based on that. Returns cleaned file.
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
    df <-  subset(df, select = -c(File))
    if (file_type == 2) {
      cleaned_df <-
        diannCleaner(df)

    }
    else if (file_type == 1) {
      cleaned_df <-
        fraggerCleanerIon(df, intensity, lfq)
    }
    else if (file_type == 3) {
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

#The next four functions perform essentially the same processes but for different
#types of input dataframes. At its core they are combining specific columns that
#have pertinent information and removing the others. They are then transposed and
#have their zero values removed. The final dataframe is returned in a long format.

diannCleaner <- function(df) {
  #combine the ProteinGroup and Genes
  df <-
    unite(
      df,
      Protein,
      c(
        Protein.Ids,
        Genes,
        Stripped.Sequence,
        Modified.Sequence,
        Precursor.Charge
      ),
      sep = "|",
      remove = FALSE
    )
  #Proteins to rownames
  df <- tibble::column_to_rownames(df, "Protein")

  #deselect the not-user stuff
  df <- dplyr::select(
    df,
    -Protein.Group,
    -Protein.Names,
    -Proteotypic,
    -Stripped.Sequence,
    -Precursor.Charge,
    -Modified.Sequence,
    -First.Protein.Description,
    -Protein.Group,
    -Genes,
    -Precursor.Id,
    -Protein.Ids
  )

  #Clean the Sample Names
  colnames(df) <- sapply(strsplit(colnames(df), "\\", fixed = TRUE), tail, 1)
  colnames(df) <- stringr::str_remove(colnames(df), ".d")

  #transpose, remove zeros
  df_transposed <-  data.frame(t(df))
  rownames(df_transposed) <- colnames(df)
  colnames(df_transposed) <- rownames(df)
  df <- df_transposed

  df <- removeZero(df)


  #Long Format it
  df <- setDT(df, keep.rownames = "Sample")
  df <- tidyr::gather(df, colnames(df)[2:ncol(df)], key = "Protein", value = "Intensity")
  df$Intensity <- as.numeric(df$Intensity)
  return(df)


}


fraggerCleanerIon <-
  function(df, intensity, lfq) {
    # change Protein name
    if (grep('Peptide Sequence', names(df))| grep('Peptide.Sequence', names(df))) {
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
      unite(
        df,
        Protein,
        c(Protein.Ids, Gene, Sequence, Modified),
        sep = "|",
        remove = FALSE
      )
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
        dplyr::select(df, Protein, contains("Intensity"), -contains("MaxLFQ"))
    }
    #col to rowname
    df <- data.frame(df[, -1], row.names = df[, 1])


    #Get rid of substring part we do not need

    colnames(df) <-
      stringr::str_remove(colnames(df), removeSubString)

    #transpose, remove zeros
    df_transposed <-  data.frame(t(df))
    rownames(df_transposed) <- colnames(df)
    colnames(df_transposed) <- rownames(df)
    df <- df_transposed

    df <- removeZero(df)

    #Long Format it
    df <- setDT(df, keep.rownames = "Sample")
    df <- tidyr::gather(df, colnames(df)[2:ncol(df)], key = "Protein", value = "Intensity")

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
            c(Protein.Ids, Gene, Sequence),
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
        dplyr::select(df, Protein, contains("Intensity"), -contains("MaxLFQ"))
    }

    #col to rowname
    df <- data.frame(df[, -1], row.names = df[, 1])


    #Get rid of substring part we do not need

    colnames(df) <-
      stringr::str_remove(colnames(df), removeSubString)

    #transpose, remove zeros
    df_transposed <-  data.frame(t(df))
    rownames(df_transposed) <- colnames(df)
    colnames(df_transposed) <- rownames(df)
    df <- df_transposed

    df <- removeZero(df)

    #Long Format it
    df <- setDT(df, keep.rownames = "Sample")
    df <- tidyr::gather(df, colnames(df)[2:ncol(df)], key = "Protein", value = "Intensity")

    return(df)
  }


genericCleaner <- function(df) {


  #combine the ProteinGroup and Genes
  df <-
    unite(
      df,
      Protein,
      c(Protein.Ids, Genes, Sequence, Unique),
      sep = "|",
      remove = TRUE
    )
  #Proteins to rownames
  df <- tibble::column_to_rownames(df, "Protein")

  #transpose, remove zeros
  df_transposed <-  data.frame(t(df))
  rownames(df_transposed) <- colnames(df)
  colnames(df_transposed) <- rownames(df)
  df <- df_transposed

  df <- removeZero(df)

  #Long Format it
  df <- setDT(df, keep.rownames = "Sample")
  df <- tidyr::gather(df, colnames(df)[2:ncol(df)], key = "Protein", value = "Intensity")
  df$Protein <- as.character(df$Protein)
  df$Intensity <- as.numeric(df$Intensity)
  return(df)


}

#Removes all zeros from the dataframe to allow for log2 to be performed at a latter step.
#Returns same dataframe but with all zeros replaced by NA
removeZero <- function(df) {
  for (j in seq_len(ncol(df))) {
    set(df, which((df[[j]]) == 0), j, NA)
  }
  return(df)
}

#Top level function that calls dataframe splitting functions and outlier removal functions.
#After both case and control have been split and outlier removed the two dataframes
#are merged back into a single dataframe with all outliers removed. Returns inputted
#dataframe with outliers removed.
caseSplit <-
  function(combined_cutoff_df,
           ctr,
           case,
           std,
           updateProgress = NULL) {
    if (is.function(updateProgress)) {
      text <- "Creating Cases"
      updateProgress(detail = text)
    }
    case_df <- data.frame()
    temp_ctr <- splitCtr(combined_cutoff_df, ctr)
    temp_ctr <- removeOutlier(temp_ctr, std, updateProgress)
    temp_dx <- splitDisease(combined_cutoff_df, case)
    temp_dx <- removeOutlier(temp_dx, std, updateProgress)

    case_df <- rbind(temp_ctr, temp_dx)
    if (is.function(updateProgress)) {
      text <- "Case Created"
      updateProgress(detail = text)
    }

    return(as.data.frame(case_df))


  }
#Outliers are removed by this function. The control and case dataframe are entered
#at a later step. Returns dataframe with outliers removed.
removeOutlier <- function(df, std, updateProgress = NULL) {
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
    dplyr::select(Protein, Sample, Intensity, File, Case)
  if (is.function(updateProgress)) {
    text <- "Removed Outliers"
    updateProgress(detail = text)
  }
  return(returndf)
}

#Splits dataframe into a disease only dataframe for outlier removal. Returns
#dataframe with only case values.

splitDisease <- function(df, case) {
  df <-
    subset(df, grepl(paste(case, collapse = "|"), df$Sample))

  df$Case <- case

  if (nrow(df) < 0) {
    return(NULL)
  }
  return(df)

}
#Splits the dataframe into a control only dataframe for outlier removal. Returns
#dataframe with only control values/

splitCtr <- function(df, case) {
  df <-
    subset(df, grepl(paste(case, collapse = "|"), df$Sample))

  df$Case <- "CTR"

  if (nrow(df) < 0) {
    return(NULL)
  }
  return(df)
}


#
# Keeping functons for possible later implementation
#
#
# getKruskal <- function(case_df, updateProgress = NULL) {
#   if (is.function(updateProgress)) {
#     text <- "Performing Kruskal(This may take a while)"
#     updateProgress(detail = text)
#   }
#
#
#
#   kruskal_results <-    case_df %>%
#     group_by(Protein) %>%
#     kruskal_test(Intensity ~ Case) %>%
#     adjust_pvalue(method = "BH") %>%
#     separate(Protein, sep = '\\|', c("sp", "Protein.ID", "Accession", "Peptide"))
#
#   kruskal_results <-
#     cbind(kruskal_results,
#           start_seq = NA,
#           end_seq = NA)
#
#   if (is.function(updateProgress)) {
#     text <- "Kruskal Complete"
#     updateProgress(detail = text)
#   }
#
#   return(as.data.frame(kruskal_results))
# }
#
#
#
# getANOVA <- function(case_df, updateProgress = NULL) {
#   if (is.function(updateProgress)) {
#     text <- "Performing ANOVA (This may take a while)"
#     updateProgress(detail = text)
#   }
#
#
#   anova_results <- case_df %>%
#     drop_na(Intensity) %>%
#     group_by(Protein) %>%
#     rstatix::anova_test(Intensity ~ Case, detailed = T) %>%
#     adjust_pvalue(method = "BH") %>%
#     as_tibble() %>%
#     separate(Protein, sep = '\\|', c("sp", "Protein.ID", "Accession", "Peptide")) %>%
#     magrittr::set_class(c("anova_test", "rstatix_test", "data.frame"))
#
#
#
#   if (is.function(updateProgress)) {
#     text <- "ANOVA Complete"
#     updateProgress(detail = text)
#   }
#
#   return(anova_results)
#
#
# }



#The two following functions runs statistical tests using the data inputted by
#the user after it has been cleaned outliers removed from the dataframe. Returns
#dataframe with statistical results.

getTtest <- function(case_df, p, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    text <- "Performing Ttest (This may take a while)"
    updateProgress(detail = text)
  }

  ttest_results <-  case_df %>%


    #Run the T-test and adjustments
    group_by(Protein) %>%
    t_test(Intensity ~ Case, detailed = T) %>%
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
    cbind(ttest_results, start_seq = NA, end_seq = NA)
  ttest_results$File <- unique(case_df$File)

  if (ttest_results$group1[1] == "CTR") {
    ttest_results %<>% rename(group1 = "CTR", group2 = "Case")

  }

  else{
    ttest_results %<>% rename(group2 = "CTR", group1 = "Case")

  }

  if (is.function(updateProgress)) {
    text <- "Ttest Complete"
    updateProgress(detail = text)
  }

  return(as.data.frame(ttest_results))
}





getMannWhit <- function(case_df, p, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    text <- "Performing MannWhitney (This may take a while)"
    updateProgress(detail = text)
  }

  mannwhit_results <-  case_df %>%


    #Run the T-test and adjustments
    group_by(Protein) %>%
    wilcox_test(Intensity ~ Case, detailed = T) %>%
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
    cbind(mannwhit_results, start_seq = NA, end_seq = NA)
  mannwhit_results$File <- unique(case_df$File)


  if (mannwhit_results$group1[1] == "CTR") {
    mannwhit_results %<>% rename(group1 = "CTR", group2 = "Case")

  }

  else{
    mannwhit_results %<>% rename(group2 = "CTR", group1 = "Case")

  }

  if (is.function(updateProgress)) {
    text <- "MannWhit Complete"
    updateProgress(detail = text)
  }

  return(as.data.frame(mannwhit_results))
}




#Find the peptide location in the larger protein string and then marks the start
#and end point into two separate columns. These are then used for plotting in the
#plotting loop. Returns dataframe with peptide locations appeneded to main dataframe.

locatePeptides <- function(df, updateProgress = NULL) {
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

# cleanPCAVol <- function(df,conditionA, conditionB, labelA, labelB){
#
# }

runPCA <- function(df,condition, label){
  df <- df[!is.na(df$Protein.Existence),]
  all_labels <- c(label)


  selected_data <- df %>%
    dplyr::select(2:8, dplyr::contains("MaxLFQ", ignore.case = TRUE)) ## Update for other file input types

  long_data <- selected_data %>%
    pivot_longer(
      cols = contains("MaxLFQ", ignore.case = TRUE),
      names_to = "Sample",
      values_to = "Intensity"
    )
if(length(label) == 6){
  long_data <- long_data %>%
    mutate(
      Group = case_when(

        str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
        str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
        str_detect(Sample, paste(condition[[3]])) ~ label[[3]],
        str_detect(Sample, paste(condition[[4]])) ~ label[[4]],
        str_detect(Sample, paste(condition[[5]])) ~ label[[5]],
        str_detect(Sample, paste(condition[[6]])) ~ label[[6]],
        TRUE                                        ~ "Other"
      ))

}
 else if(length(label) == 5){
    long_data <- long_data %>%
      mutate(
        Group = case_when(

          str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
          str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
          str_detect(Sample, paste(condition[[3]])) ~ label[[3]],
          str_detect(Sample, paste(condition[[4]])) ~ label[[4]],
          str_detect(Sample, paste(condition[[5]])) ~ label[[5]],

          TRUE                                        ~ "Other"
        ))

 }
  else if(length(label) == 4){
    long_data <- long_data %>%
      mutate(
        Group = case_when(

          str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
          str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
          str_detect(Sample, paste(condition[[3]])) ~ label[[3]],
          str_detect(Sample, paste(condition[[4]])) ~ label[[4]],

          TRUE                                        ~ "Other"
        ))

  }
  else if(length(label) == 3){
    long_data <- long_data %>%
      mutate(
        Group = case_when(

          str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
          str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
          str_detect(Sample, paste(condition[[3]])) ~ label[[3]],

          TRUE                                        ~ "Other"
        ))

  }
else{
  long_data <- long_data %>%
    mutate(
      Group = case_when(

        str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
        str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
        TRUE                                        ~ "Other"
      ))
}

  long_data_log <- long_data %>%
    mutate(log2Intensity = log2(long_data$Intensity)) %>%
    filter(!is.infinite(log2Intensity))

  filtered_proteins <- long_data_log %>%
    group_by(`Protein.ID`, Group) %>%                    # group by protein & condition
    summarise(valid_n = sum(!is.na(log2Intensity)), .groups = "drop") %>%
    group_by(`Protein.ID`) %>%
    summarise(max_valid = max(valid_n), .groups = "drop") %>%
    filter(max_valid >= 2)


  filtered_matrix_4min <- long_data_log %>%
    filter(`Protein.ID` %in% filtered_proteins$`Protein.ID`)

  #(optional): Pivot back to wide format ---
  filtered_matrix_4min_wide <- filtered_matrix_4min %>%
    dplyr::select(`Protein.ID`, Gene, Sample, log2Intensity) %>%
    pivot_wider(names_from = Sample, values_from = log2Intensity)


  expr <- filtered_matrix_4min_wide %>% select(-`Protein.ID`, -Gene)
  sample_names <- colnames(expr)


  if(length(label) == 6){
    group_map <- tibble(Sample = sample_names) %>%
      mutate(
        Group = case_when(

          str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
          str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
          str_detect(Sample, paste(condition[[3]])) ~ label[[3]],
          str_detect(Sample, paste(condition[[4]])) ~ label[[4]],
          str_detect(Sample, paste(condition[[5]])) ~ label[[5]],
          str_detect(Sample, paste(condition[[6]])) ~ label[[6]],
          TRUE                                        ~ "Other"
        ))

  }
  else if(length(label) == 5){
    group_map <- tibble(Sample = sample_names) %>%
      mutate(
        Group = case_when(

          str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
          str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
          str_detect(Sample, paste(condition[[3]])) ~ label[[3]],
          str_detect(Sample, paste(condition[[4]])) ~ label[[4]],
          str_detect(Sample, paste(condition[[5]])) ~ label[[5]],

          TRUE                                        ~ "Other"
        ))

  }
  else if(length(label) == 4){
    group_map <- tibble(Sample = sample_names) %>%
      mutate(
        Group = case_when(

          str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
          str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
          str_detect(Sample, paste(condition[[3]])) ~ label[[3]],
          str_detect(Sample, paste(condition[[4]])) ~ label[[4]],

          TRUE                                        ~ "Other"
        ))

  }
  else if(length(label) == 3){
    group_map <- tibble(Sample = sample_names) %>%
      mutate(
        Group = case_when(

          str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
          str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
          str_detect(Sample, paste(condition[[3]])) ~ label[[3]],

          TRUE                                        ~ "Other"
        ))

  }
  else{
    group_map <- tibble(Sample = sample_names) %>%
      mutate(
        Group = case_when(

          str_detect(Sample, paste(condition[[1]])) ~ label[[1]],
          str_detect(Sample, paste(condition[[2]])) ~ label[[2]],
          TRUE                                        ~ "Other"
        ))
  }

  subset_samples <- group_map %>%
    filter(Group %in% all_labels) %>%
    pull(Sample)


  expr_sub <- expr[, subset_samples, drop = FALSE]

  # 4) Impute NAs by per-protein (row) median and transpose (samples x proteins)
  expr_mat <- as.matrix(expr_sub)
  mode(expr_mat) <- "numeric"

  row_meds <- apply(expr_mat, 1, function(x)
    median(x, na.rm = TRUE))
  keep_rows <- is.finite(row_meds)
  expr_mat  <- expr_mat[keep_rows, , drop = FALSE]
  row_meds  <- row_meds[keep_rows]

  na_idx <- which(is.na(expr_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0)
    expr_mat[na_idx] <- row_meds[na_idx[, 1]]

  expr_t <- t(expr_mat)
  rownames(expr_t) <- colnames(expr_mat)

  # 5) Remove zero-variance proteins (needed for scale.=TRUE)
  sds <- apply(expr_t, 2, sd, na.rm = TRUE)
  keep_cols <- is.finite(sds) & sds > 1e-8
  expr_t2 <- expr_t[, keep_cols, drop = FALSE]

  # 6) PCA
  pca <- prcomp(expr_t2, center = TRUE, scale. = TRUE)
  var_expl <- (pca$sdev ^ 2 / sum(pca$sdev ^ 2)) * 100
  pc1_lab <- paste0("PC1 (", round(var_expl[1], 1), "%)")
  pc2_lab <- paste0("PC2 (", round(var_expl[2], 1), "%)")

  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$Sample <- rownames(scores)
browser()

  plot_df <- scores %>%
    left_join(group_map, by = "Sample") %>%
    filter(Group %in% all_labels)

  # 8) Publication-style PCA plot
  p_pca_all <- ggplot(plot_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 3,
               alpha = 0.95,
               shape = 16) +
    stat_ellipse(
      aes(group = Group, color = Group),
      level = 0.95,
      linewidth = 0.8,
      linetype = 2
    ) +

    labs(x = pc1_lab, y = pc2_lab, title = "PCA") +
    theme_bw(base_size = 14) +
    coord_cartesian() +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      axis.text  = element_text(color = "black", size = 8),
      axis.title = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 8),
      plot.margin = margin(10, 10, 10, 10)
    )

  p_pca_all

}


runVol <- function(df,conditionA, conditionB, labelA, labelB, genes){



  df <- df[!is.na(df$Protein.Existence),]
  all_labels <- c(labelA,labelB)


  selected_data <- df %>%
    dplyr::select(2:8, dplyr::contains("MaxLFQ", ignore.case = TRUE)) ## Update for other file input types

  long_data <- selected_data %>%
    pivot_longer(
      cols = contains("MaxLFQ", ignore.case = TRUE),
      names_to = "Sample",
      values_to = "Intensity"
    )

  long_data <- long_data %>%
    mutate(
      Group = case_when(

        str_detect(Sample, paste(conditionA)) ~ labelA,
        str_detect(Sample, paste(conditionB)) ~ labelB,
        TRUE                                        ~ "Other"
      ))


  long_data_log <- long_data %>%
    mutate(log2Intensity = log2(long_data$Intensity)) %>%
    filter(!is.infinite(log2Intensity))

  filtered_proteins <- long_data_log %>%
    group_by(`Protein.ID`, Group) %>%                    # group by protein & condition
    summarise(valid_n = sum(!is.na(log2Intensity)), .groups = "drop") %>%
    group_by(`Protein.ID`) %>%
    summarise(max_valid = max(valid_n), .groups = "drop") %>%
    filter(max_valid >= 2)


  filtered_matrix_4min <- long_data_log %>%
    filter(`Protein.ID` %in% filtered_proteins$`Protein.ID`)

  #(optional): Pivot back to wide format ---
  filtered_matrix_4min_wide <- filtered_matrix_4min %>%
    dplyr::select(`Protein.ID`, Gene, Sample, log2Intensity) %>%
    pivot_wider(names_from = Sample, values_from = log2Intensity)


  expr <- filtered_matrix_4min_wide %>% select(-`Protein.ID`, -Gene)
  sample_names <- colnames(expr)



  group_map <- tibble(Sample = sample_names) %>%
    mutate(
      Group = case_when(

        str_detect(Sample, paste(conditionA)) ~ labelA,
        str_detect(Sample, paste(conditionB)) ~ labelB,
        TRUE                                        ~ "Other"
      ))

  subset_samples <- group_map %>%
    filter(Group %in% all_labels) %>%
    pull(Sample)


  expr_sub <- expr[, subset_samples, drop = FALSE]

  # 4) Impute NAs by per-protein (row) median and transpose (samples x proteins)
  expr_mat <- as.matrix(expr_sub)
  mode(expr_mat) <- "numeric"

  row_meds <- apply(expr_mat, 1, function(x)
    median(x, na.rm = TRUE))
  keep_rows <- is.finite(row_meds)
  expr_mat  <- expr_mat[keep_rows, , drop = FALSE]
  row_meds  <- row_meds[keep_rows]

  na_idx <- which(is.na(expr_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0)
    expr_mat[na_idx] <- row_meds[na_idx[, 1]]

  expr_t <- t(expr_mat)
  rownames(expr_t) <- colnames(expr_mat)

  # 5) Remove zero-variance proteins (needed for scale.=TRUE)
  sds <- apply(expr_t, 2, sd, na.rm = TRUE)
  keep_cols <- is.finite(sds) & sds > 1e-8
  expr_t2 <- expr_t[, keep_cols, drop = FALSE]


  # 1) Select A & B columns from your 4-min wide matrix (log2)
  all_samples <- dplyr::setdiff(colnames(filtered_matrix_4min_wide), c("Protein.ID", "Gene"))
  A_cols  <- grep(labelA, all_samples, ignore.case = TRUE, value = TRUE) #Label
  B_cols <- grep(labelB, all_samples, ignore.case = TRUE, value = TRUE) #Label
  stopifnot(length(A_cols) > 0, length(B_cols) > 0)

  sel_cols <- c(A_cols, B_cols)

  wide_bd <- filtered_matrix_4min_wide %>%
    dplyr::select(`Protein.ID`, Gene, dplyr::all_of(sel_cols))

  #  2) Impute NAs by per-protein (row) median across selected samples (still log2)
  expr_mat <- as.matrix(wide_bd %>% dplyr::select(-`Protein.ID`, -Gene))
  mode(expr_mat) <- "numeric"

  row_medians <- apply(expr_mat, 1, function(x)
    median(x, na.rm = TRUE))
  keep_rows   <- is.finite(row_medians)
  expr_mat    <- expr_mat[keep_rows, , drop = FALSE]
  row_medians <- row_medians[keep_rows]

  na_idx <- which(is.na(expr_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0)
    expr_mat[na_idx] <- row_medians[na_idx[, 1]]

  expr_imp <- as.data.frame(expr_mat) %>%
    dplyr::mutate(`Protein.ID` = wide_bd$`Protein.ID`[keep_rows],
                  Gene         = wide_bd$Gene[keep_rows])

  # 3) Long format + group labels (A vs B only)
  long_imp <- expr_imp %>%
    tidyr::pivot_longer(
      cols = -c(`Protein.ID`, Gene),
      names_to = "Sample",
      values_to = "log2Intensity"
    ) %>%
    dplyr::mutate(Group = dplyr::case_when(
      str_detect(Sample, conditionA) ~ "A", #Label A
      str_detect(Sample, conditionB) ~ "B", #Label B
      TRUE ~ "Other"
    )) %>%
    dplyr::filter(Group %in% c("A", "B"))

  # 4) Per-protein stats: log2FC + pval (Welch; Wilcoxon fallback on zero variance)
  stats_long <- long_imp %>%
    dplyr::group_by(`Protein.ID`, Gene, Group) %>%
    dplyr::summarise(
      vals = list(log2Intensity),
      mean = mean(log2Intensity),
      n    = n(),
      .groups = "drop"
    )
  stats_4 <- stats_long %>%
    tidyr::pivot_wider(
      id_cols = c(`Protein.ID`, Gene),
      names_from = Group,
      values_from = c(vals, mean, n),
      names_sep = ".") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pval = {
      v_A  <- tryCatch(
        `vals.A`,
        error = function(e)
          numeric(0)
      )
      v_B <- tryCatch(
        `vals.B`,
        error = function(e)
          numeric(0)
      )
      if (length(v_A) >= 2 && length(v_B) >= 2) {
        sd_A <- stats::sd(v_A)
        sd_B <- stats::sd(v_B)
        if (is.finite(sd_A) &&
            is.finite(sd_B) && sd_A > 0 && sd_B > 0) {
          stats::t.test(v_A, v_B, var.equal = FALSE)$p.value
        } else if (length(unique(c(v_A, v_B))) > 1) {
          suppressWarnings(stats::wilcox.test(v_A, v_B, exact = FALSE)$p.value)
        } else
          1
      } else
        NA_real_
    }, # A - B (positive = up in A)
    log2FC = {
      m_A  <- tryCatch(
        `mean.A`,
        error = function(e)
          NA_real_
      )
      m_B <- tryCatch(
        `mean.B`,
        error = function(e)
          NA_real_
      )
      m_A - m_B
    }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(`Protein.ID`, Gene, log2FC, pval)

  stats_4$Gene <- toupper(stats_4$Gene)
  stats_4 <- stats_4 %>% adjust_pvalue("pval", "adjusted","BH") #add input for correction type
  top_genes <- as.data.frame(head(stats_4[order(stats_4$adjusted), "Gene"], 50)) #input for n
  # 5) Build plotting df with consistent calls & labels
  df <- stats_4 %>%
    dplyr::mutate(
      diffexpressed = dplyr::case_when(
        is.finite(log2FC) &
          is.finite(adjusted) & adjusted < 0.05 & log2FC >  0 ~ "Up",
        is.finite(log2FC) &
          is.finite(adjusted) & adjusted < 0.05 & log2FC <  0 ~ "Down",
        TRUE ~ "NS"
      ),
      # Edit the set below to choose which genes to label:
      delabel = ifelse(stats_4$Gene %in% top_genes$Gene, stats_4$Gene, NA)
    )


  #  6) Nature-style volcano with boxed labels + internal legend
  # y-limits with headroom
  y_max <- max(-log10(df$adjusted), na.rm = TRUE)
  y_lim <- c(0, y_max + 0.6)

  # factor & palette
  df$diffexpressed <- factor(df$diffexpressed, levels = c("Up", "Down", "NS"))
  volc_cols <- c(
    "Up" = "#B2182B",
    "Down" = "#2166AC",
    "NS" = "#7A7A7A"
  )

  nature_theme <- theme_classic(base_size = 9) +
    theme(
      panel.border   = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.6
      ),
      panel.grid     = element_blank(),
      axis.line      = element_blank(),
      axis.ticks     = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(2.2, "pt"),
      axis.text      = element_text(color = "black", size = 14),
      axis.title     = element_text(
        color = "black",
        size = 14,
        margin = margin(4, 4, 4, 4)
      ),
      plot.title     = element_text(hjust = 0.5, size = 14),
      plot.margin    = margin(6, 6, 6, 6)
    )

  volcano_plot_nature <- ggplot(df, aes(
    x = log2FC,
    y = -log10(adjusted),
    color = diffexpressed
  )) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey40"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey40"
    ) +
    geom_point(size = 2, alpha = 0.5) +
    scale_color_manual(values = volc_cols, drop = FALSE) +
    coord_cartesian(
      ylim = y_lim,
      xlim = c(-4.2, 2),
      expand = FALSE
    ) +
    scale_x_continuous(breaks = seq(-10, 10, 2)) +
    labs(
      title = paste(labelA, "vs.", labelB),
      x = expression("log"[2] * "FC"),
      y = expression("-log"[10] * "p")
    ) +
    nature_theme

  label_df <- subset(df, delabel != "" &
                       is.finite(adjusted) & is.finite(log2FC))

  volcano_plot_nature_box_inside <- volcano_plot_nature +
    ggrepel::geom_label_repel(
      data = label_df,
      aes(label = delabel, color = diffexpressed),
      fill = alpha("white", 0.9),
      label.size = 0,
      size = 4,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.25,
      segment.color = "grey60",
      segment.size = 0.25,
      segment.curvature = 0.1,
      segment.ncp = 2,
      segment.angle = 20,
      max.overlaps = Inf,
      na.rm = TRUE,
      show.legend = FALSE   # << stops labels from changing the legend glyph
    ) +
    guides(
      color = guide_legend(
        title = "Expression",
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(
          shape = 16,
          size = 3,
          alpha = 0.9
        )  # << keeps colored circles
      )
    ) +
    theme(
      legend.position = c(0.95, 0.93),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80"),
      legend.key.size = unit(8, "pt"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8, face = "bold")
    )

  volcano_plot_nature_box_inside



}


runCluster <- function(df){

}
