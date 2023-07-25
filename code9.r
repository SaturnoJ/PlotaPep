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


##Read the deasease/control files and read the fasta files
dataframe <- read_delim("C:\\Users\\Jan\\Downloads\\combined_peptide.tsv")
fasta <- as.data.frame(readAAStringSet("C:\\Users\\Jan\\combined_fasta.fasta"))
fasta <- tibble::rownames_to_column(fasta,"accession")
colnames(fasta) <- c("accession","seq")


##Extract the top n most frequent proteins
extract_top <- function(dataset, n) {
    top_n(dataset, n)
}

##Make a dataframe with the top n proteins and extract the top 5 most frequent ones.
protein <- dataframe$`Protein ID`
protein_table <- table(protein)
df_protein <- protein_table %>% 
    as.data.frame() %>% 
    arrange(desc(Freq))
top_n_df <- extract_top(df_protein, 5)

##Create a function that cleans the data and produces and stores a plot of the data
Visualize_data <- function(data, protein, fasta){
    # clean the data
    clean_combined_peptide_fragger <- function(dataset, Intensity, MaxLFQ, cutoff){
        
        
        #change Protein name
        names(dataset)[names(dataset) == 'Protein ID'] <- 'ProteinID'
        names(dataset)[names(dataset) == 'Peptide Sequence'] <- 'Sequence'
        
        #Select columns of interest
        dataset <- dplyr::mutate(dataset, Protein = paste0(ProteinID, "_", Gene, "_", Sequence))
        
        # Select the specified Intensity. 
        # We also start creating the substring we will remove from sample names.
        if (Intensity == "Intensity" || Intensity == "intensity"){
            dataset <- dplyr::select(dataset, Protein, contains("Intensity"), -contains("Total"), -contains("Unique"))
            removeSubString <- paste0(" ", "Intensity")
        } else if (Intensity == "Total" || Intensity == "total"){
            dataset <- dplyr::select(dataset, Protein, contains("Total Intensity"))
            removeSubString <- paste0(" ", "Total", " Intensity")
        } else if (Intensity == "Unique" || Intensity == "unique"){
            dataset <- dplyr::select(dataset, Protein, contains("Unique Intensity"))
            removeSubString <- paste0(" ", "Unique", " Intensity")
        } else {
            print("Something went wrong with selecting which Intensity. You sure you wrote it right? Choose from Total, Unique or Intensity")
        }
        
        #select or not select MaxLFQ
        if (MaxLFQ == "MaxLFQ" || MaxLFQ == "maxlfq" || MaxLFQ == "MAXLFQ" || MaxLFQ == "max lfq"){
            dataset <- dplyr::select(dataset, Protein, contains("MaxLFQ"))
            removeSubString <- paste0(" MaxLFQ", removeSubString)
        } else if(MaxLFQ == "None" || MaxLFQ == "none"){
            dataset <- dplyr::select(dataset, Protein, contains("Intensity"), -contains("MaxLFQ"))
        } else {
            print("Something went wrong. You have to write MaxLFQ or None(for no MAXLFQ selection)")
            print("Alternativly, you might have selected MaxLFQ but you used old version of Fragpipe (its been there since 17.0/17.1")
        }
        
        #col to rowname
        dataset <- tibble::column_to_rownames(dataset, "Protein")
        
        #Get rid of substring part we do not need
        colnames(dataset) <- stringr::str_remove(colnames(dataset), removeSubString)
        
        #transpose, remove zeros, log2
        dataset <- data.frame(t(dataset))
        dataset[dataset == 0] <- NA
        dataset <- log2(dataset)
        
        #apply cutoff on the protein level
        if (cutoff %in% 0:100){
            dataset <- purrr::discard(dataset, ~sum(is.na(.x))/length(.x)* 100 >= (100-cutoff))
        } else {
            print("You selected a cut off filter out of 0 until 100.....")
            print("Or maybe you wrote it as a string (if it has quotation marks around it)")
        }
        
        #Long Format it
        dataset <- tibble::rownames_to_column(dataset, "Sample")
        dataset <- tidyr::gather(dataset, colnames(dataset)[2:ncol(dataset)], key = "Protein", value = "Intensity") 
        
        return(dataset)
    }
    ##Run the code to clean the data
    outlier_removed_df <- clean_combined_peptide_fragger(data, "intensity", "MaxLFQ", 95)
    
    ##Remove the outliers
    OutlierRemover <- function(dataset){
        
        #Determine correlations between sample and theoretical
        SampleCorrelations <- dataset%>%
            group_by(Protein) %>%
            dplyr::mutate(medianProtein = median(Intensity, na.rm = T)) %>%
            ungroup() %>%
            group_by(Sample) %>%
            dplyr::mutate(SampleCor = stats::cor(Intensity, medianProtein, method = "pearson", use="pairwise.complete.obs")) %>%
            ungroup() %>%
            distinct(Sample, .keep_all = T)
        
        SampleCorrelations <- SampleCorrelations %>%
            dplyr::mutate(Outlier = SampleCor < (1 - 1*sd(unique(SampleCorrelations$SampleCor), na.rm = T)))
        
        #Make the DF that we will return (i.e. without outliers)
        returnDataset <- dataset %>%
            dplyr::filter(Sample %in% dplyr::filter(SampleCorrelations, Outlier == F)$Sample) %>%
            dplyr::select(Protein, Sample, Intensity)  
        
        #Half minimum value (per protein) imputation
        PCA_DF <- dataset %>% 
            group_by(Protein) %>%
            dplyr::mutate(IntensImputed = replace_na(Intensity, mean(Intensity, na.rm = T)/2)) %>%
            dplyr::select(Protein, Sample, IntensImputed) %>%
            spread(key = "Protein", value = "IntensImputed") %>%
            column_to_rownames("Sample")
        
        #scale and PCA
        PCA_DF_Results <- prcomp(scale(as.matrix(PCA_DF)))
        eigs <- PCA_DF_Results$sdev^2
        variance_percentage <- (eigs / sum(eigs))*100
        pc1var <- round(variance_percentage[1],digits=0)
        pc2var <- round(variance_percentage[2],digits=0)
        
        return(returnDataset)
    }
    
    ##Run the code that removes outliers
    outlier_removed_df <- OutlierRemover(outlier_removed_df)
    
    ##Apply protein level cutoff
    combined_cutoff_df <- outlier_removed_df %>% 
        mutate(uniqueSamplesInDF = length(unique(Sample))) %>%
        group_by(Protein) %>%
        mutate(countMissingValues = sum(!is.na(Intensity))) %>%
        mutate(percentageMissing = countMissingValues / uniqueSamplesInDF) %>%
        filter(percentageMissing > 0.95) %>% # Apply cutoff of 95%
        dplyr::select(Sample, Protein, Intensity) #clean up by selecting the original columns
    combined_cutoff_df <- as.data.frame(combined_cutoff_df)
    
    ##extracts the alzheimers dataset from the merged dataset
    split_ad <- function(data, ctr) {
        df <- subset(combined_cutoff_df, !grepl("CTR", combined_cutoff_df$Sample))
        if (isEmpty(df)) {
            return(NULL)
        }
        return(df)
        
    }
    
    ##extracts the control dataset from the merged dataset
    split_ctr <- function(combined_cutoff_df) {
        df <- subset(combined_cutoff_df, grepl("CTR", combined_cutoff_df$Sample))
        if (isEmpty(df)) {
            return(NULL)
        }
        return(df)
    }
    
    ##Add a row to the control and alzheimers dataset that shows if it is alzheimers or control and also takes out any points and underscores that interfere with splitting the protein gene and peptide sequence
    control_df <- split_ctr(combined_cutoff_df)
    control_df <- cbind(control_df, "CTR")
    colnames(control_df) <- c('Sample','Protein','Intensity', 'DX')
    alzheimers_df <- split_ad(combined_cutoff_df, control_df)
    alzheimers_df <- cbind(alzheimers_df, "Alz")
    colnames(alzheimers_df) <- c('Sample','Protein','Intensity', 'DX')
    combined_df <- rbind(control_df, alzheimers_df)
    combined_df$Protein <- gsub('\\.','A',combined_df$Protein)
    for (i in 1:nrow(combined_df)){
        fn <- "_"
        rp <- "A"
        n <- 2
        a <- str_count(combined_df$Protein[i],"_")
        if (a == 3){
            regmatches(combined_df$Protein[i], gregexpr(fn, combined_df$Protein[i])) <- list(c(rep(fn,n-1),rp))
        }
        i <- i -1
    }
    
    ##Apply a peptide level cutoff of 95%
    ttest_results <- combined_df %>%
        
        
        #Run the T-test and adjustments
        group_by(Protein) %>%
        t_test(Intensity ~ DX, detailed = T) %>%
        adjust_pvalue(method = "BH") %>%
        
        #Split the Protein name in Uniprot and Gene
        separate(Protein, c("UniprotID", "Gene","Peptide")) %>%
        
        #Determine Fold change. Since we work with log-transformed values we can just substract
        mutate(FC = estimate1 - estimate2) %>%
        
        #Create log10 p-vals
        mutate(log10adjustP = -1*log10(p.adj)) %>%
        
        #Determine if up or down regulated
        mutate(Direction = ifelse(p.adj > 0.05, "NotSignificant", ifelse(FC < 0, "Down", "Up")))
    
    ##Locate the peptides in the protein sequence and add the start and end point to the dataframe 
    located_peptides <- ttest_results %>%
        filter(UniprotID == protein)
    located_peptides <- cbind(located_peptides, start_seq=NA, end_seq=NA)
    locate_peptides <- function(cleavages, fasta) {
        for (i in 1:nrow(cleavages)) {
            locate <- as.data.frame(str_locate(toString(fasta$seq), toString(cleavages$Peptide[i])))
            cleavages$start_seq[i] <- locate$start
            cleavages$end_seq[i] <- locate$end
        }
        return(cleavages)
    }
    located_peptides <- locate_peptides(located_peptides, fasta)
    
    
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
        geom_vline(xintercept = 1, lwd = 2, alpha = .5) +
        geom_vline(xintercept = protein_length, lwd = 2, alpha = 0.5) +
        
        #set a horizontal line to inform about FC = 0
        geom_hline(yintercept = 0, color = "black") +
        
        #Here we build the peptide "blocks". Do note that all column-info goes INSIDE the aes() part.
        geom_rect(
            aes(
                xmin = located_peptides$start_seq,
                xmax = (located_peptides$start_seq + (located_peptides$end_seq - located_peptides$start_seq)),
                ymin = filtered_results$FC-0.05,
                ymax = filtered_results$FC+0.05,
                fill = isSignificant
            ),
            
            #Here I specify some stuff that will be universal for all blocks irrespective of column info.
            col = "black",
            alpha = 0.75) +
        
        #Set the Ggplot theme, limit the y-axis for FC.
        theme_bw() +
        ylim(-2, 2) +
        
        #Specify the colours I want to use for the isSignificant column
        scale_fill_manual(values = c("yes" = "red","no" = "grey")) +
        theme(legend.position = "bottom") +
        
        #x and yaxis titles
        xlab("Protein Sequence") +
        ylab("FC") 
}
##Make an array to store the data and run the function for all the top n proteins
array_plots <- c()
top_n_df <- top_n_df[-c(1,3),]
for (i in 1:nrow(top_n_df))
    Visualize_data(dataframe, top_n_df$protein[i], fasta[i])
array_plots <- append(y,Visualize_data)
