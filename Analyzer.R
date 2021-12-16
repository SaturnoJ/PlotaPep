#Function removes trypic cleavages from full dataset by using the subset function and grepl for boolean values. Subset and grepl are funcitons of base R
remove_tryp <- function(symboled) {
  return(subset(
    symboled,
    (
      !grepl('\\R$', symboled$Sequence) &
        !grepl('\\K$', symboled$Sequence) #&
      #symboled$`Peptide Next AA` != "P"
    ) |
      (
        symboled$`Peptide Previous AA`  != "R" &
          symboled$`Peptide Previous AA` != "K"
      )
  ))
  
}
#Function finds semi-trypic cleavages in dataset by looping through comparing the amino acid before the sequence and ending amino acid in the sequence.
find_semi <- function(filtered) {
  for (i in 1:nrow(filtered)) {
    if (filtered$`Peptide Previous AA`[i]  == "R" |
        filtered$`Peptide Previous AA`[i] == "K") {
      filtered$N_terminus[i] <- "True"
    } else {
      filtered$N_terminus[i] <- "False"
    }
    
    if ((grepl('\\R$', filtered$Sequence[i]) ||
         grepl('\\K$', filtered$Sequence[i])))
    {
      filtered$C_terminus[i] <- "True"
    } else{
      filtered$C_terminus[i] <- "False"
      
    }
    
    
  }
  return(filtered)
}
#Function removes the control cohort from the main dataset.
split_ctr <- function(semi_tryp) {
  return(subset(semi_tryp, grepl("CTR", semi_tryp$file)))
  
}
#Uses the control file to intersect with main dataset to get Alzheimer's Dementia cohort from set.
split_ad <- function(semi_tryp, ctr) {
  return(subset(semi_tryp,!(semi_tryp$file %in% ctr$file)))
  
}
#Loads already found semi-trypic values into aggregated dataset.
load_semi <- function(cleavage_x, semi_tryp) {
  for (i in 1:nrow(cleavage_x)) {
    for (j in 1:nrow(semi_tryp)) {
      if (cleavage_x$Var1[i] == semi_tryp$Sequence[j]) {
        cleavage_x$N_terminus[i] <- semi_tryp$N_terminus[j]
        cleavage_x$C_terminus[i] <- semi_tryp$C_terminus[j]
        break
      }
      
    }
    
  }
  return(cleavage_x)
}
#Loads protein symbols into aggregate dataset
load_protein <- function(cleavage_x, semi_tryp) {
  for (i in 1:nrow(semi_tryp)) {
    for (j in 1:nrow(cleavage_x)) {
      if (semi_tryp$Sequence[i] == cleavage_x$Var1[j] &&
          !is.null(cleavage_x$protein[j])) {
        cleavage_x$protein[j] <- semi_tryp$Accession[i]
        
      }
      
    }
    
  }
  return(cleavage_x)
  
}
#Loads where missed cleavage location occurred into aggregate dataset
cut_at <- function(cleavage_x) {
  for (i in 1:nrow(cleavage_x)) {
    if (cleavage_x$N_terminus[i] == 'False') {
      cleavage_x$cleavage_loc[i] <- cleavage_x$start_seq[i]
      
    }
    
    else if (cleavage_x$C_terminus[i] == 'False')
      
      cleavage_x$cleavage_loc[i] <- cleavage_x$end_seq[i]
    
  }
  
  return(cleavage_x)
}
#Filters out proteins above cutoff and not being searched for
filterer <- function(df, cutoff, symbols) {
  filtered <- filter(df, cutoff > df$Expect)
  x<-NULL
  for(i in 1:nrow(symbols)){
  y <- filter(filtered,
                     grepl(symbols$UNIPROT[i],filtered$Accession))
  x <- rbind(x,y)
  }
  filtered <-x
  
  return(filtered)
}
#Adds symbols to Accession column in dataframe
add_symbols <- function(filtered, symbols) {
  for (i in 1:nrow(symbols)) {
    for (j in 1:nrow(filtered)) {
      if (grepl(symbols$UNIPROT[i], filtered$Accession[j])) {
        filtered$Accession[j] <- symbols$SYMBOL[i]
        
      }
    }
  }
  return(filtered)
}

locate <- function(cleavage_x, fasta) {
  
  
  for (i in 1:nrow(fasta)) {
    for(j in 1:nrow(cleavage_x)){
    if (grepl(cleavage_x$protein[j], fasta$accession[i])) {
      locate <-
        as.data.frame(str_locate(fasta$seq[i], cleavage_x$Var1[j]))
      cleavage_x$start_seq[j] <- locate$start
      cleavage_x$end_seq[j] <- locate$end
    }
  }
}
  return(cleavage_x)
  
}

setwd("V:/Jason/HendriksFilesSemiTrytpic/CSV")
path <- "V:/Jason/HendriksFilesSemiTrytpic/CSV"
file.names <- list.files(path, pattern = ".csv")

library("tidyverse")
library("org.Hs.eg.db")
library("biomartr")
library("magrittr")
library("Biostrings")

ptm <- proc.time()

# file <-
#   read_delim(
#     file.names[1],
#     delim = "\t",
#     escape_double = FALSE,
#     trim_ws = TRUE,
#     show_col_types = FALSE
#   )
# file$file <- file.names[1]
# df <- file
# 
# for (i in 2:length(file.names)) {
#   file <-
#     read_delim(
#       file.names[i],
#       delim = "\t",
#       escape_double = FALSE,
#       trim_ws = TRUE,
#       show_col_types = FALSE
#     )
#   file$file <- file.names[i]
#   df <- rbind(df, file)
# }


fasta <- as.data.frame(readAAStringSet("protiens.fasta"))

fasta <- tibble::rownames_to_column(fasta, "accession")
colnames(fasta) <- c("accession", "seq")

uniprot <- c("P05067", "P02649", "P10636", "P10909", "P04156")


symbols <- select(org.Hs.eg.db, uniprot, "SYMBOL", "UNIPROT")



filtered <- filterer(df, 0.1, symbols)

symboled <- add_symbols(filtered, symbols)


symboled$N_terminus <- 'NULL'
symboled$C_terminus <- 'NULL'

removed <- remove_tryp(symboled)
semi_tryp <- find_semi(removed)
ctr <- split_ctr(semi_tryp)
ad <- split_ad(semi_tryp, ctr)


ad_seq <- as.data.frame(table(ad$Sequence))
ad_seq$Cohort <- c("AD")
ctr_seq <- as.data.frame(table(ctr$Sequence))
ctr_seq$Cohort <- c("CTR")

all_seq <- (rbind(ad_seq, ctr_seq))
all_seq$Var1 <- as.character(all_seq$Var1)

start_seq <- vector("numeric", nrow(all_seq))
N_terminus <- vector("character", nrow(all_seq))
end_seq <- vector("numeric", nrow(all_seq))
C_terminus <- vector("character", nrow(all_seq))
protein <- vector("character", nrow(all_seq))
cleavage_loc <- vector("numeric",nrow(all_seq))


cleavage_x <-
  data.frame(protein,all_seq,cleavage_loc, start_seq, end_seq, N_terminus, C_terminus)

cleavage_x <- load_semi(cleavage_x, semi_tryp)

cleavage_x <- load_protein(cleavage_x, semi_tryp)

cleavage_x <- locate(cleavage_x, fasta)

cleavage_x <- cut_at(cleavage_x)

# mapt_com <- ggplot(cleavage_x[cleavage_x$protein == 'MAPT', ],
#                    aes(
#                      x = cleavage_loc,
#                      y = Freq,
#                      group_by = cleavage_loc,
#                      fill = Cohort
#                    )) + geom_col() + theme_minimal() + labs(x = "Cleavage Site", y = "Frequency", title = "MAPT Cleavage Tabulation: AD vs CTRL") +  geom_text(
#                      data = subset(cleavage_x, Freq > 35 &
#                                      cleavage_x$protein == "MAPT"),
#                      aes(label = cleavage_loc),
#                      vjust = -1.0,
#                      hjust = -0.2,
#                      check_overlap =  TRUE
#                    ) + facet_wrap(vars(Cohort), nrow = 2,)
# app_com <- ggplot(cleavage_x[cleavage_x$protein == 'APP', ],
#                   aes(
#                     x = cleavage_loc,
#                     y = Freq,
#                     group_by = cleavage_loc,
#                     fill = Cohort
#                   )) + geom_col(width = 0.10) + theme_minimal() + labs(x = "Cleavage Site", y = "Frequency", title = "APP Cleavage Tabulation: AD vs CTRL") + geom_text(
#                     data = subset(cleavage_x, Freq > 5 &
#                                     cleavage_x$protein == "APP"),
#                     aes(label = cleavage_loc),
#                     vjust = -1.0,
#                     hjust = -0.2,
#                     check_overlap =  TRUE
#                   )  + facet_wrap(vars(Cohort), nrow = 2,)
# 
# apoe_com <- ggplot(cleavage_x[cleavage_x$protein == 'APOE', ],
#                    aes(
#                      x = cleavage_loc,
#                      y = Freq,
#                      group_by = cleavage_loc,
#                      fill = Cohort
#                    )) + geom_col() + theme_minimal() + labs(x = "Cleavage Site", y = "Frequency", title = "APOE Cleavage Tabulation: AD vs CTRL") + geom_text(
#                      data = subset(cleavage_x, Freq > 5 &
#                                      cleavage_x$protein == "APOE"),
#                      aes(label = cleavage_loc),
#                      vjust = -1.0,
#                      hjust = -0.2,
#                      check_overlap =  TRUE
#                    ) + facet_wrap(vars(Cohort), nrow = 2,)
# prnp_com <- ggplot(cleavage_x[cleavage_x$protein == 'PRNP', ],
#                    aes(
#                      x = cleavage_loc,
#                      y = Freq,
#                      group_by = cleavage_loc,
#                      fill = Cohort
#                    )) + geom_col(width = 0.25) + theme_minimal() + labs(x = "Cleavage Site", y = "Frequency", title = "PRNP Cleavage Tabulation: AD vs CTRL") + geom_text(
#                      data = subset(cleavage_x, Freq > 5 &
#                                      cleavage_x$protein == "PRNP"),
#                      aes(label = cleavage_loc),
#                      vjust = -1.0,
#                      hjust = -0.2,
#                      check_overlap =  TRUE
#                    ) + facet_wrap(vars(Cohort), nrow = 2,)
# clu_com <- ggplot(cleavage_x[cleavage_x$protein == 'CLU', ],
#                   aes(
#                     x = cleavage_loc,
#                     y = Freq,
#                     group_by = cleavage_loc,
#                     fill = Cohort
#                   )) + geom_col() + theme_minimal() + labs(x = "Cleavage Site", y = "Frequency", title = "CLU Cleavage Tabulation: AD vs CTRL") + geom_text(
#                     data = subset(cleavage_x, Freq > 5 &
#                                     cleavage_x$protein == "CLU"),
#                     aes(label = cleavage_loc),
#                     vjust = -1.0,
#                     hjust = -0.2,
#                     check_overlap =  TRUE
#                   ) + facet_wrap(vars(Cohort), nrow = 2,)
# 
# 
# ggsave(
#   file = "prnp_com.svg",
#   plot = prnp_com,
#   width = 10,
#   height = 8
# )
# ggsave(
#   file = "clu_com.svg",
#   plot = clu_com,
#   width = 10,
#   height = 8
# )
# 
# ggsave(
#   file = "apoe_com.svg",
#   plot = apoe_com,
#   width = 10,
#   height = 8
# )
# ggsave(
#   file = "app_com.svg",
#   plot = app_com,
#   width = 10,
#   height = 8
# )
# ggsave(
#   file = "mapt_com.svg",
#   plot = mapt_com,
#   width = 10,
#   height = 8
# )

proc.time() - ptm
end

