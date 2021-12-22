#Function removes trypic cleavages from full dataset by using the subset function and grepl for boolean values. Subset and grepl are funcitons of base R
remove_tryp <- function(df) {
  return(subset(
    df,
    (
      !grepl('\\R$', df$Sequence) &
        !grepl('\\K$', df$Sequence) #&
      #df$`Peptide Next AA` != "P"
    ) |
      (
        df$`Peptide Previous AA`  != "R" &
          df$`Peptide Previous AA` != "K"
      )
  ))
  
}
#Function finds semi-trypic cleavages in dataset by looping through comparing the amino acid before the sequence and ending amino acid in the sequence.
find_semi <- function(df) {
  for (i in 1:nrow(df)) {
    if (df$`Peptide Previous AA`[i]  == "R" |
        df$`Peptide Previous AA`[i] == "K") {
      df$N_terminus[i] <- "True"
    } else {
      df$N_terminus[i] <- "False"
    }
    
    if ((grepl('\\R$', df$Sequence[i]) ||
         grepl('\\K$', df$Sequence[i])))
    {
      df$C_terminus[i] <- "True"
    } else{
      df$C_terminus[i] <- "False"
      
    }
    
    
  }
  return(df)
}
#Function removes the control cohort from the main dataset.
split_ctr <- function(df) {
  return(subset(df, grepl("CTR", df$file)))
  
}
#Uses the control file to intersect with main dataset to get Alzheimer's Dementia cohort from set.
split_ad <- function(df, ctr) {
  return(subset(df,!(df$file %in% ctr$file)))
  
}
#Loads already found semi-trypic values into aggregated dataset.
load_semi <- function(cleavages, df) {
  for (i in 1:nrow(cleavages)) {
    for (j in 1:nrow(df)) {
      if (cleavages$Var1[i] == df$Sequence[j]) {
        cleavages$N_terminus[i] <- df$N_terminus[j]
        cleavages$C_terminus[i] <- df$C_terminus[j]
        break
      }
      
    }
    
  }
  return(cleavages)
}
#Loads protein symbols into aggregate dataset
load_protein <- function(cleavages, df) {
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(cleavages)) {
      if (df$Sequence[i] == cleavages$Var1[j] &&
          !is.null(cleavages$protein[j])) {
        cleavages$protein[j] <- df$Accession[i]
        
      }
      
    }
    
  }
  return(cleavages)
  
}
#Loads where missed cleavage location occurred into aggregate dataset
cut_at <- function(cleavages) {
  for (i in 1:nrow(cleavages)) {
    if (cleavages$N_terminus[i] == 'False') {
      cleavages$cleavage_loc[i] <- cleavages$start_seq[i]
      
    }
    
    else if (cleavages$C_terminus[i] == 'False')
      
      cleavages$cleavage_loc[i] <- cleavages$end_seq[i]
    
  }
  
  return(cleavages)
}
#Filters out proteins above cutoff and not being searched for
filterer <- function(df, cutoff, symbols) {
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
#Adds symbols to Accession column in dataframe
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
          as.data.frame(str_locate(fasta$seq[i], cleavages$Var1[j]))
        print(locate$start)
        print(is.na(locate))
        if(!is.na(locate)){
          cleavages$start_seq[j] <- locate$start
          cleavages$end_seq[j] <- locate$end
        }
        else{
       
        }
      }
    }
  }
  return(cleavages)
  
}

remove_false <- function(df) {
  return(subset(df,!grepl("XXX", df$Accession)))
}
quant <- function(proteins, x) {
  if (missing(x)) {
    return(subset(proteins, Freq > quantile(proteins$Freq, probs = 0.95)))
  } else{
    return(subset(proteins, Freq > quantile(proteins$Freq, probs = x)))
  }
}

top_n <- function(df, x) {
  df <- df[order(-df$Freq),]
  if (missing(x)) {
    return(df <- head(df, -(nrow(df) - 100)))
  } else{
    return(df <- head(df, -(nrow(df) - x)))
    
  }
}

library("tidyverse")
library("org.Hs.eg.db")
library("biomartr")
library("magrittr")
library("Biostrings")
library("dplyr")
library("tibble")
library("diffdf")

setwd("V:/Jason/HendriksFilesSemiTrytpic/CSV")
path <- "V:/Jason/HendriksFilesSemiTrytpic/CSV"
file.names <- list.files(path, pattern = ".csv")

ptm <- proc.time()


# 
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
# 
# 
# 
# #used for testing
# dr<- df
df<- dr


df <- remove_false(df)

ctr <- split_ctr(df)
ad <- split_ad(df, ctr)

ad <- as.data.frame(table(ad$Accession))
ad$Cohort <- c("AD")
ctr <- as.data.frame(table(ctr$Accession))
ctr$Cohort <- c("CTR")

percentile_ctr <- quant(ctr)
percentile_ad <- quant(ad)

top_ctr <- top_n(ctr)
top_ad <- top_n(ad)

fasta <- as.data.frame(readAAStringSet("proteins.fasta"))

fasta <- tibble::rownames_to_column(fasta, "accession")
colnames(fasta) <- c("accession", "seq")

uniprot <- c("P05067", "P02649", "P10636", "P10909", "P04156")

symbols <- select(org.Hs.eg.db, uniprot, "SYMBOL", "UNIPROT")



df <- filterer(df, 0.1, symbols)

 

df$N_terminus <- 'NULL'
df$C_terminus <- 'NULL'

df <- remove_tryp(df)
df <- find_semi(df)
ctr <- split_ctr(df)
ad <- split_ad(df, ctr)


seq_ad <- as.data.frame(table(ad$Sequence))
seq_ad$Cohort <- c("AD")
seq_ctr <- as.data.frame(table(ctr$Sequence))
seq_ctr$Cohort <- c("CTR")

seq_all <- (rbind(seq_ad, seq_ctr))
seq_all$Var1 <- as.character(seq_all$Var1)

start_seq <- vector("numeric", nrow(seq_all))
N_terminus <- vector("character", nrow(seq_all))
end_seq <- vector("numeric", nrow(seq_all))
C_terminus <- vector("character", nrow(seq_all))
protein <- vector("character", nrow(seq_all))
cleavage_loc <- vector("numeric", nrow(seq_all))


cleavages <-
  data.frame(protein,
             seq_all,
             cleavage_loc,
             start_seq,
             end_seq,
             N_terminus,
             C_terminus)

cleavages <- load_semi(cleavages, df)

cleavages <- load_protein(cleavages, df)

cleavages <- locate(cleavages, fasta)

cleavages <- cut_at(cleavages)



proc.time() - ptm
end
