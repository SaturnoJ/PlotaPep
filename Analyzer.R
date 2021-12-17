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
load_semi <- function(cleavage_x, df) {
  for (i in 1:nrow(cleavage_x)) {
    for (j in 1:nrow(df)) {
      if (cleavage_x$Var1[i] == df$Sequence[j]) {
        cleavage_x$N_terminus[i] <- df$N_terminus[j]
        cleavage_x$C_terminus[i] <- df$C_terminus[j]
        break
      }
      
    }
    
  }
  return(cleavage_x)
}
#Loads protein symbols into aggregate dataset
load_protein <- function(cleavage_x, df) {
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(cleavage_x)) {
      if (df$Sequence[i] == cleavage_x$Var1[j] &&
          !is.null(cleavage_x$protein[j])) {
        cleavage_x$protein[j] <- df$Accession[i]
        
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

locate <- function(cleavage_x, fasta) {
  for (i in 1:nrow(fasta)) {
    for (j in 1:nrow(cleavage_x)) {
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
ptm <- proc.time()

setwd("V:/Jason/HendriksFilesSemiTrytpic/CSV")
path <- "V:/Jason/HendriksFilesSemiTrytpic/CSV"
file.names <- list.files(path, pattern = ".csv")

file <-
  read_delim(
    file.names[1],
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE
  )
file$file <- file.names[1]
df <- file

for (i in 2:length(file.names)) {
  file <-
    read_delim(
      file.names[i],
      delim = "\t",
      escape_double = FALSE,
      trim_ws = TRUE,
      show_col_types = FALSE
    )
  file$file <- file.names[i]
  df <- rbind(df, file)
}

df <- remove_false(df)


dr<- df
df<- dr

ctr <- split_ctr(df)
ad <- split_ad(df, ctr)

ad <- as.data.frame(table(ad$Accession))
ad$Cohort <- c("AD")
ctr <- as.data.frame(table(ctr$Accession))
ctr$Cohort <- c("CTR")

percentile_ctr <- quant(ctr)
percentile_ad <- quant(ad)

ctr_top <- top_n(ctr)
ad_top <- top_n(ad)

fasta <- as.data.frame(readAAStringSet("protiens.fasta"))

fasta <- tibble::rownames_to_column(fasta, "accession")
colnames(fasta) <- c("accession", "seq")

uniprot <- c("P05067", "P02649", "P10636", "P10909", "P04156")

symbols <- select(org.Hs.eg.db, uniprot, "SYMBOL", "UNIPROT")



df <- filterer(df, 0.1, symbols)

df <- add_symbols(df, symbols)


df$N_terminus <- 'NULL'
df$C_terminus <- 'NULL'

df <- remove_tryp(df)
df <- find_semi(df)
ctr <- split_ctr(df)
ad <- split_ad(df, ctr)


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
cleavage_loc <- vector("numeric", nrow(all_seq))


cleavage_x <-
  data.frame(protein,
             all_seq,
             cleavage_loc,
             start_seq,
             end_seq,
             N_terminus,
             C_terminus)

cleavage_x <- load_semi(cleavage_x, df)

cleavage_x <- load_protein(cleavage_x, df)

cleavage_x <- locate(cleavage_x, fasta)

cleavage_x <- cut_at(cleavage_x)



proc.time() - ptm
end
