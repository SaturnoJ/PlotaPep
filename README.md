# Tally

name is subject to change.
Tally is a proteolytic cleavage analysis tool. Tally takes in matched spectra from CSV files and determines what specific peptide cleavages occured in that spectra analysis.

## Installation

Uses Biostrings, tidyverse, and org.Hs.eg.db R packages. 

```R

    install.packages("BiocManager")
    BiocManager::install("Biostrings")
    BiocManager::install("org.Hs.eg.db")
    install.packages("tidyverse")
library("Biostrings")
library("org.Hs.eg.db")
library("tidyverse")

```
