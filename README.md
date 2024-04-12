# PlotaPeP
For use with MSfragger peptide outputs. Upload a fasta files with proteins you wish to plot, upload fragger data, identify cohorts as labled in fragger out put, select colors, std, modifications, run.

## Installation

Uses Biostrings, tidyverse, and org.Hs.eg.db R packages. 

```R

    install.packages("BiocManager")
    BiocManager::install("Biostrings")
    BiocManager::install("org.Hs.eg.db")
    install.packages("tidyverse")
    install.packages("shiny")
install.packages("shinyBS")
install.packages("shinydisconnect")
install.packages("shinyjs")
install.packages("shinydashboard")

```
