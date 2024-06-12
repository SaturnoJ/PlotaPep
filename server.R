library(tidyr)
library(shiny)
library(shinyjs)
library(shinydashboard)
library("tidyverse")
library(org.Hs.eg.db) # (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("org.Hs.eg.db")
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

source("functions.R", local = T)


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
  j <- reactiveVal(1)
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
  colors_reactive <- reactiveVal()
  std <- reactiveVal()
  svg <- reactiveVal()
  less_than <- reactiveVal()
  file_type <- reactiveVal()
  #######
  #Fasta Input
  #######
  
  observeEvent(input$fastaInput, {
    fasta_reactive(readAAStringSet(input$fastaInput$datapath))
    
  })
  
  
  observeEvent(input$uniprotInput, {
    uniprotTemp <- input$uniprotInput
    uniprotTemp <- strsplit(uniprotTemp, ", |,| , ")
    uniprot(uniprotTemp)
  })
  
  
  #######
  
  #File Input
  ######
  
  
  
  
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
  
  observeEvent(input$file_typeInput, {
    file_type(input$file_typeInput)
    
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
  
  observeEvent(input$pInput, {
    p_cutoff(input$pInput)
    
  })
  
  observeEvent(input$stdInput, {
    std(input$stdInput)
    
  })
  
  observeEvent(input$svgInput, {
    svg(input$svgInput)
    
  })
  
  observeEvent(input$comparativeInput, {
    comparative(input$comparativeInput)
    
  })
  
  observeEvent(input$ctrInput, {
    name_temp <- input$ctrInput
    # name_temp <- strsplit(name_temp, ", |,| , ")
    ctr_name(unlist(name_temp))
    
  })
  
  observeEvent(input$cohortInput, {
    name_temp <- input$cohortInput
    name_temp <- strsplit(name_temp, ", |,| , ")
    cohort_name(unlist(name_temp))
    
  })
  
  observeEvent(input$colorsInput, {
    colors_temp <- input$colorsInput
    colors_temp <- strsplit(colors_temp, ", |,| , ")
    colors_reactive(unlist(colors_temp))
    
  })
  
  
  
  ########
  
  
  
  
  
  
  observeEvent(input$confirmFile, {
    #UI LOGIC
    # tryCatch({
      #######
      #Progress bar
      #######
      progress <- shiny::Progress$new()
      progress$set(message = "Preparing Data", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      updateProgress <- function(value = NULL,
                                 detail = NULL) {
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
      less_than_two <- data.frame()
      uniprot_ids <- as.data.frame(uniprot())
      colnames(uniprot_ids)  <- c("Protein.ID")
      protein_plot <- list()
      
      cohorts <- toupper(cohort_name())
      parametric_type <- parametric()
      ctr <- toupper(ctr_name())
      protein_id_loop <- list
      colors_df <- as.list(colors_reactive())
      fasta <- as.data.frame(fasta_reactive())
      fasta <- tibble::rownames_to_column(fasta, "accession")
      suppressWarnings(
        fasta %<>%
          separate(accession, sep = '\\|', c("sp", "Protein.ID", "info")) %>%
          separate(info, sep = 'GN=', c("info", "Gene")) %>%
          separate(Gene, sep = " ", c("Gene", "other"))
      )
      fasta <- fasta[-c(1, 5)]
      df <- files()
      #######
      #Dateframe creation
      #######


        for (i in df) {
          #Comparative
          #########
          for (j in 1:length(cohorts)) {
            temp <-
              cleanData(
                i,
                as.integer(intensity()),
                as.integer(lfq()),
                as.integer(cutoff()),
                as.integer(file_type()),
                updateProgress
              )
            temp[, 1] <- toupper(temp[, 1])
            if (!isEmpty(temp[temp$Sample %like% toupper(cohorts[j]),])) {
              cohort_df <-
                cohortSplit(temp,
                            ctr,
                            cohorts[j],
                            as.integer(std()),
                            updateProgress)
              
            }
            else{
              next
            }
            
            possible <-
              cohort_df %>%  spread(Cohort, Intensity) %>% group_by(Protein) %>% summarise(ctr = sum(!is.na(get("CTR"))),
                                                                                           cohort1 = sum(!is.na(get(
                                                                                             cohorts[j]
                                                                                           )))) %>%
              mutate(possible = ifelse(ctr < 2 |
                                         cohort1 < 2, FALSE, TRUE)) %>%
              filter(possible)
            
            cohort_df <-
              cohort_df %>% filter(Protein %in% possible$Protein)
            cohort_df$Intensity <- log2(cohort_df$Intensity)
            
            
            if (parametric_type == 0) {
              statisical_test <-
                getMannWhit(cohort_df, p_cutoff(), updateProgress)
              
              
              
            }
            else{
              statisical_test <- getTtest(cohort_df, p_cutoff(), updateProgress)
              
              
            }
            temp <-  anti_join(cohort_df, temp, by = "Protein")
            less_than_two <- rbind(less_than_two, temp)
            comparative_combined <-
              rbind(comparative_combined, statisical_test)
          }
        }
        
        if (is.function(updateProgress)) {
          text <- "Removing Outliers"
          updateProgress(detail = text)
        }
      file_names <- unique(comparative_combined$File)
      #Plotting loop
      ######
      if (comparative() == 0) {
        leftovers <-
          anti_join(statisical_test, fasta, by = "Protein.ID")
        
        statisical_test %<>% inner_join(fasta,
                                        by = c("Protein.ID", "Gene"),
                                        multiple = "all") %>% arrange(Protein.ID)
        
        located_peptides(locatePeptides(statisical_test, updateProgress))
        #######
        if (nrow(uniprot_ids) > 0) {
          protein_id_loop <- uniprot_ids
        }
        else{
          protein_id_loop <-  unique(statisical_test$Protein.ID)
        }
        for (i in protein_id_loop) {
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
            mutate(isSignificant = ifelse(filtered_results()$p.adj < p_cutoff(), "yes", "no")) %>%
            
            #Start plotting
            ggplot() +
            
            #We take here the 'mean' but this is of course X-times the same value
            scale_y_continuous(limits = c(NA, NA)) +
            xlim(0, protein_length + 1) +
            
            #Create some lines to help visualise the start and end of the protein.
            geom_vline(xintercept = 0,
                       lwd = 2,
                       alpha = .5) +
            geom_vline(xintercept = protein_length + 1,
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
              linewidth = .15,
            ) +
            
            #Set the Ggplot theme, limit the y-axis for FC.
            theme_bw() +
            
            #Specify the colours I want to use for the isSignificant column
            scale_color_manual(values = c("yes" = "black", "no" = "red")) +
            theme(legend.position = "bottom") +
            
            #x and yaxis titles
            xlab("Protein Sequence") +
            ylab("FC") +
            labs(
              subtitle = paste0(
                as.character(i),
                " p == ",
                p_cutoff(),
                " Standard Deviations ==",
                std()
              ),
              title = as.character(filtered_results()$Gene)
            )
          
          
          
          
        }
      }
      ######
      else{
        leftovers <-
          anti_join(comparative_combined, fasta, by =c('Protein.ID',"Gene"))
        comparative_combined %<>% inner_join(fasta, by = c('Protein.ID',"Gene"), multiple = "all") %>% arrange(Protein.ID) %>% mutate(
          colors = case_when(
            File == as.character(file_names[1]) ~  as.character(colors_df[1]),
            File == as.character(file_names[2]) ~ as.character(colors_df[2]),
            File == as.character(file_names[3]) ~ as.character(colors_df[3]),
            File == as.character(file_names[4]) ~ as.character(colors_df[4])
          )
        )
        
        
        located_peptides(locatePeptides(comparative_combined, updateProgress))
        colors <- distinct(comparative_combined, Cohort, colors)
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
          plotting_protein(
            plotting_protein() %>%
              
              #Lets make a column based on significance
              mutate(
                isSignificant = ifelse(filtered_results()$p.adj < p_cutoff(), "yes", "no")
              )
            %>% mutate(color = ifelse(
              isSignificant == "yes",
              darken(plotting_protein()$color, 0.4),
              lighten(plotting_protein()$color, 0.2)
            ))
            
            %>% mutate(
              File = ifelse(
                isSignificant == "yes",
                paste0(plotting_protein()$File, "_Significant"),
                plotting_protein()$File
              )
            )
           
          )
          
          # names(pal_temp) <- colors$Cohort
          # pal(pal_temp)
          protein_plot[[i]] = plotting_protein() %>%
            
            #Start plotting
            ggplot() +
            
            #We take here the 'mean' but this is of course X-times the same value
            scale_y_continuous(limits = c(NA, NA)) +
            xlim(0, protein_length + 1) +
            
            #Create some lines to help visualise the start and end of the protein.
            geom_vline(xintercept = 0,
                       lwd = 2,
                       alpha = .5) +
            geom_vline(xintercept = protein_length + 1,
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
                #fill = isSignificant,
                fill = color,
              ),
              alpha = 0.75,
              linewidth = .25,
              
            ) +
            
            #Set the Ggplot theme, limit the y-axis for FC.
            theme_bw() +
            
            #Specify the colours I want to use for the isSignificant column
            #scale_fill_manual(values = c("yes" = darken(plotting_protein()$color,0.4), "no" = lighten(plotting_protein()$color,0.2))) +
            scale_fill_identity("File",
                                labels = setNames(
                                  unique(plotting_protein()$File),
                                  unique(plotting_protein()$color)
                                ),
                                guide = "legend") +
            theme(legend.position = "bottom",
                  legend.key = element_rect(colour = "white")) +
            guides(fill = guide_legend(order = 1),
                   color = guide_legend(order = 2)) +
            #x and yaxis titles
            xlab("Protein Sequence") +
            ylab("FC") +
            
            labs(
              subtitle = paste0(
                as.character(i),
                " p == ",
                p_cutoff(),
                " standard deviations ==",
                std()
              ),
              title = as.character(filtered_results()$Gene)
            )
          
          
          
        }
        
        
        
      }
      #Reactive Variables
      ######
      protein_plotr(protein_plot)
      plot_fasta(fasta)
      if (nrow(uniprot_ids) == 0) {
        # browser()
        if (comparative() == 0) {
          temp_id <- as.data.frame(unique(statisical_test$Protein.ID))
          final_dataframe(located_peptides())
          
        }
        else{
          temp_id <- as.data.frame(unique(comparative_combined$Protein.ID))
          final_dataframe(located_peptides())
          
        }
        colnames(temp_id) <- c("Protein.ID")
        protein_ids(temp_id)
      }
      else{
        protein_ids(uniprot_ids)
      }
      removed_proteins(leftovers)
      less_than(less_than_two)
      final_dataframe(located_peptides())
      ######
      
      
      shinyjs::enable("previousPlot")
      shinyjs::enable("nextPlot")
      
      shinyjs::enable("downloadPlot")
      
    # },
    
    # warning = function(warn) {
    #   showNotification(paste0("Warning!Warning!Warning!"), type = 'warning')
    # },
    # error = function(err) {
    #   showNotification(paste0("Check input parameters. Ensure that the cohorts identified are represented in the dataset. Ensure that dataset type is correct for file and if more than one file is being used that comparative is selected."), type = 'err')
    # })
    # 
    
    
  })
  
  output$downloadPlot <- downloadHandler(
    filename = 'plots.zip',
    content = function(file) {
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      saved <- data.frame(Protein.ID = character())
      plots <- as.list(protein_plotr())
      if (as.integer(svg()) == 0) {
        for (i in unique(protein_ids()$Protein.ID)) {
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          if (nrow(plotting_protein()) > 0) {
            filtered_results(filter(located_peptides(), Protein.ID == i))
            filtered_results(inner_join(
              plot_fasta(),
              plotting_protein(),
              by  = c('Protein.ID', 'x', 'Gene', 'info'),
              multiple = "all"
            ))
            
            protein_length <- unique(nchar(filtered_results()$x))
            
            
            ggsave(
              paste(
                as.character(filtered_results()$Gene[1]),
                "_",
                as.character(filtered_results()$Protein.ID[1]),
                ".png",
                sep = ""
              ),
              plot = plots[[i]],
              device = "png"
            )
            
          }
        }
      }
      
      else{
        for (i in unique(protein_ids()$Protein.ID)) {
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          if (nrow(plotting_protein()) > 0) {
            filtered_results(filter(located_peptides(), Protein.ID == i))
            filtered_results(inner_join(
              plot_fasta(),
              plotting_protein(),
              by  = c('Protein.ID', 'x', 'Gene', 'info'),
              multiple = "all"
            ))
            
            protein_length <- unique(nchar(filtered_results()$x))
            
            
            ggsave(
              paste(
                as.character(filtered_results()$Gene[1]),
                "_",
                as.character(filtered_results()$Protein.ID[1]),
                ".svg",
                sep = ""
              ),
              plot = plots[[i]],
              device = "svg"
            )
          }
        }
        
      }
      
      write.csv(final_dataframe()
                , file = "final_dataframe.csv"
                , row.names = F)
      write.csv(removed_proteins()
                ,
                file = "removed_proteins_fasta.csv"
                ,
                row.names = F)
      write.csv(less_than()
                ,
                file = "removed_proteins_less_two.csv"
                ,
                row.names = F)
      
      zip_files <-
        list.files(path = getwd(), pattern = ".svg$|.csv$|.png$")
      
      zip::zip(file, files = zip_files)
      files_to_delete <-
        dir(path = getwd() , pattern = "*.png$|*.svg$|*.csv$")
      file.remove(file.path(getwd(), files_to_delete))
      
    }
    
  )
  
  
  
  
  observeEvent(input$nextPlot, {
    plots <- as.list(protein_plotr())
    uniprot_ids <- as.data.frame(uniprot())
    colnames(uniprot_ids)  <- c("Protein.ID")
    if (nrow(uniprot_ids) > 0) {
      if (j() < nrow(uniprot_ids)) {
        output$plot <- renderUI({
          print(j())
          i <-  unique(uniprot_ids[j(), 'Protein.ID'])
          
          
          print(i)
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          filtered_results(filter(located_peptides(), Protein.ID == i))
          filtered_results(inner_join(
            plot_fasta(),
            plotting_protein(),
            by  = c('Protein.ID', 'x', 'Gene', 'info'),
            multiple = "all"
          ))
          
          
          protein_length <- unique(nchar(filtered_results()$x))
          
          
          renderPlot({
            plots[[i]]
            
          })
          
        })
        j(j() + 1)
        print("addd")
        print((j))
      }
      else {
        output$plot <- renderUI({
          i <-  unique(uniprot_ids[j(), 'Protein.ID'])
          print(i)
          print(j())
          
          
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          filtered_results(filter(located_peptides(), Protein.ID == i))
          filtered_results(inner_join(
            plot_fasta(),
            plotting_protein(),
            by  = c('Protein.ID', 'x', 'Gene', 'info'),
            multiple = "all"
          ))
          
          protein_length <- unique(nchar(filtered_results()$x))
          
          
          renderPlot({
            plots[[i]]
          })
          
        })
        j(1)
        
      }
    }
    else{
      if (j() < nrow(protein_ids())) {
        output$plot <- renderUI({
          i <-  unique(protein_ids()[j(), 'Protein.ID'])
          
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          filtered_results(filter(located_peptides(), Protein.ID == i))
          filtered_results(inner_join(
            plot_fasta(),
            plotting_protein(),
            by  = c('Protein.ID', 'x', 'Gene', 'info'),
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
          i <-  unique(protein_ids()[j(), 'Protein.ID'])
          
          
          
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          filtered_results(filter(located_peptides(), Protein.ID == i))
          filtered_results(inner_join(
            plot_fasta(),
            plotting_protein(),
            by  = c('Protein.ID', 'x', 'Gene', 'info'),
            multiple = "all"
          ))
          
          protein_length <- unique(nchar(filtered_results()$x))
          
          
          renderPlot({
            plots[[i]]
          })
          
        })
        j(1)
        
      }
    }
    
    
  })
  
  observeEvent(input$previousPlot, {
    plots <- as.list(protein_plotr())
    uniprot_ids <- as.data.frame(uniprot())
    colnames(uniprot_ids)  <- c("Protein.ID")
    if (nrow(uniprot_ids) > 0) {
      if (j() < nrow(uniprot_ids)  & j() != 1) {
        output$plot <- renderUI({
          print(j())
          i <-  unique(uniprot_ids[j(), 'Protein.ID'])
          
          
          print(i)
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          filtered_results(filter(located_peptides(), Protein.ID == i))
          filtered_results(inner_join(
            plot_fasta(),
            plotting_protein(),
            by  = c('Protein.ID', 'x', 'Gene', 'info'),
            multiple = "all"
          ))
          
          
          protein_length <- unique(nchar(filtered_results()$x))
          
          
          renderPlot({
            plots[[i]]
            
          })
          
        })
        j(j() - 1)
        
      }
      else {
        output$plot <- renderUI({
          i <-  unique(uniprot_ids[j(), 'Protein.ID'])
          print(i)
          print(j())
          
          
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          filtered_results(filter(located_peptides(), Protein.ID == i))
          filtered_results(inner_join(
            plot_fasta(),
            plotting_protein(),
            by  = c('Protein.ID', 'x', 'Gene', 'info'),
            multiple = "all"
          ))
          
          protein_length <- unique(nchar(filtered_results()$x))
          
          
          renderPlot({
            plots[[i]]
          })
          
        })
        j(1)
        
      }
    }
    else{
      if (j() <= nrow(protein_ids()) & j() != 1) {
        output$plot <- renderUI({
          i <-  unique(protein_ids()[j(), 'Protein.ID'])
          
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          filtered_results(filter(located_peptides(), Protein.ID == i))
          filtered_results(inner_join(
            plot_fasta(),
            plotting_protein(),
            by  = c('Protein.ID', 'x', 'Gene', 'info'),
            multiple = "all"
          ))
          
          
          protein_length <- unique(nchar(filtered_results()$x))
          
          
          renderPlot({
            plots[[i]]
            
          })
          
        })
        
        j(j() - 1)
        
      }
      else {
        output$plot <- renderUI({
          i <-  unique(protein_ids()[j(), 'Protein.ID'])
          
          
          
          plotting_protein(filter(located_peptides(), Protein.ID == i))
          filtered_results(filter(located_peptides(), Protein.ID == i))
          filtered_results(inner_join(
            plot_fasta(),
            plotting_protein(),
            by  = c('Protein.ID', 'x', 'Gene', 'info'),
            multiple = "all"
          ))
          
          protein_length <- unique(nchar(filtered_results()$x))
          
          
          renderPlot({
            plots[[i]]
          })
          
        })
        j(as.integer(nrow(protein_ids())))
        
      }
    }
    
    
  })
  
  
  
  
}