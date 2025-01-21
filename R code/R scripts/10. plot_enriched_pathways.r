# Filter and plot enriched pathways 

# input data
# - .tsv files with enriched pathways
# - top_match_between_IMC_to_RNA_clusters_original.csv (generated in 9. RNAseq_IMC_clusters_matching.R)

#### 0. Load packages #####

#install.packages("devtools")
library(devtools)
#BiocManager::install(version = "3.14")
#BiocManager::install("ComplexHeatmap")
#install_github('jokergoo/ComplexHeatmap') # does not work
library(ComplexHeatmap)
library(Matrix)
library(data.table)
library(stringr)
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)

setwd("~/Dropbox (UK Dementia Research Institute)/Alessia Lab Dropbox/Projects/UKDRI project/Image analysis/HistoCat:SIMPLI analysis/051022 TREM2 cohort/R code")

#### 1. define folder to analyse and comparisons >> OPEN << ####

# define which input folder to analyse
subfolder <- "all_samples_comparing_clusters"
subfolder <- "CTRL_vs_AD_revision"
subfolder <- "CTRL_vs_AD"

# define comparisons
if_comparison <- TRUE # put to FALSE if there are no 
comparison <- c("amyloid_beta", "diagnosis", "PHF1")

#### 2. create output folders ####
path_plots <- "plots/RNAseq/pathway_analysis/"
path_tables <- "tables/RNAseq/pathway_analysis/"

dir.create(path_tables)
dir.create(path_plots)

dir.create(paste0(path_tables, subfolder))
dir.create(paste0(path_plots, subfolder))

input_path <- paste0("input/pathway_analysis/", subfolder)

#### 3. filter and plot pathways #####

# pathways to exclude
text_out <- "epithelial|lens fiber cell|erythrocyte|crest|duct|aortic|barrier|granulocyte|osteoblast|ruffle|vascular|angiogenesis|cholerae|pulmonary|myeloid cell|endothelial|fat cell|ventricular|Taste|Insulin|Salivary|Bacterial|Legionellosis|Leishmaniasis|leukemia|glioma|histamine|platelet|cardiac|positive|negative|angiotensin|morphogenesis|renal|alcoholism|cardio|gastric"

if (if_comparison){
  for (comp in 1:length(comparison)){
    ##### 3.1 Read in and filter pathways ####
    
    # get file list and open in massive multiple dataframe (large list)
    file_list <- list.files(input_path, 
                            recursive = TRUE, 
                            pattern = "tsv", 
                            full.names = T)
    
    file_list <- file_list[grep(comparison[comp], file_list)]
    
    celltype <- gsub('^[^\\/]*\\/([^\\/]*)', "", file_list)
    celltype <- sub('.', '', celltype)
    celltype <- sub("/.*", "", celltype)
    
    dt_list<- lapply(file_list, function(file){
      dt_l <- list()
      for(j in file){
        #j=file_list[1]
        file_name <- tools::file_path_sans_ext(basename(j))
        dt_l[[file_name]] <- read.delim(j) %>% dplyr::mutate(database = file_name)
      }
      dt <- do.call(rbind, dt_l) 
    })
    
    names(dt_list) <- file_list
    
    dt_list_mod <- lapply(names(dt_list), function(x){
      dt <- dt_list[[x]]
      dt <- dt %>%
        mutate(direction = purrr::map_chr(x,
                                          ~strsplit(., "/")[[1]][4]) %>% gsub(sprintf("pathway_analysis_%s_", comparison[comp]), "", .),
               celltype = purrr::map_chr(x,
                                         ~strsplit(., "/")[[1]][5]),
               trem2 = purrr::map_chr(x,
                                      ~strsplit(., "/")[[1]][6])
        )
    } )
    
    # filter pathways for FDR, overlap and odds_ratio
    dt_list_select <- lapply(dt_list_mod, function(dt) {dt <- dt %>%
      dplyr::filter(!duplicated(tolower(description))) %>%
      dplyr::filter(FDR < 0.04, 
                    overlap >= 3,
                    odds_ratio > 8) %>%
      dplyr::filter(!grepl(text_out, description, ignore.case = TRUE)) %>%
      dplyr::filter(!duplicated(clusters)) %>%
      as.data.frame()})
    
    dt_list_select <- do.call(rbind, dt_list_select) 
    
    dotplot <- dt_list_select %>%
      mutate(log10_fdr = ifelse(direction == "Up", -log10(FDR), log10(FDR)))
    
    # keep only RNA clusters matched with IMC clusters
    candidates <- read.csv("tables/RNAseq/clusters_matching/top_match_between_IMC_to_RNA_clusters_original.csv")
    candidates <- unique(candidates$RNA_cluster)
    
    dotplot <- dotplot %>% 
      subset(dotplot$celltype %in% candidates) %>%
      mutate(comparison = comparison[comp])
    
    dotplot <- dotplot %>% 
      subset(dotplot$trem2 %in% c("CV", "TREM2var"))
    
    write.csv(dotplot, paste0(path_tables, subfolder, "/pathway_analysis_", comparison[comp], "_matched_cluster.csv"))
    
    
    ##### 3.2 Plot ####
    
    dat <- dotplot
    
    p <- ggplot(dat, 
                aes(x = trem2, 
                    y = description,
                    color = direction, 
                    #fill = col_padj,
                    alpha = FDR,
                    size = odds_ratio
                    #shape = group
                )
    ) + 
      geom_point() +
      facet_wrap(~celltype,
                 # scales = "free_y",
                 ncol = 14
      ) +
      xlab(label = "") +
      ylab(label = "") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                       size = 10, color = "black"),
            strip.text = element_text(size = 12, hjust = -0,
                                      angle = 90,
                                      color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12, color = "black"),
            strip.placement = "outside",
            strip.background = element_blank(),
            axis.ticks.x = element_blank(), 
            axis.ticks.y = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_line(colour = "lightgrey", linewidth = 0.1),
            legend.text = element_text(size = 12, color = "black"),
            legend.title = element_text(size = 12, color = "black"),
      ) +
      scale_color_manual(values = c("Up" = "red",
                                    "Down" ="blue")
      ) +
      scale_alpha(range = c(1, 0.1))
    
    p 
    
    pdf(paste0(path_plots, subfolder, "/pathways_in_", comparison[comp], ".pdf"),
        height = 30, 
        width = 14)
    print(p)
    dev.off()
    
  }
  #### 3.3. merge comparisons and plot #####
  
  # merge dataframes into one
  # get file list and open in massive multiple dataframe (large list)
  file_list_new <- list.files(paste0(path_tables, subfolder), 
                              recursive = TRUE, 
                              pattern = "csv", 
                              full.names = T)
  
  # Initialize an empty list to store data frames
  data_frames <- list()
  
  # Loop through each filename, read the data frame, add the filename as a new column, and store in the list
  for(i in seq_along(file_list_new)) {
    # Read the data frame (assuming CSV files; adjust read function as needed)
    df <- read.csv(file_list_new[i])
    
    # Store the modified data frame in the list
    data_frames[[i]] <- df
  }
  
  # Merge all data frames into a single data frame
  merged_df <- bind_rows(data_frames)
  
  # View the merged data frame
  write.csv(merged_df, paste0(path_tables, subfolder, "/pathway_analysis_combined_matched_cluster.csv"))
  
  # Plot
  
  p <- ggplot() +
    geom_point(data = merged_df, 
               aes(x = comparison, 
                   y = description,
                   color = direction, 
                   #fill = col_padj,
                   alpha = FDR,
                   size = odds_ratio
                   #shape = group
               )
    ) + 
    scale_y_discrete(drop=FALSE) +
    facet_grid(~celltype
               # scales = "free_y",
               #ncol = 14
    ) +
    xlab(label = "") +
    ylab(label = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                     size = 10, color = "black"),
          strip.text = element_text(size = 12, hjust = -0,
                                    angle = 90,
                                    color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          strip.placement = "outside",
          strip.background = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_line(colour = "lightgrey", linewidth = 0.1),
          legend.text = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
    ) +
    scale_color_manual(values = c("Up" = "red",
                                  "Down" ="blue")
    ) +
    scale_alpha(range = c(1, 0.1))
  
  
  p$layers = rev(p$layers)
  
  pdf(paste0(path_plots, subfolder, "/merged_pathways.pdf"), 
      height = 35, 
      width = 14)
  print(p)
  dev.off()
  
} else {
  
  ##### 3.4 Read in and filter pathways ####
  
  # get file list and open in massive multiple dataframe (large list)
  file_list <- list.files(input_path, 
                          recursive = TRUE, 
                          pattern = "tsv", 
                          full.names = T)
  
  celltype <- gsub('.tsv', "", file_list)
  celltype <- sapply(strsplit(celltype, "/"), "[", 4)
  
  dt_list<- lapply(file_list, function(file){
    dt_l <- list()
    for(j in file){
      file_name <- tools::file_path_sans_ext(basename(j))
      dt_l[[file_name]] <- read.delim(j) %>% dplyr::mutate(database = file_name)
    }
    dt <- do.call(rbind, dt_l) 
  })
  
  names(dt_list) <- file_list
  
  dt_list_mod <- lapply(names(dt_list), function(x){
    dt <- dt_list[[x]]
    dt <- dt %>%
      mutate(celltype = purrr::map_chr(x,
                                       ~strsplit(., "/")[[1]][4])
      )
  } )
  
  # filter pathways for FDR, overlap and odds_ratio
  dt_list_select <- lapply(dt_list_mod, function(dt) {dt <- dt %>%
    dplyr::filter(!duplicated(tolower(description))) %>%
    filter(FDR <= 0.05, 
           overlap >= 3, 
           odds_ratio > 8) %>%
    dplyr::filter(!grepl(text_out, description, ignore.case = TRUE)) %>%
    filter(!duplicated(clusters)) %>%
    # dplyr::mutate(celltype = x) %>%
    as.data.frame()})
  
  dt_list_select <- do.call(rbind, dt_list_select) 
  
  dotplot <- dt_list_select
  
  # keep only RNA clusters matched with IMC
  candidates <- read.csv("tables/RNAseq/clusters_matching/top_match_between_IMC_to_RNA_clusters_original.csv")
  candidates <- unique(candidates$RNA_cluster)
  
  dotplot <- dotplot %>% 
    subset(dotplot$celltype %in% candidates)
  
  write.csv(dotplot, paste0(path_tables, subfolder, "/pathway_analysis_only_matched_cluster.csv"))
  
  ##### 3.5 Plot ####
  
  dat <- dotplot
  
  p <- ggplot() +
    geom_point(data = dat, 
               aes(x = celltype, 
                   y = description,
                   #color = direction, 
                   #fill = col_padj,
                   alpha = FDR,
                   size = odds_ratio
                   #shape = group
               )
    ) + 
    scale_y_discrete(drop=FALSE) +
    xlab(label = "") +
    ylab(label = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                     size = 10, color = "black"),
          strip.text = element_text(size = 12, hjust = -0,
                                    angle = 90,
                                    color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          strip.placement = "outside",
          strip.background = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_line(colour = "lightgrey", linewidth = 0.1),
          legend.text = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
    ) +
    scale_color_manual(values = c("Up" = "red",
                                  "Down" ="blue")
    ) +
    scale_alpha(range = c(1, 0.1))
  
  p$layers = rev(p$layers)
  
  pdf(paste0(path_plots, subfolder, "/merged_pathways.pdf"), 
      height = 30, 
      width = 14)
  print(p)
  dev.off()
}












