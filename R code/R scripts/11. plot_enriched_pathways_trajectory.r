# Filter and plot enriched pathways from trajectory analysis

# input data
# - .tsv files with enriched pathways
# - top_match_between_IMC_to_RNA_clusters_original.csv (generated in 9. RNAseq_IMC_clusters_matching.R)

#### 0. Load packages, create folders #####

library(dplyr)
#install.packages("devtools")
#library(devtools)
#install_github('jokergoo/ComplexHeatmap')
library(ComplexHeatmap)
library(Matrix)
library(data.table)
library(stringr)
library(ggplot2)
library(tidyr)
library(plyr)
library(vctrs)
library(ggpubr)
library(scales)
library(viridis)

setwd("~/R code/")

# create output folders
subfolder <- "trajectory_analysis"

path_plots <- "plots/RNAseq/pathway_analysis/"
path_tables <- "tables/RNAseq/pathway_analysis/"

dir.create(path_tables)
dir.create(path_plots)

dir.create(paste0(path_tables, subfolder))
dir.create(paste0(path_plots, subfolder))

input_path <- paste0("input/pathway_analysis/", subfolder)

#### 1. filter and plot pathways ####

# pathways to exclude
text_out <- "Sarcolemma|ossification|PC12|epithelial|lens fiber cell|erythrocyte|crest|duct|aortic|barrier|granulocyte|osteoblast|ruffle|vascular|angiogenesis|cholerae|pulmonary|myeloid cell|endothelial|fat cell|ventricular|Taste|Insulin|Salivary|Bacterial|Legionellosis|Leishmaniasis|leukemia|glioma|histamine|platelet|cardiac|positive|negative|angiotensin|morphogenesis|renal|alcoholism|cardio|gastric"

# find all cluster names
clusters <- list.dirs(path = input_path, full.names = TRUE, recursive = FALSE)
clusters <- gsub(paste0(input_path, "/trajectory_"), "", clusters)

for (clu in 1:length(clusters)){
  #clu=2
  
  # create folders
  dir.create(paste0(path_tables, subfolder, '/', clusters[clu]))
  dir.create(paste0(path_plots, subfolder, '/', clusters[clu]))
  
  ##### 3.1 Read in and filter pathways ####
  
  # get file list and open in massive multiple dataframe (large list)
  file_list <- list.files(path = paste0(input_path, "/trajectory_", clusters[clu]), 
                          recursive = TRUE, 
                          pattern = "tsv", 
                          full.names = T)
  
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
                                        ~strsplit(., "/")[[1]][4]) %>% gsub("trajectory_", "", .),
             module = purrr::map_chr(x,
                                       ~strsplit(., "/")[[1]][5])
      )
  } )
  
  # select only pathways with FDR < 0.05 and overlap > 4
  dt_list_select <- lapply(dt_list_mod, function(dt) {dt <- dt %>%
    dplyr::filter(!duplicated(tolower(description))) %>%
    filter(FDR < 0.05,
           overlap >= 3, 
           odds_ratio > 8
           ) %>%
    dplyr::filter(!grepl(text_out, description, ignore.case = TRUE)) %>%
    filter(!duplicated(clusters)) %>%
    as.data.frame()})
  
  dt_list_select <- do.call(rbind, dt_list_select) 
  
  write.csv(dt_list_select, paste0(path_tables, subfolder, "/", clusters[clu], "/trajectory_analysis_", clusters[clu], ".csv"))
  
  
  ##### 3.2 Plot####
  
  dat <- dt_list_select
  
  # define module order in the plot
  
  if (clusters[clu] == "Inh-L3-4-PVALB-HOMER3"){
    module_order <- c("Module4", "Module6", "Module2", "Module1", "Module3", "Module5")
  }
  
  if (clusters[clu] == "Exc-L5-RORB-LINC01202"){
    module_order <- c("Module3", "Module4", "Module6", "Module7", "Module2", "Module1", "Module5")
  }
  
  if (clusters[clu] == "Exc-L6-THEMIS-LINC00343"){
    module_order <- c("Module1", "Module2", "Module3", "Module4")
  }
  
  if (clusters[clu] == "Exc-L4-6-RORB-LCN15"){
    module_order <- c("Module2", "Module6", "Module5", "Module3", "Module4", "Module1")
  }
  
  # make module category the level factor for x axis
  dat <- dat %>% mutate(
    module = factor(module, levels = module_order)
  )

  # define order of pathways for plot
  order <- list(
    synapses = unique(dat$description[grep("synapses|synaptic|spine|presynapse|postsynaptic|dendrite|synapse", dat$description, ignore.case=TRUE)]),
    neurotransmitter = unique(dat$description[grep("Calcium|channel|receptor|neurotransmitter|potentiation|depression| ion |chloride", dat$description, ignore.case=TRUE)]),
    #secondary_modification = unique(dat$description[grep("Mannosylation|SUMO|glycosylation|phosphorylation|dephosphorylation|kinase|tyrosine", dat$description, ignore.case=TRUE)]),
    stress = unique(dat$description[grep("Homologous Recombination|Bcl-2|apoptosis|Double−Strand Break|Interleukin-6|granule|oxygen|stress|death|amyloid|Longevity|apoptotic", dat$description, ignore.case=TRUE)]), 
    metabolism = unique(dat$description[grep("Gluconeogenesis|Glycolytic|Proton|glycosylation|Sialyltransferase|glucose|metabolic|biosynthesis|metabolism|aerobic|ATP|mitochondrial|mitochondrion|glycolysis|alcohol", dat$description, ignore.case=TRUE)]),
    adhesion = unique(dat$description[grep("Migration|Netrin|Glycosaminoglycan|Proteoglycan|Glycosphingolipid|catenin|Oligosaccharide|Nerve Growth|Semaphorin|cell-cell|adhesion|matrix|junction|microtubule|Axon|actin|axonogenesis|projection|invagination", dat$description, ignore.case=TRUE)]),
    Ca_regulation = unique(dat$description[grep("cAMP|cGMP|Phosphodiesterase", dat$description, ignore.case=TRUE)]),
    autophagy = unique(dat$description[grep("Tau Protein Binding|MHC Class II Protein Complex Binding|Sema3A/PAK-dependent axon repulsion|Tau-Protein|Ubiquitination|Autoubiquitination|Macroautophagy|autophagy|Phagocytic|Unfolded", dat$description, ignore.case=TRUE)])
  )
  
  # add all remaining pathways
  order$other <- setdiff(unique(dat$description), unique(unlist(lapply(order, rbind))))

  # remove empty lists
  order <- Filter(function(k) length(k)>0, order)
  
  # turn into dataframe
  order2 <- plyr::ldply(order, cbind)
  names(order2) <- c("category", "pathway")
  
  # remove duplicated pathways
  order2 <- order2[!duplicated(order2$pathway), ]
  
  # add pathway category to dat
  dat$category <- order2$category[match(dat$description, order2$pathway)]
  
  # filter for out "others" pathways
  dat <- dat[!dat$category == "other",]
  
  # add pseudotime
  module_order_dat <- as.data.frame(module_order)
  module_order_dat$pseudotime <- rownames(module_order_dat)
  
  dat$pseudotime <- module_order_dat$pseudotime[match(dat$module, module_order_dat$module_order)]
  
  # define color palette
  COLS <- viridis(nrow(module_order_dat), option = "B") #“magma” (or “A”), “inferno” (or “B”), “plasma” (or “C”), “viridis” 
  
  # plot
  p <- ggplot() +
    geom_point(data = dat, 
               aes(x = module, 
                   y = description,
                   color = module, 
                   #fill = module,
                   alpha = FDR,
                   size = odds_ratio
                   #shape = group
               ), 
               stat = "identity"
    ) + 
    scale_color_manual(values=COLS) + 
    scale_y_discrete(drop=FALSE) +
    xlab(label = "") +
    ylab(label = "") +
    theme_bw() +
    guides(color = "none") + # remove the module legend
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
    scale_alpha(range = c(1, 0.1))

  p2 <- ggplot(dat, 
                  aes(x=module, y=1, fill = pseudotime)
                  ) +
    geom_tile() +
    labs(y = "Pseudotime") +
    theme_void() +
    scale_fill_manual(values=COLS) +
    theme(legend.position = 'none', 
          axis.title.y=element_text(angle = 0, hjust = 1, size = 20)
    )
    
  blank <- ggplot() + geom_blank() + theme_void()
  
  p2_half_size <- ggarrange(
    blank, p2,
    nrow = 1,
    ncol = 2,
    common.legend = FALSE,
    widths = c(8, 4)
  )
  
  p3 <- ggarrange(p2_half_size, 
                      p,
                      nrow = 2,
                      #hjust = 10,
                      # vjust = 1.5,
                      heights = c(0.2, 10)
                      #widths = c(2),
                      #align = 'v'
                      )
  
  width_plot <- 10
  height_plot <- 8

  pdf(paste0(path_plots, subfolder, "/", clusters[clu], "/trajectory_analysis_", clusters[clu], ".pdf"), 
      height = height_plot, 
      width = width_plot)
  print(p3)
  dev.off()
  
}
  
#### 2. merge pathway dataframe into one #####

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

# change names of lists in list
names(data_frames) <- file_list_new
names(data_frames) <- lapply(strsplit(names(data_frames), "/"), `[`, 5)

# Merge all data frames into a single data frame
merged_df <- bind_rows(data_frames)

# View the merged data frame
write.csv(merged_df, paste0(path_tables, subfolder, "/trajectory_analysis_combined_cluster.csv"))

  
  
  