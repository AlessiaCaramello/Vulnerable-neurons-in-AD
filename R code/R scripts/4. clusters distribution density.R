# Distribution of cells from single clusters

# input data:
# - clustered_cells.csv (from SIMPLI)
# - vars_per_sample.csv (metadata)
# - clusters_labels.csv (metadata)

##### 0. Load libraries ######

#install.packages("tidyverse")
library(tidyverse)
library(reshape)
library(data.table)
library(dplyr)
library(stringr)
library(tibble)
#install.packages("progress")
library(progress)
library(ggplot2)
library(ggsignif)
#install.packages("ggpubr")
library(ggpubr)
library(rstatix)
library(plyr)
#install.packages("DirichletReg")
library(DirichletReg)
library(reshape2)
library(pheatmap)
#install.packages("writexl")
library(writexl)
library(readxl)
library(openxlsx)

# set working directory
setwd("~/R code")

# create output folder
dir.create("plots/clusters/distribution")
dir.create("tables/clusters/distribution")

##### 1. Import files and define cluster resolution to use #####

# import files
clustered_cells <- data.table::fread("input/clustered_cells.csv")

# define resolution to analyse
res <- 'res_0.9_ids'

##### 2. Import/Calculate size/proportion of each layer in the cortex, and max cortex length ######

# import file to calculate layer proportion
layers_sample_replicates <- read_xlsx("tables/clusters/distribution/layers_proportion_for_sample_replicates.xlsx")

# otherwise run the following to generate file layers_sample_replicates

# ## Calculate the proportion of each layer for our acquired image
# 
# # get ROIs names
# sample_replicates <- unique(clustered_cells_new$Metadata_sample_name)
# #cells_position_clean_new <- cell_position
# 
# # Import reference layer thickness
# layers_proportion <- read_xlsx("input/layers_proportions.xlsx")
# 
# # make a table for layers length for each sample replicate
# layers_sample_replicates <- data.frame("sample_replicates" = sample_replicates,
#                                        "sample_names" = stringr::str_extract(sample_replicates, "[^_]*"),
#                                        "max_length" = 0,
#                                        "layer_1" = 0,
#                                        "layer_2" = 0,
#                                        "layer_3" = 0,
#                                        "layer_4" = 0,
#                                        "layer_5" = 0,
#                                        "layer_6" = 0,
#                                        "max_width" = 0)
# 
# # calculate max length for each sample replicate
# for (i in 1:length(sample_replicates)){
#   # extract rows of sample replicate to process
#   replicate_name <- sample_replicates[i]
# 
#   sample_rows <- clustered_cells_new %>% filter(
#     clustered_cells_new$Metadata_sample_name == replicate_name)
# 
#   # calculate max Y and insert in layers_sample_replicates.csv
#   layers_sample_replicates$max_length[i] = max(sample_rows$Y)
# 
#   # calculate max X in layers_sample_replicates.csv
#   layers_sample_replicates$max_width[i] = max(sample_rows$X)
# }
# 
# # check that max length extracted matches max length of original file
# max(layers_sample_replicates$max_length) == max(clustered_cells_new$Y)
# 
# # calculate size of each layer in the experiment
# for (rep in 1:nrow(layers_sample_replicates)){
#   for (layer in 1:nrow(layers_proportion)){
#     if (layer == 1) {
# 
#       # calculate end of layer 1
#       layer_perc_exp <- ((layers_sample_replicates$max_length[rep]/100)
#                          *layers_proportion[[layer,2]])
# 
#       # add layer length to layers_sample_replicates
#       layers_sample_replicates[[rep, layer+3]] <- layer_perc_exp
#     }
#     if (layer > 1) {
#       # calculate end of layers 2 to 6
#       layer_perc_exp <- layer_perc_exp + ((layers_sample_replicates$max_length[rep]/100)
#                                           *layers_proportion[[layer,2]])
#       # add layer length to layers_sample_replicates
#       layers_sample_replicates[[rep, layer+3]] <- layer_perc_exp
#     }
#   }
# }
# 
# 
# # save layers proportions for these samples
# write_csv(layers_sample_replicates, "tables/SIMPLI_analyses/clusters/distribution/layers_proportion_for_sample_replicates.csv")


##### 3. Calculate layers area ####
# # calculate layers' length
layers_sample_replicates$layer_1_l <- layers_sample_replicates$layer_1
layers_sample_replicates$layer_2_l <- layers_sample_replicates$layer_2 - layers_sample_replicates$layer_1
layers_sample_replicates$layer_3_l <- layers_sample_replicates$layer_3 - layers_sample_replicates$layer_2
layers_sample_replicates$layer_4_l <- layers_sample_replicates$layer_4 - layers_sample_replicates$layer_3
layers_sample_replicates$layer_5_l <- layers_sample_replicates$layer_5 - layers_sample_replicates$layer_4
layers_sample_replicates$layer_6_l <- layers_sample_replicates$layer_6 - layers_sample_replicates$layer_5

# calculate layers area
layers_sample_replicates$layer_1_area <- layers_sample_replicates$layer_1_l*layers_sample_replicates$max_width
layers_sample_replicates$layer_2_area <- layers_sample_replicates$layer_2_l*layers_sample_replicates$max_width
layers_sample_replicates$layer_3_area <- layers_sample_replicates$layer_3_l*layers_sample_replicates$max_width
layers_sample_replicates$layer_4_area <- layers_sample_replicates$layer_4_l*layers_sample_replicates$max_width
layers_sample_replicates$layer_5_area <- layers_sample_replicates$layer_5_l*layers_sample_replicates$max_width
layers_sample_replicates$layer_6_area <- layers_sample_replicates$layer_6_l*layers_sample_replicates$max_width
layers_sample_replicates$tot_ROI_area <- layers_sample_replicates$max_length*layers_sample_replicates$max_width

# calculate mm^2 area
layers_sample_replicates$tot_ROI_mm2 <- layers_sample_replicates$tot_ROI_area/1000000
layers_sample_replicates$layer_1_mm2 <- layers_sample_replicates$layer_1_area/1000000
layers_sample_replicates$layer_2_mm2 <- layers_sample_replicates$layer_2_area/1000000
layers_sample_replicates$layer_3_mm2 <- layers_sample_replicates$layer_3_area/1000000
layers_sample_replicates$layer_4_mm2 <- layers_sample_replicates$layer_4_area/1000000
layers_sample_replicates$layer_5_mm2 <- layers_sample_replicates$layer_5_area/1000000
layers_sample_replicates$layer_6_mm2 <- layers_sample_replicates$layer_6_area/1000000

# save layers proportions for these samples
write_csv(layers_sample_replicates, "tables/clusters/distribution/layers_proportion_for_sample_replicates.csv")

##### 4. Calculate number of cell_types (cluster) per layer #####

# import layers sizes
#layers_sample_replicates <- read_xlsx("tables/clusters/distribution/layers_proportion_for_sample_replicates.csv")

#import clustered_cells_new
clustered_cells_new <- read_csv("input/clustered_cells_new.csv")

#subset
clustered_cells_new_position <- clustered_cells_new[, colnames(clustered_cells_new) %in% c('Metadata_sample_name', res, 'X', 'Y')]
clustered_cells_new_position$cell_type <-clustered_cells_new_position$res_0.1_ids 
clustered_cells_new_position$res_0.1_ids <- NULL

colnames(clustered_cells_new_position)[colnames(clustered_cells_new_position) == 'Metadata_sample_name'] <- 'sample_replicates'
clustered_cells_new_position$sample_name <- stringr::str_extract(clustered_cells_new_position$sample_replicates, "[^_]*")

# sort samples names and assigned cells names
cell_types <- sort(unique(clustered_cells_new_position$cell_type))
replicates <- sort(unique(clustered_cells_new_position$sample_replicates))


# create a table to store data
cells_per_layer <- data.frame(matrix(ncol = 4))
colnames(cells_per_layer) <- c("sample_replicates", "layer", "cell_type", "cell_number")
cells_per_layer_all <- cells_per_layer

# calculate cell per layer
good = 0
wrong = 0
cells_missing = 0

for (rep in 1:length(replicates)){
  # filter for 1 sample replicate (ROI) at a time
  
  sample_rows <- clustered_cells_new_position[which(clustered_cells_new_position$sample_replicates == replicates[rep]),]
  
  for (cell in 1:length(cell_types)){
  
    # filter for 1 assigned cell at a time
    cells_rows <- sample_rows[which(sample_rows$cell_type == cell_types[cell]),]
    
    for (layer in 1:6){
    
      if (layer == 1) {
        # pick layer length
        layer_1_length <- as.numeric(layers_sample_replicates[grep(replicates[rep], 
                                                        layers_sample_replicates$sample_replicates), layer+3]) 
        
        # filter for layer length
        cells <- nrow(subset(cells_rows, Y <= layer_1_length))
        
        # add values in the temporary table
        cells_per_layer$sample_replicates <- replicates[rep]
        cells_per_layer$layer <- layer
        cells_per_layer$cell_type <- cell_types[cell]
        cells_per_layer$cell_number <- cells
        
        # append to final table
        cells_per_layer_all <- rbind(cells_per_layer_all, cells_per_layer)
        
        # change layer before
        layer_before <- layer_1_length
      }
      if (layer > 1) {
        # pick layer length
        layer_length <- as.numeric(layers_sample_replicates[grep(replicates[rep], 
                                                      layers_sample_replicates$sample_replicates), layer+3])  
        
        # filter for layer length
        cells <- nrow(subset(cells_rows, Y <= layer_length & Y > layer_before))
        
        # add values in the temporary table
        cells_per_layer$sample_replicates <- replicates[rep]
        cells_per_layer$layer <- layer
        cells_per_layer$cell_type <- cell_types[cell]
        cells_per_layer$cell_number <- cells
        
        # append to final table
        cells_per_layer_all <- rbind(cells_per_layer_all, cells_per_layer)
        
        # change layer before
        layer_before <- layer_length
      }
    }
  }
  
  temp1 <- filter(cells_per_layer_all, sample_replicates == replicates[rep])
  #temp2 <- filter(clusters_position, sample_replicates == replicates[rep])
  temp2 <- filter(clustered_cells_new_position, sample_replicates == replicates[rep])
  
  if (sum(temp1$cell_number) == nrow(temp2)){
    good = good + 1
  } 
  if (sum(temp1$cell_number) > nrow(temp2)){
    print(paste('replicate ', replicates[rep], ' has ', (sum(temp1$cell_number) - nrow(temp2)), 'more cells than expected'))
    wrong = wrong + 1
    cells_missing = cells_missing + abs(sum(temp1$cell_number) - nrow(temp2))
  } 
  if (sum(temp1$cell_number) < nrow(temp2)){
    print(paste('replicate ', replicates[rep], ' has ', (sum(temp1$cell_number) - nrow(temp2)) ,'less cells than expected'))
    wrong = wrong + 1
    cells_missing = cells_missing + abs(sum(temp1$cell_number) - nrow(temp2))
  }
}

# check if the first row is empty and remove it
if(is.na(cells_per_layer_all[1, 1])) {
  cells_per_layer_all <- cells_per_layer_all[-1, ]
}

# Save final table
write_csv(cells_per_layer_all, "tables/clusters/distribution/cluster_per_layer.csv")

rm(temp1, temp2, cell, cell_types, cells, cells_missing, good, wrong, 
   cells_rows, sample_rows, cells_per_layer, layer, layer_1_length, layer_before, layer_length, rep)

##### 5  Calculate cell density per layer #####

# load input data
#cells_per_layer_all <- read_xlsx("tables/clusters/distribution/cluster_per_layer.xlsx")
#layers_sample_replicates <- read_xlsx("tables/clusters/distribution/layers_proportion_for_sample_replicates.xlsx")

# adapt layers name
layers_sample_replicates <- layers_sample_replicates[, colnames(layers_sample_replicates) %in% 
                                                       c("sample_replicates", "layer_1_mm2", "layer_2_mm2", "layer_3_mm2", 
                                                         "layer_4_mm2", "layer_5_mm2", "layer_6_mm2")]
layers_sample_replicates <- layers_sample_replicates %>% 
  pivot_longer(cols = c("layer_1_mm2", "layer_2_mm2", "layer_3_mm2", 
                        "layer_4_mm2", "layer_5_mm2", "layer_6_mm2"), names_to = 'layer', values_to = 'mm2_area')
layers_sample_replicates$layer <- sapply(strsplit(layers_sample_replicates$layer, "_"), "[[", 2)

# add layers area
layers_sample_replicates$layer <- as.character(layers_sample_replicates$layer)
cells_per_layer_all$layer <- as.character(cells_per_layer_all$layer)

cells_per_layer_all <- left_join(cells_per_layer_all, layers_sample_replicates, by = c("sample_replicates", "layer"))

# calculate density
cells_per_layer_all$cell_density <- cells_per_layer_all$cell_number / cells_per_layer_all$mm2_area

# Save final table
write_xlsx(cells_per_layer_all, "tables/clusters/distribution/cluster_per_layer.xlsx")


##### 6. Calculate average number and density of cells per layer per sample and per group #####

# import file
#cells_per_layer_all <- read_xlsx("tables/clusters/distribution/cluster_per_layer.xlsx") 

# add column with sample name
cells_per_layer_all <- add_column(cells_per_layer_all, "sample_name" = 
                                    stringr::str_extract(cells_per_layer_all$sample_replicates, 
                                                         "[^_]*"), .before = "sample_replicates")

# calculate average cell number per layer, per SAMPLE, per cluster 
cells_per_layer_sample_avg <- aggregate(cells_per_layer_all$cell_number, 
                                        list(cells_per_layer_all$sample_name,
                                             cells_per_layer_all$layer, 
                                             cells_per_layer_all$cell_type), mean)

# calculate average density cell number per layer, per SAMPLE, per cluster 
density_cells_per_layer_sample_avg <- aggregate(cells_per_layer_all$cell_density, 
                                        list(cells_per_layer_all$sample_name,
                                             cells_per_layer_all$layer, 
                                             cells_per_layer_all$cell_type), mean)

# calculate sum cell number per layer, per SAMPLE, per cluster 
cells_per_layer_sample_sum <- aggregate(cells_per_layer_all$cell_number, 
                                        list(cells_per_layer_all$sample_name,
                                             cells_per_layer_all$layer, 
                                             cells_per_layer_all$cell_type), sum)
# change column names
colnames(cells_per_layer_sample_avg) <- c("sample_name", "layer", "cluster", "cell_number_avg")
colnames(density_cells_per_layer_sample_avg) <- c("sample_name", "layer", "cluster", "cell_density_avg")
colnames(cells_per_layer_sample_sum) <- c("sample_name", "layer", "cluster", "cell_number_sum")

# check in which cluster cells are lost 
sample_names <- unique(cells_per_layer_all$sample_name)

cell_types <- unique(cells_per_layer_all$cell_type)

missing_cells = 0
for (cell in 1:length(cell_types)){
  test <- cells_per_layer_sample_sum[which(cells_per_layer_sample_sum$cluster == as.integer(cell_types[cell])),]
  test2 <- clustered_cells_new_position[which(clustered_cells_new_position$cell_type == as.integer(cell_types[cell])),]
  
  if (nrow(test2) == sum(test$cell_number_tot) & # check if all cells were counted
      length(sample_names)*6 == nrow(test)) { # check if all layers are there
    print(paste0("all good for assigned cell ", cell_types[cell]))
  } else {
    print(paste0('We are missing ', abs(nrow(test2) - sum(test$cell_number_tot)), ' cells in cluster #', cell_types[cell]))
    missing_cells = missing_cells + abs(nrow(test2) - sum(test$cell_number_tot))
  }
}

#check total and percentage of cells are lost 
print(paste0('we are missing ', missing_cells, ' cells in total, out of ', sum(cells_per_layer_all$cell_number), 
             ' (', round(missing_cells/nrow(clustered_cells_new_position)*100, 3) ,'% of total cells)'))

# add column with group (if you have more than one mouse per group)
cells_per_layer_sample_avg <-
  add_column(cells_per_layer_sample_avg,
             "group" = stringr::str_extract(cells_per_layer_sample_avg$sample_name,
                                            "[^_]*"), .before = "sample_name")

# calculate average per group
cells_per_layer_group_avg <- aggregate(cells_per_layer_sample_avg$cell_number_avg,
                                       list(cells_per_layer_sample_avg$group,
                                            cells_per_layer_sample_avg$layer,
                                            cells_per_layer_sample_avg$cluster), mean)
# change column names
colnames(cells_per_layer_group_avg) <- c("sample_name", "layer", "cluster", "cell_number_avg")

# add column with group
cells_per_layer_sample_sum <-
  add_column(cells_per_layer_sample_sum,
             "group" = stringr::str_extract(cells_per_layer_sample_sum$sample_name,
                                            "[^_]*"), .before = "sample_name")

# calculate AVERAGE per group (from sum)
cells_per_layer_group_avg_sum <- aggregate(cells_per_layer_sample_sum$cell_number_tot,
                                       list(cells_per_layer_sample_sum$group,
                                            cells_per_layer_sample_sum$layer,
                                            cells_per_layer_sample_sum$cluster), mean)

# change column names
colnames(cells_per_layer_group_avg_sum) <- c("sample_name", "layer", "cluster", "cell_number_avg")


# calculate SUM per group
cells_per_layer_group_sum <- aggregate(cells_per_layer_sample_sum$cell_number_tot,
                                           list(cells_per_layer_sample_sum$group,
                                                cells_per_layer_sample_sum$layer,
                                                cells_per_layer_sample_sum$cluster), sum)

# change column names
colnames(cells_per_layer_group_sum) <- c("sample_name", "layer", "cluster", "cell_number_sum")

# calculate SUM per group per layer (not split in clusters) 
tot_cells_per_layer <- aggregate(cells_per_layer_group_sum$cell_number_sum,
                                 list(cells_per_layer_group_sum$sample_name,
                                      cells_per_layer_group_sum$layer), sum)
# change column names
colnames(tot_cells_per_layer) <- c("group", "layer", "cell_number_sum")

# save
write_xlsx(cells_per_layer_sample_avg, "tables/clusters/distribution/clusters_per_layer_sample_avg.xlsx")
write_xlsx(density_cells_per_layer_sample_avg, "tables/clusters/distribution/clusters_density_per_layer_sample_avg.xlsx")

#merge dfs
cells_per_layer_sample_avg_final <- left_join(cells_per_layer_sample_avg, 
                                              density_cells_per_layer_sample_avg, by = c('sample_name', 'layer', 'cluster'))
write_xlsx(cells_per_layer_sample_avg_final, "tables/clusters/distribution/clusters_per_layer_sample_avg_final.xlsx")

##### 7. Statistical analysis + plotting ######
##### 7.1 Define analysis to do --> OPEN <-- #####

analysis <- 'density' # here choose bw cell density or average cell number
#analysis <- 'cell_number'

##### 7.2 Define variables ######

if( analysis == 'density') {
  column_oi <- "cell_density_avg"
  folder <- "clusters"
  subfolder <- "distribution"
  subsubfolder <- "density"
  is_percentage <- FALSE
  plot_title <- "nuclei per sample / mm^2"
} else if (analysis == 'cell_number') {
  column_oi <- "cell_number_avg"
  folder <- "clusters"
  subfolder <- "distribution"
  subsubfolder <- "cell_number"
  is_percentage <- FALSE
  plot_title <- "nuclei per sample"
}

##### 7.3 Add variants ##### 

# read in list of sample/variants
variants <- read_csv("input/vars_per_sample.csv")

# add empty columns in dat and dat_long
cells_per_layer_sample_den$TREM2var <- 0
cells_per_layer_sample_den$ApoE4pos <- 0
cells_per_layer_sample_den$AD_stage <- 0
cells_per_layer_sample_den$TREM2var_ADstage <- 0
cells_per_layer_sample_den$ApoE4pos_TREM2var <- 0

# get all sample names
sample_names <- unique(cells_per_layer_sample_den$sample_name)

# assign variants to samples
for (sample in 1:length(sample_names)){
  # assign ApoE variant
  cells_per_layer_sample_den[cells_per_layer_sample_den$sample_name == sample_names[sample],]$ApoE4pos <-
    variants[variants$sample_name == sample_names[sample],]$ApoE4pos

  # assign TREM2 variant
  cells_per_layer_sample_den[cells_per_layer_sample_den$sample_name == sample_names[sample],]$TREM2var <-
    variants[variants$sample_name == sample_names[sample],]$TREM2var

  # assign AD_stage variant
  cells_per_layer_sample_den[cells_per_layer_sample_den$sample_name == sample_names[sample],]$AD_stage <-
    variants[variants$sample_name == sample_names[sample],]$AD_stage

  # assign TREM2var_ADstage variant
  cells_per_layer_sample_den[cells_per_layer_sample_den$sample_name == sample_names[sample],]$TREM2var_ADstage <-
    variants[variants$sample_name == sample_names[sample],]$TREM2var_ADstage

  # assign ApoE4pos_TREM2var variant
  cells_per_layer_sample_den[cells_per_layer_sample_den$sample_name == sample_names[sample],]$ApoE4pos_TREM2var <-
    variants[variants$sample_name == sample_names[sample],]$ApoE4pos_TREM2var
}


##### 7.4 For loop: 3 conditions, stat analysis and plot - CELLS PER CLUSTER PER LAYER #####

# define columns to keep
countedl <- cells_per_layer_sample_avg_final[,which(colnames(cells_per_layer_sample_avg_final) %in% c("group", "sample_name", "cluster", "layer" , column_oi))]

# TEMPORARLY MODIFICATION (only one singe mouse per group)
countedl$group <- countedl$sample_name
countedl$sample_name <- NULL

# make directories
dir.create(paste0("plots/", folder, "/", subfolder))
dir.create(paste0("tables/", folder, "/", subfolder))

dir.create(paste0("plots/", folder, "/", subfolder, "/", subsubfolder))
dir.create(paste0("tables/", folder, "/", subfolder, "/", subsubfolder))

folder_tables <- paste0("tables/", folder, "/", subfolder, "/", subsubfolder, "/") 
folder_plots <- paste0("plots/", folder, "/", subfolder, "/", subsubfolder, "/")

# Create table with conditions for for loop 
conditions <- data.frame(matrix(ncol = 3, nrow = 3))
conditions <- data.frame(
  cond_name = c(1, 2, 3),
  order = I(list(c("CtrlCV", "AlzCV"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzTREM2"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzR62H", "Alz472H"))),
  group_column = c("group", "group", "TREM2var")
)

##### Start for loop 

for (cond in 1:nrow(conditions)) {
  #cond = 1
  
  print(paste0(" ----- START CONDITION #", cond, " ------"))
  
  dat <- as.data.frame(countedl)
  x <- unlist(conditions$order[cond])
  layers_name <- unique(dat$layer)
  
  if ((conditions$group_column[cond]) != "group") {
    dat$group <- NULL
    colnames(dat)[which(colnames(dat) == conditions$group_column[cond])] <- 'group'
  }
  
  dat <- dat[dat$group %in% x,]

  # create. new directory for output
  dir.create(paste0(folder_plots, paste(x, collapse = "_")))
  dir.create(paste0(folder_tables, paste(x, collapse = "_")))
  
  # ##### Save file for paper #####
  # 
  # # select only columns of interest
  # #dat_paper <- dat[, which(colnames(dat) %in% c('group', "sample_name",'layer', 'cluster', column_oi))]
  # colnames(dat_paper)[which(colnames(dat) == column_oi)] <- "count"
  # 
  # dat_paper_wide <- dat_paper %>% pivot_wider(names_from = 'cluster',
  #                                             values_from = "count")
  # 
  # colnames(dat_paper_wide)[4:39] <- paste0("cluster_", colnames(dat_paper_wide)[4:39])
  # 
  # # save
  # write_xlsx(dat_paper_wide, paste0(folder_tables,  paste(x, collapse = "_"), "/fig_4.cell_density_per_layer_in_", paste(x, collapse = "_") , ".xlsx"))
  # 
  # # calculate average and SD for sample groups
  # 
  # #dat <- dat[,colnames(dat) %in% c("group", "cell_type", "prop_avg")]
  # 
  # dat_paper_avg <- dat_paper %>% 
  #   dplyr::group_by(group, cluster, layer) %>%
  #   dplyr::summarise(count = paste(round(mean(count), 2), 
  #                                  round(sd(count), 2), 
  #                                  sep = ' ± ')
  #   )
  # 
  # dat_paper_avg_wide <- dat_paper_avg %>% pivot_wider(names_from = 'cluster',
  #                                                     values_from = "count")
  # 
  # colnames(dat_paper_avg_wide)[3:ncol(dat_paper_avg_wide)] <- paste0("cluster_", colnames(dat_paper_avg_wide)[3:ncol(dat_paper_avg_wide)])
  # 
  # # save
  # write_xlsx(dat_paper_avg_wide, paste0(folder_tables, paste(x, collapse = "_"), "/fig_4.avg_cell_density_per_layer_in_", paste(x, collapse = "_") , ".xlsx"))
  # 
  ##### Preparative tasks #####
  
  # select only columns of interest
  dat <- dat[, which(colnames(dat) %in% c('group', 'layer', 'cluster', column_oi))]
  colnames(dat)[which(colnames(dat) == column_oi)] <- "count"
  
  # get groups and cell types analysed
  groups <- unique(dat$group)
  
  dat <- as.data.frame(dat)
  dat$count <- as.numeric(dat$count)
  
  ##### Transform percentage ########
  # explanation arcsin: https://www.programmingr.com/tutorial/arcsine-transformation/
  # explanation arcsin and logit with nice plots: http://strata.uga.edu/8370/rtips/proportions.html
  
  if (is_percentage) {
    # add column with arcsin, log and radarcsin
    dat$arcsin_perc <- asin(sqrt(dat$count/100))
    
    # rename arcsin_perc to count, so that statistical analysis will be done on that column
    # but keep count as percentage as it will be used for plotting
    colnames(dat)[which(colnames(dat) == 'count')] <- "percentage"
    colnames(dat)[which(colnames(dat) == 'arcsin_perc')] <- "count"
    
    # dat$log_perc <- log((((dat$count)+0.00001)/100)/
    #                       (1-((dat$count+0.00001)/100)))
    # 
    # dat$radarcsin_perc <- rad2deg(asin(sqrt(((dat$count)+0.00001)/100)))
    
    print("Cell number is a percentage and has been transformed")
  } else {
    print("Cell number is **NOT** percentage and has **NOT** been transformed")
  }
  
  ##### Test normality ########
  # testing normality with a Shapiro-Wilk Test.
  # source: https://www.statology.org/test-for-normality-in-r/
  
  #prepare data table for result of normality test
  tests_todo <- data.frame(matrix(ncol = 4))
  colnames(tests_todo) <- c('group', 'cluster', 'layer', 'dist')
  
  tests_todo_temp <- tests_todo
  
  # test normality 
  # group_by: group, layer and cell_type 
  for (group in 1:length(groups)){
    #group = 1
    skip_to_next <- FALSE
    tryCatch({
      #group = 1
      samples <- dat[dat$group == groups[group], ]
      
      for (cluster in unique(dat$cluster)){
        #cluster = unique(dat$cluster)[3]
        samples1 <- samples[samples$cluster == cluster, ]
        
        for (layer in 1:6){
          #layer = 1
          samples2 <- samples1[samples1$layer == layer, ]
          
          if ((shapiro.test(samples2$count)$p.value) < 0.05) {
            #print(paste0("Cluster #", cluster, " cells in layer ", layer, " in ", groups[group], " are **NOT** normally distributed"))
            
            tests_todo_temp$group <- groups[group]
            tests_todo_temp$cluster <- cluster
            tests_todo_temp$layer <- layer
            tests_todo_temp$dist <- 'kruskal.wilcox'
            
            tests_todo <- rbind(tests_todo, tests_todo_temp)
            
          } 
          
          if ((shapiro.test(samples2$count)$p.value) > 0.05) {
            #print(paste0("Cluster #", cluster, " cells in layer ", layer, " in ", groups[group], " are normally distributed"))
            
            tests_todo_temp$group <- groups[group]
            tests_todo_temp$cluster <- cluster
            tests_todo_temp$layer <- layer
            tests_todo_temp$dist <- 'anova.tukey'
            
            tests_todo <- rbind(tests_todo, tests_todo_temp)
          }
        } 
        
      }
      
    }, error = function(e) {skip_to_next <<- TRUE})
    if(skip_to_next) { next }     
  }
  
  # remove first row
  if (is.na(tests_todo[1,1])){
    tests_todo <- tests_todo[-1,] 
  }
  tests_todo <- as.data.frame(tests_todo)
  
  #prepare data table for stat tests to use per cell_type
  tests_todo_final <- data.frame(matrix(ncol = 3))
  colnames(tests_todo_final) <- c('cluster', 'layer', 'test')
  
  tests_todo_final_temp <- tests_todo_final
  
  # check which layers are normally distributed and define stat test to use for analysis
  for (layer in 1:6){
    #layer = 1 
    samples <- tests_todo[tests_todo$layer == layer, ]
    
    for (cluster in unique(dat$cluster)){
      
      #cluster <- unique(dat$cluster)[1]
      samples1 <- samples[samples$cluster == cluster, ]
      
      if (nrow(samples) > 1) {
        if (all(samples$dist == "anova.tukey")){
          tests_todo_final_temp$cluster <- cluster
          tests_todo_final_temp$layer <- layer
          tests_todo_final_temp$test <- "anova.tukey"
        } else {
          tests_todo_final_temp$cluster <- cluster
          tests_todo_final_temp$layer <- layer
          tests_todo_final_temp$test <- "kruskal.wilcox"
        }
      } else {
        tests_todo_final_temp$cluster <- cluster
        tests_todo_final_temp$layer <- layer
        tests_todo_final_temp$test <- "kruskal.wilcox"
      }
      tests_todo_final <- rbind(tests_todo_final, tests_todo_final_temp)
      
    }
    
  }
  
  # remove first row
  if (is.na(tests_todo_final[1,1])){
    tests_todo_final <- tests_todo_final[-1,] 
  }
  tests_todo_final <- as.data.frame(tests_todo_final)
  rm(tests_todo_final_temp)
  rm(tests_todo_temp)
  
  # only some data are normally distributed:
  # --> ANOVA + Tukey on normally distributed data
  # --> Kruskal-wallis + Wilcoxon test on NOT normally distributed data
  
  ##### Wilcox test (when only two groups are present) #####
  
  # generate stat.test file to add tukey p.adj of count, for plotting
  stat.test <- data.frame(matrix(ncol = 8))
  colnames(stat.test) <- c('group1', 'group2', 'cluster', 'layer', 'p.adj', 'method', 'anova_krus', 'y.position')
  
  # independent 2-group Mann-Whitney U Test
  # wilcox.test(y~A)
  # where y is numeric and A is A binary factor
  
  if (length(groups) == 2){
    
    stat.test_temp <- stat.test

    for (layer in 1:length(layers_name)) {
      #layer = 1
      samples <- dat[dat$layer == layers_name[layer], ]
      
      for (cluster in unique(dat$cluster)){
        #cluster = unique(dat$cluster)[1]
        samples1 <- samples[samples$cluster == cluster, ]
        
        wt_res <-  wilcox.test(samples1$count ~ samples1$group)
        
        stat.test_temp$group1 <- groups[1]
        stat.test_temp$group2 <- groups[2]
        stat.test_temp$cluster <- cluster
        stat.test_temp$layer <- layers_name[layer]
        stat.test_temp$p.adj <- wt_res$p.value
        stat.test_temp$method <- "wilcox.test"
        stat.test_temp$anova_krus <- "N/A"
        stat.test_temp$y.position <- max(dat$count)*0.2
        
        stat.test <- rbind(stat.test, stat.test_temp)
        
      }
      
    }
    print("Wilcox test was performed as only two groups are being tested")
    rm(stat.test_temp)
  } else {
    print("Wilcox test not applied as groups analysed are more than 2")
  }
  
  # remove first row
  if (is.na(stat.test[1,1])){
    stat.test <- stat.test[-1,] 
  }
  stat.test <- as.data.frame(stat.test)
  
  ##### ANOVA + Tukey (for normally distributed data) #####
  # source: https://www.scribbr.com/statistics/anova-in-r/#:~:text=We%20can%20perform%20an%20ANOVA,levels%20of%20the%20independent%20variable.
  
  # remember to change layer names as you did for Wilcox
  
  if (length(groups) > 2){
    
    # filter samples for ANOVA
    temp <- filter(tests_todo_final, test == 'anova.tukey')
    
    layers_anova <- temp$layer
    
    samples_anova <- dat %>% filter(layer %in% layers_anova)
    
    # check we have all cells
    sample_names <- unique(dat$sample_name)
    nrow(samples_anova) == length(layers_anova)*length(sample_names)*unique(dat$cluster)
    
    # filter samples per layer
    if (nrow(samples_anova) > 0){
      
      for (cluster in unique(dat$cluster)){
        #cluster = unique(dat$cluster)[1]
        samples_anova1 <- samples_anova[samples_anova$cluster == cluster, ]
        
        for (layer in layers_anova){
          #layer = 2
          samples_anova2 <- samples_anova[samples_anova$layer == layer, ]
          
          #### ANOVA + Tukey test
          # save ANOVA results
          one.way <- aov(count ~ group, data = samples_anova2)
          
          # perform tukey
          tukey.one.way <- TukeyHSD(one.way)
          
          # export into dataframe
          tukey.one.way.df <- as.data.frame(tukey.one.way[1])
          tukey.one.way.df <- rownames_to_column(tukey.one.way.df)
          
          # extract group1 and group2 of comparison
          tukey.one.way.df$group1 <- sub(".*\\-", "", tukey.one.way.df$rowname)
          tukey.one.way.df$group2 <- sub("\\-.*", "", tukey.one.way.df$rowname)
          
          # modify tukey.one.way.df to fit stat.test
          colnames(tukey.one.way.df)[colnames(tukey.one.way.df) == 'group.p.adj'] <- 'p.adj'
          tukey.one.way.df$rowname <- NULL
          tukey.one.way.df$group.diff <- NULL
          tukey.one.way.df$group.lwr <- NULL
          tukey.one.way.df$group.upr <- NULL
          
          tukey.one.way.df <- tukey.one.way.df %>% relocate(p.adj, .after = group2)
          tukey.one.way.df$method <- 'anova.tukey'
          tukey.one.way.df$layer <- layer
          tukey.one.way.df <- tukey.one.way.df %>% relocate(layer, .after = group2)
          tukey.one.way.df$cluster <- cluster
          tukey.one.way.df <- tukey.one.way.df %>% relocate(cluster, .after = group2)
          
          # add y.position and Pr(>F)
          y_increase <- max(samples_anova$count)*0.2
          for (n in 1:nrow(tukey.one.way.df)){
            tukey.one.way.df$y.position[n] <- max(samples_anova$count)+(y_increase+((y_increase/2)*(n-1)))
            tukey.one.way.df$anova_krus[n] <- summary(one.way)[[1]][["Pr(>F)"]][1]
          }
          stat.test <- rbind(stat.test, tukey.one.way.df)
        }
        
      }
      
      print("ANOVA + Tukey tests performed on multiple groups")
      
    } else {
      print("No samples for ANOVA")
    }
  } else {
    print("ANOVA + Tukey tests not performed as only two groups are tested")
  }
  
  ##### Kruskal + Wilcoxon (for not normally distributed data) #####
  
  
  # remember to change layer names as you did for Wilcox
  
  if (length(groups) > 2){
    
    # filter samples for Kruskal
    temp <- filter(tests_todo_final, test == 'kruskal.wilcox')
    
    layers_kruskal <- unique(temp$layer)
    
    samples_kruskal <- dat %>% filter(layer %in% layers_kruskal)
    
    # check we have all cells
    sample_names <- unique(dat$sample_name)
    nrow(samples_kruskal) == length(layers_kruskal)*length(sample_names)
    
    if (nrow(samples_kruskal) > 0){
      
      for (cluster in unique(dat$cluster)){
        #cluster = unique(dat$cluster)[1]
        samples_kruskal1 <- samples_kruskal[samples_kruskal$cluster == cluster, ]
        
        for (layer in layers_kruskal){
          #layer = layers_kruskal[1]
          samples_kruskal2 <- samples_kruskal1[samples_kruskal1$layer == layer, ]
          
          # perform Kruskal test and save
          kruskal.test <- kruskal.test(count ~ group, data = samples_kruskal2)
          
          # perform Wilcox test and save
          wilcox.test <- pairwise.wilcox.test(samples_kruskal2$count, samples_kruskal2$group,
                                              p.adjust.method = "BH", 
                                              exact = FALSE)
          
          # modify wilcox.test to fit stat.test
          wilcox.test <- rownames_to_column(as.data.frame(wilcox.test$p.value))
          wilcox.test <- melt(wilcox.test)
          colnames(wilcox.test) <- c('group2', 'group1', 'p.adj')
          wilcox.test$method <- 'kruskal.wilcox'
          
          wilcox.test <- wilcox.test %>% relocate(group2, .after = group1)
          wilcox.test$layer <- layer
          wilcox.test <- wilcox.test %>% relocate(layer, .after = group2)
          wilcox.test$cluster <- cluster
          wilcox.test <- wilcox.test %>% relocate(cluster, .after = group2)
          
          # add y.position and Pr(>F) and remove rows from same samples
          y_increase <- max(samples_kruskal1$count)*0.2
          for (n in 1:nrow(wilcox.test)){
            wilcox.test$y.position[n] <- max(samples_kruskal1$count)+(y_increase+((y_increase/2)*(n-1)))
            wilcox.test$anova_krus[n] <- kruskal.test$p.value
            if (is.nan(wilcox.test$p.adj[n])){
              wilcox.test$p.adj[n] <- 0.999
            } 
          }  
          # remove rows with NA
          wilcox.test <- na.omit(wilcox.test)
          
          # bind to previous tibble
          stat.test <- rbind(stat.test, wilcox.test)
          
        } 
        
      }
      
      print("Kruskal + Wilcoxon tests performed on multiple groups")
      
    } else {
      print("No samples for Kruskal")
    }
    
  } else {
    print("Kruskal + Wilcoxon tests not performed as only two groups are tested")
  }
  
  # remove first row
  if (is.na(stat.test[1,1])){
    stat.test <- stat.test[-1,] 
  }
  stat.test <- as.data.frame(stat.test)
  rm(stat.test_temp)
  
  ##### Add label and save ########
  
  # add label
  stat.test$label <- ""
  stat.test <- mutate(stat.test, label = dplyr::case_when(
    stat.test$p.adj <= 0.001 ~ "***",
    stat.test$p.adj <= 0.01 ~ "**",
    stat.test$p.adj <= 0.05 ~ "*",
    stat.test$p.adj > 0.05 ~ "",
    is.na(stat.test$p.adj) ~ "",
  )
  )

  ##### PLOT with p.adj generated elsewhere #####
  # source: https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
  
  # load clusters labels
  cluster_labels <- read.xlsx("input/clusters_labels.xlsx")
  cluster_labels <- cluster_labels[, c('cluster', 'label', 'full_label', 'cell_type')]
  cluster_labels <- na.omit(cluster_labels)
  
  # select only clusters of interest
  cluster_labels <- cluster_labels[gsub("cluster_", "", cluster_labels$cluster) %in% unique(dat$cluster),]
  
  dat_plot <- dat
  
  # re-order groups order
  dat_plot <- dat_plot %>% arrange(sapply(group, function(y) which(y == x)))

  
  # filter p-adj that are not significant !!!1 DO NOT RUN UNTIL YPU HAVE SIGNIFICANT CLUSTER !!!
  #stat.test.signif <- filter(stat.test, p.adj <= 0.05)
  stat.test.signif <- stat.test # use this in
  
  # ylim
  ylim_plot <- max(stat.test.signif$y.position)
  
  # set as dataframe 
  dat_plot <- as.data.frame(dat_plot)
  
  layers_label <- as_labeller(
    c("1" = "L1", "2" = "L2", "3" = "L3", "4" = "L4", "5" = "L5", "6" = "L6")
  )
  
  groups_label <- c(
    "AlzCV" = "#D02D27",
    "AlzTREM2" = "#D39898",
    "AlzR47H" = "#D39898",
    "AlzR62H" = "#D8C5C5",
    "CtrlCV" = "#314B93",
    "CtrlTREM2" = "#8297C4"
  )
  
  # # change y.positions
  # stat_temp2 <- data.frame(matrix(ncol = ncol(stat.test.signif)))
  # colnames(stat_temp2) <- colnames(stat.test.signif)
  # 
  # if (nrow(stat.test.signif) > 0) {
  #   
  #   for (i in 1:length(unique(stat.test.signif$layer))){
  #     samples_temp <- filter(dat_plot, dat_plot$layer == unique(stat.test.signif$layer)[i])
  #     stat_temp <- filter(stat.test.signif, stat.test.signif$layer == unique(stat.test.signif$layer)[i])
  #     
  #     for (a in 1:nrow(stat_temp)){
  #       stat_temp$y.position[a] <- max(samples_temp$count) + ((ylim_plot/30)*a)
  #     }
  #     stat_temp2 <- rbind(stat_temp2, stat_temp)
  #   }
  #   
  #   # remove first row
  #   if (is.na(stat_temp2[1,1])){
  #     stat_temp2 <- stat_temp2[-1,]
  #   }
  #   
  #   stat.test.signif <- stat_temp2
  #   rm(stat_temp2, stat_temp, samples_temp)
  # }
  
  # generate one BOX plot per cluster
  clusters <- unique(dat$cluster)
  
  for (cluster in clusters){
    skip_to_next <- FALSE
    tryCatch({
      #cluster=8
      # filter per cluster
      samples <- dat_plot[dat_plot$cluster == cluster, ]
      stat <- stat.test.signif[stat.test.signif$cluster == cluster, ]
      stat <- droplevels(stat)
      
      stat_temp2 <- data.frame(matrix(ncol = ncol(stat)))
      colnames(stat_temp2) <- colnames(stat)
      
      if (nrow(stat) > 0) {
        # change y.positions
        for (i in 1:length(unique(stat$layer))){
          samples_temp <- filter(samples, samples$layer == unique(stat$layer)[i])
          stat_temp <- filter(stat, stat$layer == unique(stat$layer)[i])
          
          for (a in 1:nrow(stat_temp)){
            stat_temp$y.position[a] <- max(samples_temp$count) + ((ylim_plot/50)*a)
          }
          stat_temp2 <- rbind(stat_temp2, stat_temp)
        }
        
        # remove first row
        if (is.na(stat_temp2[1,1])){
          stat_temp2 <- stat_temp2[-1,]
        }
        
        stat <- stat_temp2
        rm(stat_temp2, stat_temp, samples_temp)
      }
      
      # Box plot facetted by "marker"
      p <- ggboxplot(samples, x = "group", y = "count",
                     color = "group", 
                     #palette = c("#314B93","#8297C4","#D02D27","#D39898"),
                     #palette = c("#314B93", "#D02D27"),
                     #palette = c("#314B93", "#8297C4","#D02D27", "#D39898", "#D8C5C5"),
                     palette = groups_label,
                     add = "jitter", 
                     ylab = cluster_labels$full_label[grep(paste0("cluster_", cluster), cluster_labels$cluster)],
                     short.panel.labs = TRUE, 
                     legend = NULL) +
        facet_wrap(~ layer, 
                   ncol = 6, 
                   switch = "x", 
                   labeller = layers_label) +
        theme(legend.title = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), 
              strip.text.x = element_text(angle = 90), 
              axis.title.x = element_blank()) +
        #ylim(0,ylim_plot) +
        rotate_y_text(angle = 90) +
        rremove("legend") +
        border()
      
      # Use only p.format as label. Remove method name.
      p <- p + stat_pvalue_manual(stat,
                                  label = "label",
                                  tip.length = 0,
                                  label.size = 7
      )
      
      pdf(paste0(folder_plots, 
                 paste(x, collapse = "_"), 
                 "/density_cluster_", 
                 cluster, 
                 "_distribution_in_", 
                 paste(x, collapse = "_"), 
                 ".pdf"), 
          width = 7, 
          height = 3)
      print(p)
      dev.off()
      
      write_csv(stat, paste0(folder_tables, 
                             paste(x, collapse = "_"), 
                             "/stats_average_cluster", 
                             cluster, 
                             "_distirbution_in_", 
                             paste(x, collapse = "_"), 
                             ".csv"))
      
      # print plot done
      print(paste0("Plot for cluster_", cluster, " - condition #", cond, " was generated"))
      
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }     
  }
  
  print(paste0(" ----- END CONDITION #", cond, " ------"))
  
}
