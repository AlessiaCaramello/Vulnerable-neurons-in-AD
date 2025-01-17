# Calculate cells positive for markers of interest

# input files:
# - clustered_cells.csv (SIMPLI output)
# - markers_threshold.csv (metadata)
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
#install.packages("REdaS")
library(REdaS)
library(RColorBrewer)
#install.packages("viridis")
library(viridis)
#install.packages("writexl")
library(writexl)
library(readxl)
library(openxlsx)

# set working directory
setwd("~/UK Dementia Research Institute Dropbox/Alessia Caramello/Alessia Lab Dropbox/Projects/UKDRI project/Image analysis/HistoCat:SIMPLI analysis/051022 TREM2 cohort/R code")

# create output folder
#dir.create("plots/clusters/cells_positive_for_markers")
#dir.create("tables/clusters/cells_positive_for_markers")

##### 1. Import files, define cluster resolution and markers of interest #####

# import files
clustered_cells_new <- read.csv("input/clustered_cells.csv")

# define resolution to analyse
res <- 'res_0.9_ids'

# define markers of interest
markers_oi <- c("Ab", "pTau")

# create folder where to store 
dir.create(paste0("plots/clusters/cells_positive_for_markers/", paste0(markers_oi, collapse = '_'), '/'))
dir.create(paste0("tables/clusters/cells_positive_for_markers/", paste0(markers_oi, collapse = '_'), '/'))

# define paths for savings tables and plots
path_tables <- paste0("tables/clusters/cells_positive_for_markers/", paste0(markers_oi, collapse = '_'), '/')
path_plots <- paste0("plots/clusters/cells_positive_for_markers/", paste0(markers_oi, collapse = '_'), '/')

##### 2. Keep only columns with intensity of markers of interest #####

# keep only useful columns
clusters_intensity <- clustered_cells_new %>% 
  select(c("Metadata_sample_name", 
           as.character(res), 
           'X', 'Y',
           paste0("Intensity_MeanIntensity_", markers_oi)
           ))

colnames(clusters_intensity)[5:ncol(clusters_intensity)] <- gsub("Intensity_MeanIntensity_", "", colnames(clusters_intensity)[5:ncol(clusters_intensity)])

# add rows with sample_name and group
colnames(clusters_intensity)[1:2] <- c('sample_replicate', 'cluster')

# extract sample name and group
clusters_intensity$sample_name = str_extract(clusters_intensity$sample_replicate, "[^_]*_[^_]*")
clusters_intensity$group = sub("_.*", "", clusters_intensity$sample_replicate)

# exception: run when you have only 1 sample per group
clusters_intensity$sample_name <- paste0(clusters_intensity$group, '_', sapply(strsplit(clusters_intensity$sample_replicate, "_"), "[[", 3))

##### 3. Count number of cells with intensity above custom threshold #####

# generate dataframe to store markers and thresholds values
# markers_threshold <- as.data.frame(matrix(nrow = length(markers_oi), ncol = 2))
# colnames(markers_threshold) <- c('markers', 'threshold')
# markers_threshold$markers <- markers_oi
#write_xlsx(markers_threshold, 'input/markers_threshold.xlsx')

# fill in thresholds and load file again
markers_threshold <- read_csv('input/markers_threshold.csv')

############### for loop for counting positive cells for EACH marker of interest
for (i in 1:nrow(markers_threshold)) {
  #i=1
  marker <- as.character(markers_threshold$markers[i])
  threshold <- as.numeric(markers_threshold$threshold[i])
  
  count_cells <- clusters_intensity %>% 
    group_by(cluster, sample_name) %>% 
    dplyr::reframe(marker_pos = count(.data[[marker]] > threshold)) %>% 
    as.data.frame() %>%
    ungroup()
  
  count_cells <- count_cells[!(count_cells$marker_pos$x == "FALSE"), ] 
  count_cells$marker_pos2 <- count_cells$marker_pos$freq
  count_cells$marker_pos <- NULL
  
  # add missing samples
  samples <- unique(clusters_intensity$sample_name)
  clusters <- unique(clusters_intensity$cluster)
  
  for (c in 1:length(clusters)){
    temp1 <- count_cells %>% subset(count_cells$cluster == clusters[c])
    if (length(setdiff(samples, unique(temp1$sample_name))) > 0){
      diff <- setdiff(unique(clusters_intensity$sample_name), unique(temp1$sample_name))
      #print(paste0(paste(diff, collapse = ", "), " samples are missing in cluster_", clusters[c]))
      for (i in 1:length(diff)){
        temp2 <- as.data.frame(matrix(ncol = ncol(count_cells)))
        colnames(temp2) <- colnames(count_cells)
        temp2$sample_name <- diff[i]
        temp2$marker_pos2 <- as.numeric(0)
        temp2$cluster <- clusters[c]
        count_cells <- rbind(count_cells, temp2)
      }
    }
  }
  
  # add missing clusters
  for (s in 1:length(samples)){
    temp1 <- count_cells %>% subset(count_cells$sample_name == samples[s])
    if (length(setdiff(clusters, unique(temp1$cluster))) > 0){
      diff <- setdiff(clusters, unique(temp1$cluster))
      #print(paste0(paste(diff, collapse = ", "), " clusters are missing from sample ", samples[s]))
      for (i in 1:length(diff)){
        temp2 <- as.data.frame(matrix(ncol = ncol(count_cells)))
        colnames(temp2) <- colnames(count_cells)
        temp2$sample_name <- samples[s] 
        temp2$marker_pos2 <- as.numeric(0)
        temp2$cluster <- diff[i]
        count_cells <- rbind(count_cells, temp2)
      }
    }
  }
  
  if (length(clusters)*length(samples) == nrow(count_cells)) {
    print(paste0('All good for marker ', marker))
  } else {
    print(paste0('Something is wrong with marker ', marker))
  }
  
  # save dataframe
  write.xlsx(count_cells, paste0(path_tables, 'clusters_positive_for_', marker, '.xlsx'))
}

rm(count_cells, temp1, temp2)

##### 4. Count total number of cells per cluster and calculate % of marker+ cells #####

## read into one single dataframe all of the dataframes generated above
file_list <- list.files(path = path_tables, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(file_list, read_csv)

# change names
names(data_list) <- basename(sub(".csv", "", sapply(strsplit(sub(".*/", "", file_list), "_"), "[[", 4)))

# combine data
combined_data <- bind_rows(data_list, .id = "marker")
combined_data <- pivot_wider(combined_data, names_from = marker, values_from = marker_pos2)

# calculate total cells per cluster
tot_cells <- clusters_intensity %>% 
  group_by(cluster, sample_name) %>% 
  dplyr::reframe(tot = n()) %>% 
  as.data.frame() %>%
  ungroup()

# add tot_cells to combined_data
combined_data <- left_join(combined_data, tot_cells, by = c("cluster" = "cluster", "sample_name" = "sample_name"))

# substitute NA values
combined_data$tot[which(is.na(combined_data$tot))] <- as.numeric(0)

# calculate proportion
combined_data <- combined_data %>%
  dplyr::mutate(across(all_of(markers_oi), 
                ~ .x / .data[['tot']] * 100, 
                .names = "perc_{.col}"))

# substitute NA values
combined_data[is.na(combined_data)] <- as.numeric(0)

# add group
combined_data$group <- purrr::map_chr(combined_data$sample_name, ~strsplit(., "_")[[1]][1])

# check we have all cells 
if(sum(combined_data$tot) == nrow(clustered_cells_new)){
  print("all cells are included in combined_data")
} else {
  print("cells MISSING in combined_data")
}

# save
write.xlsx(combined_data, "tables/clusters/cells_positive_for_markers/perc_clusters_positive_for_markers.xlsx")

# transform combined_data to long format
combined_data_long <- combined_data[,c(1,2,
                                       max(grep('perc', colnames(combined_data)))+1,
                                       grep('perc', colnames(combined_data))
                                       )
                                    ]


combined_data_long <- pivot_longer(
  combined_data_long,
  cols = names(combined_data_long)[grep('perc', colnames(combined_data_long))],      
  names_to = "marker",  
  values_to = "perc_cells"  
)

combined_data_long$marker <- gsub('perc_', '', combined_data_long$marker)

# save
write.xlsx(combined_data_long, "tables/clusters/cells_positive_for_markers/perc_clusters_positive_for_markers_long.xlsx")

##### 5. Statistical analysis + plotting ######
##### 5.1 Define analysis to do  --> OPEN <-- #####

#analysis <- 'density' # here choose bw cell density or average cell number
#analysis <- 'cell_number'
analysis <- 'percentage'

##### 5.2 Define variables ######

if(analysis == 'density') {
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
} else if (analysis == 'percentage') {
  column_oi <- "perc_cells"
  folder <- "clusters"
  subfolder <- "cells_positive_for_markers"
  subsubfolder <- "cell_number"
  is_percentage <- TRUE
  plot_title <- "% marker+ nuclei per sample"
}


##### 5.3 Add variants ##### 

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


##### 5.4 For loop: stat analysis and plot (boxplot and heatmap) #####

# define columns to keep
countedl <- combined_data_long

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
  cond = 1
  
  print(paste0(" ----- START CONDITION #", cond, " ------"))
  
  dat <- as.data.frame(countedl)
  x <- unlist(conditions$order[cond])
  
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
  dat <- dat[, which(colnames(dat) %in% c('group', 'marker', 'cluster', column_oi))]
  colnames(dat)[which(colnames(dat) == column_oi)] <- "count"
  
  # get groups and markers analysed
  groups <- unique(dat$group)
  clusters <- unique(dat$cluster)
  markers <- unique(dat$marker)
  
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
  colnames(tests_todo) <- c('group', 'cluster', 'marker', 'dist')
  
  tests_todo_temp <- tests_todo
  
  # test normality 
  # group_by: group, layer and cell_type 
  for (group in 1:length(groups)){
    skip_to_next <- FALSE
    tryCatch({
      #group = 1
      samples <- dat[dat$group == groups[group], ]
      
      for (cluster in clusters){
        #cluster = clusters[1]
        samples1 <- samples[samples$cluster == cluster, ]
        
        for (marker in markers){
          #marker <- markers[1]
          samples2 <- samples1[samples1$marker == marker, ]
          
          if ((shapiro.test(samples2$count)$p.value) < 0.05) {
            #print(paste0("Cluster #", cluster, " cells in layer ", layer, " in ", groups[group], " are **NOT** normally distributed"))
            
            tests_todo_temp$group <- groups[group]
            tests_todo_temp$marker <- marker
            tests_todo_temp$layer <- layer
            tests_todo_temp$dist <- 'kruskal.wilcox'
            
            tests_todo <- rbind(tests_todo, tests_todo_temp)
            
          } 
          
          if ((shapiro.test(samples2$count)$p.value) > 0.05) {
            #print(paste0("Cluster #", cluster, " cells in layer ", layer, " in ", groups[group], " are normally distributed"))
            
            tests_todo_temp$group <- groups[group]
            tests_todo_temp$marker <- marker
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
  colnames(tests_todo_final) <- c('cluster', 'marker', 'test')
  
  tests_todo_final_temp <- tests_todo_final
  
  # check which markers are normally distributed and define stat test to use for analysis
  for (marker in markers){
    #marker <- markers[1] 
    samples <- tests_todo[tests_todo$marker == marker, ]
    
    for (cluster in unique(dat$cluster)){
      
      #cluster <- unique(dat$cluster)[1]
      samples1 <- samples[samples$cluster == cluster, ]
      
      if (nrow(samples) > 1) {
        if (all(samples$dist == "anova.tukey")){
          tests_todo_final_temp$cluster <- cluster
          tests_todo_final_temp$marker <- marker
          tests_todo_final_temp$test <- "anova.tukey"
        } else {
          tests_todo_final_temp$cluster <- cluster
          tests_todo_final_temp$marker <- marker
          tests_todo_final_temp$test <- "kruskal.wilcox"
        }
      } else {
        tests_todo_final_temp$cluster <- cluster
        tests_todo_final_temp$marker <- marker
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
  colnames(stat.test) <- c('group1', 'group2', 'cluster', 'marker', 'p.adj', 'method', 'anova_krus', 'y.position')
  
  # independent 2-group Mann-Whitney U Test
  # wilcox.test(y~A)
  # where y is numeric and A is A binary factor
  
  if (length(groups) == 2){
    
    stat.test_temp <- stat.test
    
    for (marker in markers) {
      #marker <- markers[1]
      samples <- dat[dat$marker == marker, ]
      
      for (cluster in unique(dat$cluster)){
        #cluster = unique(dat$cluster)[1]
        samples1 <- samples[samples$cluster == cluster, ]
        
        wt_res <-  wilcox.test(samples1$count ~ samples1$group)
        
        stat.test_temp$group1 <- groups[1]
        stat.test_temp$group2 <- groups[2]
        stat.test_temp$cluster <- cluster
        stat.test_temp$marker <- marker
        stat.test_temp$p.adj <- wt_res$p.value
        stat.test_temp$method <- "wilcox.test"
        stat.test_temp$anova_krus <- "N/A"
        stat.test_temp$y.position <- max(dat$percentage)*1.1
        
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
    
    markers_anova <- temp$marker
    
    samples_anova <- dat %>% filter(marker %in% markers_anova)
    
    # check we have all cells
    sample_names <- unique(dat$sample_name)
    nrow(samples_anova) == length(markers_anova)*length(sample_names)*unique(dat$cluster)
    
    # filter samples per layer
    if (nrow(samples_anova) > 0){
      
      for (cluster in unique(dat$cluster)){
        #cluster = unique(dat$cluster)[1]
        samples_anova1 <- samples_anova[samples_anova$cluster == cluster, ]
        
        for (marker in markers_anova){
          #layer = 2
          samples_anova2 <- samples_anova[samples_anova$marker == marker, ]
          
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
          tukey.one.way.df$marker <- marker
          tukey.one.way.df <- tukey.one.way.df %>% relocate(marker, .after = group2)
          tukey.one.way.df$cluster <- cluster
          tukey.one.way.df <- tukey.one.way.df %>% relocate(cluster, .after = group2)
          
          # add y.position and Pr(>F)
          y_increase <- max(samples_anova$percentage)*1.1
          for (n in 1:nrow(tukey.one.way.df)){
            tukey.one.way.df$y.position[n] <- max(samples_anova$percentage)+(y_increase+((y_increase/2)*(n-1)))
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
    
    markers_kruskal <- unique(temp$marker)
    
    samples_kruskal <- dat %>% filter(marker %in% markers_kruskal)
    
    # check we have all cells
    sample_names <- unique(dat$sample_name)
    nrow(samples_kruskal) == length(markers_kruskal)*length(sample_names)
    
    if (nrow(samples_kruskal) > 0){
      
      for (cluster in unique(dat$cluster)){
        #cluster = unique(dat$cluster)[1]
        samples_kruskal1 <- samples_kruskal[samples_kruskal$cluster == cluster, ]
        
        for (marker in markers_kruskal){
          #marker = markers_kruskal[1]
          samples_kruskal2 <- samples_kruskal1[samples_kruskal1$marker == marker, ]
          
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
          wilcox.test$marker <- marker
          wilcox.test <- wilcox.test %>% relocate(marker, .after = group2)
          wilcox.test$cluster <- cluster
          wilcox.test <- wilcox.test %>% relocate(cluster, .after = group2)
          
          # add y.position and Pr(>F) and remove rows from same samples
          y_increase <- max(samples_kruskal1$percentage)*1.1
          for (n in 1:nrow(wilcox.test)){
            wilcox.test$y.position[n] <- max(samples_kruskal1$percentage)+(y_increase+((y_increase/2)*(n-1)))
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

  ##### Plot as BOXPLOT #####
  # source: https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
  
  # load clusters labels
  cluster_labels <- read.xlsx("input/clusters_labels.xlsx")
  cluster_labels <- cluster_labels[, c('cluster', 'label', 'full_label', 'cell_type')]
  cluster_labels <- na.omit(cluster_labels)
  
  # select only clusters of interest
  cluster_labels <- cluster_labels[gsub("cluster_", "", cluster_labels$cluster) %in% unique(dat$cluster),]
  
  # change cluster names
  dat_plot <- dat
  dat_plot$cluster <- cluster_labels$full_label[match(paste0("cluster_", dat_plot$cluster), cluster_labels$cluster)]

  dat_plot <- dat_plot[!dat_plot$cluster %in% c('positive_for_all', 'negative_for_all'), ]
  
  # re-order groups order
  dat_plot <- dat_plot %>% arrange(sapply(group, function(y) which(y == x)))
  
  
  # filter p-adj that are not significant !!!1 DO NOT RUN UNTIL YPU HAVE SIGNIFICANT CLUSTER !!!
  #stat.test.signif <- filter(stat.test, p.adj <= 0.05)
  stat.test.signif <- stat.test # use this when there's nothing significant
  
  # ylim
  ylim_plot <- max(stat.test.signif$y.position)
  
  # set as dataframe 
  dat_plot <- as.data.frame(dat_plot)
  
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
  
  # Box plot facetted by "marker"
  p <- ggboxplot(dat_plot, x = "group", y = "percentage",
                 color = "group", 
                 palette = groups_label,
                 add = "jitter", 
                 ylab = plot_title,
                 short.panel.labs = TRUE, 
                 legend = NULL) +
    facet_wrap(~ interaction(marker, cluster),
    #facet_wrap(~ c(marker, cluster),
               ncol = length(markers_oi),
               #switch = "x",
               #labeller = layers_label
    ) +
    theme(legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          #strip.text.x = element_text(angle = 90), 
          axis.title.x = element_blank()) +
    #ylim(0,ylim_plot) +
    border()
  
  # Use only p.format as label. Remove method name.
  p <- p + stat_pvalue_manual(stat.test.signif,
                              label = "label",
                              tip.length = 0,
                              label.size = 7
  )
  
  p
  
  pdf(paste0(folder_plots, 
             "/perc_marker+_cells_", 
             paste(markers_oi, collapse = "_"), 
             "_in_", 
             paste(x, collapse = "_"), 
             "_boxplot.pdf"), 
      width = 15, 
      height = 20)
  print(p)
  dev.off()
  
  write_csv(stat.test, paste0(folder_tables, 
                              "/stats_perc_marker+_cells_", 
                              paste(markers_oi, collapse = "_"), 
                              "_in_", 
                              paste(x, collapse = "_"), 
                              ".csv")) 
  
  # print plot done
  print(paste0("Boxplot for selected markers ", paste(markers_oi, collapse = "_"), " was generated"))
      
  ##### Plot as HEATMAP #####
  
  # select only columns of interest
  dat_plot$count <- NULL
  
  # calculate average cell number per layer, per SAMPLE, per marker
  dat_avg <- aggregate(dat_plot$percentage,
                       list(dat_plot$group,
                            dat_plot$marker,
                            dat_plot$cluster), mean)

  # rename colnames
  colnames(dat_avg) <- c('group', 'marker' , 'cluster', 'count')
  
  # for loop to plot each marker separately
  for (marker in markers){
    #marker <- markers[1]
    
    # filter by marker 
    dat_sel <- dat_avg[dat_avg$marker == marker,]
    
    # reshape into wide table
    dat_sel <- pivot_wider(dat_sel, names_from = group, values_from = count)
    
    # re-order groups order
    # dat_sel <- dat_sel[, c("cluster", x)]
    
    # keep only clusters of interest
    # dat_sel <- dat_sel[(dat_sel$cluster %in% 
    #                                 c("cluster_2 - GFAP, S100B", 
    #                                   "cluster_3 - S100B", 
    #                                   "cluster_34 - GFAP, Iba1", 
    #                                   "cluster_35 - S100B", 
    #                                   "cluster_0 - OLIG2", 
    #                                   "cluster_30 - OLIG2, S100B", 
    #                                   "cluster_15 - CD68, GAD1",
    #                                   "cluster_24 - Iba1, CD68", 
    #                                   "cluster_12 - Iba1")
    # ),]
    
    # define order of clusters (rows) 
    clusters_order <- c("cluster_0",
                        "cluster_1",
                        "cluster_5",
                        "cluster_16",
                        "cluster_15",
                        "cluster_10",
                        "cluster_8",
                        "cluster_9",
                        "cluster_17",
                        "cluster_3",
                        "cluster_7",
                        "cluster_11",
                        "cluster_14"
                        )
    
    clusters_order <- cluster_labels$full_label[match(clusters_order, cluster_labels$cluster)]
    
    # re-order rows
    dat_sel <- dat_sel %>% slice(match(clusters_order, dat_sel$cluster))
    
    # colnames(dat_sel)[3:4] <- c('FAD2m', 'FAD4m')
    # dat_sel$FAD2m <- as.numeric(round(dat_sel$FAD2m, 2))
    # dat_sel$FAD4m <- as.numeric(round(dat_sel$FAD4m, 2))
    
    # make into matrix for heatmap plotting
    m <- as.matrix(dat_sel[, -1])

    # remove cluster column
    m <- m[,2:3]
    
    # make sure the matrix is a matrix and numbers are numeric
    library(pheatmap)
    m <- apply(m, 2, as.numeric)
    
    # define rownames of heatmap
    rownames(m) <- dat_sel$cluster
    
    # generate cluster labels for annotation_row
    plot_labels <- cluster_labels[, c(3,4)] 
    plot_labels <- plot_labels[!plot_labels$full_label %in% c('positive_for_all', 'negative_for_all'), ]
    rownames(plot_labels) <- plot_labels$full_label
    plot_labels[, 1] <- NULL
    colnames(plot_labels) <- "assigned_cluster"
    
    # define order of clusters
    plot_labels$assigned_cluster <- factor(plot_labels$assigned_cluster,
                                              levels = c('microglia', "astrocytes", "oligos", 
                                                         "uncl_neurons", "ex_neurons", "inh_neurons"))
    
    # define colors for legend
    colors = list(
      assigned_cluster = c(
        ex_neurons = "#f66d9b", 
        inh_neurons = "#9561e2",
        uncl_neurons = "#e3342f", 
        astrocytes = "#38c172",
        microglia = "#ffed4a", 
        oligos = "#4dc0b5"
        )
    )
    
    
    # define legend breaks
    
    # breaks for Ab (0.03): c(seq(0, 40, by=0.1), seq(40.1, 100, by=10)) # glia
    # breaks for pTau (0.2): c(seq(0, 0.5, by=0.05), seq(0.6, 5, by=0.5)) # glia
    
    # if (analysis == "Ab"){
    #   my.breaks <- c(seq(0, 40, by=0.1), seq(40.1, 100, by=10)) # glia
    #   
    # } else if (analysis == "pTau") {
    #   my.breaks <- c(seq(0, 0.5, by=0.05), seq(0.6, 5, by=0.5)) # glia
    #   
    # } else {
    #   my.breaks <- c(seq(0, 0.5, by=0.05), seq(0.6, 5, by=0.5)) # glia
    # }
    # 
    # # my.colors <- c(colorRampPalette(colors = c("yellow", "green"))(length(my.breaks)/2), 
    # #                #colorRampPalette(colors = c("white", "yellow"))(length(my.breaks)/3),
    # #                colorRampPalette(colors = c("yellow", "red"))(length(my.breaks)/2))
    # 
    # my.colors <- c(colorRampPalette(colors = c(viridis_pal()(3)[1], viridis_pal()(3)[2]))(length(my.breaks)/2), 
    #                #colorRampPalette(colors = c("white", "yellow"))(length(my.breaks)/3),
    #                colorRampPalette(colors = c( viridis_pal()(3)[2], viridis_pal()(3)[3]))(length(my.breaks)/2))
    # 
    
    # plot heatmap
    plot <- pheatmap(m, 
                     cluster_cols = FALSE, 
                     cluster_rows = FALSE, 
                     angle_col = 0, 
                     #color = my.colors,
                     #breaks = my.breaks,
                     annotation_row = plot_labels,
                     annotation_colors = colors,
                     gaps_row = c(3, 7, 8, 10), 
                     border_color = NA,
                     main = paste0('% of ', marker, '+ cells per cluster'),
                     fontsize = 9
                     )
    plot
    
    # save heatmap
    pdf(paste0(folder_plots, marker, "+_cells_per_cluster_in_", paste(x, collapse = "_"), "_heatmap.pdf"), 
        height = 5, 
        width = 6)
    print(plot)
    dev.off()
    
  }

  print(paste0(" ----- END CONDITION #", cond, " ------"))
  
  
} 



##### 6. Calculate tot marker+ cells per group, regardless of cluster #####

# filter for clusters of interest (neuronal clusters in this case)
combined_data_sel <- combined_data[!(combined_data$cluster %in% c("2", "3", "4", "6", "7", "11", "12", "13", "14", "17")),]

# select columns with cell numbers only
combined_data_sel <- combined_data_sel %>% select(!starts_with('perc_'))

# turn long
combined_data_sel <- combined_data_sel %>% pivot_longer(
  cols = all_of(markers_oi), # Columns to pivot (e.g., all columns starting with "year")
  names_to = "marker",          # New column name for the column names
  values_to = "cell_number"         # New column name for the values
  )

# calculate total marker+ cell per sample
postive_cells_tot <- combined_data_sel %>% dplyr::group_by(sample_name, marker) %>% dplyr::summarise(cell_number = sum(cell_number),
                                                                                                     tot = sum(tot)
                                                                                                     )
# calculate percentage
postive_cells_tot$cell_perc <- (postive_cells_tot$cell_number / postive_cells_tot$tot) *100

# add group
postive_cells_tot$group <- purrr::map_chr(postive_cells_tot$sample_name, ~strsplit(., "_")[[1]][1])

# save
dir.create(paste0("tables/clusters/cells_positive_for_markers/tot_cell_number/"))
write_csv(postive_cells_tot, "tables/clusters/cells_positive_for_markers/tot_cell_number/tot_cells_positive_for_markers.csv")

##### 7. Statistical analysis + plotting ######
##### 7.1 Add variants ######

# read in list of sample/variants
variants <- read_csv("input/vars_per_sample.csv")

# add empty columns in dat and dat_long
postive_cells_tot$TREM2var <- "a"
postive_cells_tot$ApoE4pos <- "b"
postive_cells_tot$AD_stage <- "c"
postive_cells_tot$TREM2var_ADstage <- "d"
postive_cells_tot$ApoE4pos_TREM2var <- "e"

# get all sample names
sample_names <- unique(postive_cells_tot$sample_name)

# assign variants to samples
for (sample in 1:length(sample_names)){
  # assign ApoE variant
  postive_cells_tot[postive_cells_tot$sample_name == sample_names[sample],]$ApoE4pos <- 
    variants[variants$sample_name == sample_names[sample],]$ApoE4pos
  
  # assign TREM2 variant
  postive_cells_tot[postive_cells_tot$sample_name == sample_names[sample],]$TREM2var <- 
    variants[variants$sample_name == sample_names[sample],]$TREM2var
  
  # assign AD_stage variant
  postive_cells_tot[postive_cells_tot$sample_name == sample_names[sample],]$AD_stage <-
    variants[variants$sample_name == sample_names[sample],]$AD_stage
  
  # assign TREM2var_ADstage variant
  postive_cells_tot[postive_cells_tot$sample_name == sample_names[sample],]$TREM2var_ADstage <-
    variants[variants$sample_name == sample_names[sample],]$TREM2var_ADstage
  
  # assign ApoE4pos_TREM2var variant
  postive_cells_tot[postive_cells_tot$sample_name == sample_names[sample],]$ApoE4pos_TREM2var <-
    variants[variants$sample_name == sample_names[sample],]$ApoE4pos_TREM2var
}
##### 7.2 Define comparisons and folders (OPEN -->) #####

# Create table with conditions for for loop 
conditions <- data.frame(matrix(ncol = 3, nrow = 3))
conditions <- data.frame(
  cond_name = c(1, 2, 3),
  order = I(list(c("CtrlCV", "AlzCV"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzTREM2"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzR62H", "Alz472H"))),
  group_column = c("group", "group", "TREM2var")
)

# define column with value to test
column_oi <- 'cell_perc'

# define if value is a percentage
is_percentage <- TRUE

# define paths to folders
folder <- 'clusters'
subfolder <- "cells_positive_for_markers"
subsubfolder <- 'tot_cell_number'

# make directories
dir.create(paste0("plots/", folder))
dir.create(paste0("tables/", folder))

dir.create(paste0("plots/", folder, "/", subfolder))
dir.create(paste0("tables/", folder, "/", subfolder))

dir.create(paste0("plots/", folder, "/", subfolder, "/", subsubfolder))
dir.create(paste0("tables/", folder, "/", subfolder, "/", subsubfolder))

folder_tables <- paste0("tables/", folder, "/", subfolder, "/", subsubfolder, "/")
folder_plots <- paste0("plots/", folder, "/", subfolder, "/", subsubfolder, "/")

##### 7.3 For loop: stat analysis and plot ####

for (cond in 1:nrow(conditions)) {
  #cond = 1
  
  print(paste0(" ----- START CONDITION #", cond, " ------"))
  
  dat <- as.data.frame(postive_cells_tot)
  x <- unlist(conditions$order[cond])
  
  if ((conditions$group_column[cond]) != "group") {
    dat$group <- NULL
    colnames(dat)[which(colnames(dat) == conditions$group_column[cond])] <- 'group'
  }
  
  dat <- dat[dat$group %in% x,]
  
  ##### Save file for paper #####
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
  dat <- dat[, which(colnames(dat) %in% c('group', 'marker', column_oi))]
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
  
  ##### Wilcox test (when only two groups are present) #####
  
  # generate stat.test file to add tukey p.adj of count, for plotting
  stat.test <- data.frame(matrix(ncol = 7))
  colnames(stat.test) <- c('group1', 'group2', 'marker', 'p.adj', 'method', 'anova_krus', 'y.position')
  
  # independent 2-group Mann-Whitney U Test
  # wilcox.test(y~A)
  # where y is numeric and A is A binary factor
  
  if (length(groups) == 2){
    
    stat.test_temp <- stat.test
    
    for (marker in unique(dat$marker)) {
      #marker = 'Ab4G8'
      samples <- dat[dat$marker == marker, ]
      
      wt_res <-  wilcox.test(samples$count ~ samples$group)
      
      stat.test_temp$group1 <- groups[1]
      stat.test_temp$group2 <- groups[2]
      stat.test_temp$marker <- marker
      stat.test_temp$p.adj <- wt_res$p.value
      stat.test_temp$method <- "wilcox.test"
      stat.test_temp$anova_krus <- "N/A"
      stat.test_temp$y.position <- max(samples$percentage)*1.1
      
      stat.test <- rbind(stat.test, stat.test_temp)
      
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
  
  ##### Test normality (when more than two groups are present) ########
  # testing normality with a Shapiro-Wilk Test.
  # source: https://www.statology.org/test-for-normality-in-r/
  
  # test normality 
  # group_by: group and marker 
  
  if (length(groups) > 2){
    
    #prepare data table for result of normality test
    tests_todo <- data.frame(matrix(ncol = 3))
    colnames(tests_todo) <- c('group', 'marker', 'dist')
    
    tests_todo_temp <- tests_todo
    
    for (group in 1:length(groups)){
      skip_to_next <- FALSE
      tryCatch({
        #group = 1
        samples <- dat[dat$group == groups[group], ]
        
        for (marker in unique(dat$marker)){
          #marker = unique(dat$marker)[3]
          samples1 <- samples[samples$marker == marker, ]
          
          if ((shapiro.test(samples1$count)$p.value) < 0.05) {
            print(paste0("Marker ", marker,' in ', groups[group], " is **NOT** normally distributed"))
            
            tests_todo_temp$group <- groups[group]
            tests_todo_temp$marker <- marker
            tests_todo_temp$dist <- 'kruskal.wilcox'
            
            tests_todo <- rbind(tests_todo, tests_todo_temp)
            
          } 
          
          if ((shapiro.test(samples1$count)$p.value) > 0.05) {
            print(paste0("Marker ", marker,' in ', groups[group], " is normally distributed"))
            
            tests_todo_temp$group <- groups[group]
            tests_todo_temp$marker <- marker
            tests_todo_temp$dist <- 'anova.tukey'
            
            tests_todo <- rbind(tests_todo, tests_todo_temp)
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
    tests_todo_final <- data.frame(matrix(ncol = 2))
    colnames(tests_todo_final) <- c('marker', 'test')
    
    tests_todo_final_temp <- tests_todo_final
    
    # check which layers are normally distributed and define stat test to use for analysis
    for (marker in unique(dat$marker)){
      
      #marker <- unique(dat$marker)[1]
      samples1 <- samples[tests_todo$marker == marker, ]
      
      if (nrow(samples) > 1) {
        if (all(samples$dist == "anova.tukey")){
          tests_todo_final_temp$marker <- marker
          tests_todo_final_temp$test <- "anova.tukey"
        } else {
          tests_todo_final_temp$marker <- marker
          tests_todo_final_temp$test <- "kruskal.wilcox"
        }
      } else {
        tests_todo_final_temp$marker <- marker
        tests_todo_final_temp$test <- "kruskal.wilcox"
      }
      tests_todo_final <- rbind(tests_todo_final, tests_todo_final_temp)
      
    }
    
    # remove first row
    if (is.na(tests_todo_final[1,1])){
      tests_todo_final <- tests_todo_final[-1,] 
    }
    tests_todo_final <- as.data.frame(tests_todo_final)
    
    print("Normality test was performed as more than two groups are being tested")
    rm(stat.test_temp)
    
  } else {
    print("Normality test *NOT* applied as groups analysed are only 2")
  }
  
  rm(tests_todo_final_temp)
  rm(tests_todo_temp)
  
  # only some data are normally distributed:
  # --> ANOVA + Tukey on normally distributed data
  # --> Kruskal-wallis + Wilcoxon test on NOT normally distributed data
  
  ##### ANOVA + Tukey (for normally distributed data) #####
  # source: https://www.scribbr.com/statistics/anova-in-r/#:~:text=We%20can%20perform%20an%20ANOVA,levels%20of%20the%20independent%20variable.
  
  if (length(groups) > 2){
    
    # filter samples for ANOVA
    temp <- filter(tests_todo_final, test == 'anova.tukey')
    
    marker_anova <- temp$marker
    
    samples_anova <- dat %>% filter(marker %in% marker_anova)
    
    # check we have all cells
    sample_names <- unique(dat$sample_name)
    nrow(samples_anova) == length(marker_anova)*length(sample_names)
    
    # filter samples per layer
    if (nrow(samples_anova) > 0){
      
      for (marker in unique(samples_anova$marker)){
        #cluster = unique(dat$cluster)[1]
        samples_anova1 <- samples_anova[samples_anova$marker == marker, ]
        
        #### ANOVA + Tukey test
        # save ANOVA results
        one.way <- aov(count ~ group, data = samples_anova1)
        
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
        tukey.one.way.df$marker <- marker
        tukey.one.way.df <- tukey.one.way.df %>% relocate(marker, .after = group2)
        
        # add y.position and Pr(>F)
        y_increase <- max(samples_anova$percentage)*0.2
        for (n in 1:nrow(tukey.one.way.df)){
          tukey.one.way.df$y.position[n] <- max(samples_anova$percentage)+(y_increase+((y_increase/2)*(n-1)))
          tukey.one.way.df$anova_krus[n] <- summary(one.way)[[1]][["Pr(>F)"]][1]
        }
        stat.test <- rbind(stat.test, tukey.one.way.df)
      }
      
      print("ANOVA + Tukey tests performed on multiple groups")
      
    } else {
      print("No samples for ANOVA")
    }
  } else {
    print("ANOVA + Tukey tests not performed as only two groups are tested")
  }
  
  ##### Kruskal + Wilcoxon (for not normally distributed data) #####
  
  if (length(groups) > 2){
    
    # filter samples for Kruskal
    temp <- filter(tests_todo_final, test == 'kruskal.wilcox')
    
    marker_kruskal <- unique(temp$marker)
    
    samples_kruskal <- dat %>% filter(marker %in% marker_kruskal)
    
    # check we have all cells
    sample_names <- unique(dat$sample_name)
    nrow(samples_kruskal) == length(marker_kruskal)*length(sample_names)
    
    if (nrow(samples_kruskal) > 0){
      
      for (marker in unique(samples_kruskal$marker)){
        #cluster = unique(dat$cluster)[1]
        samples_kruskal1 <- samples_kruskal[samples_kruskal$marker == marker, ]
        
        # perform Kruskal test and save
        kruskal.test <- kruskal.test(count ~ group, data = samples_kruskal1)
        
        # perform Wilcox test and save
        wilcox.test <- pairwise.wilcox.test(samples_kruskal1$count, samples_kruskal1$group,
                                            p.adjust.method = "BH", 
                                            exact = FALSE)
        
        # modify wilcox.test to fit stat.test
        wilcox.test <- rownames_to_column(as.data.frame(wilcox.test$p.value))
        wilcox.test <- melt(wilcox.test)
        colnames(wilcox.test) <- c('group2', 'group1', 'p.adj')
        wilcox.test$method <- 'kruskal.wilcox'
        
        wilcox.test <- wilcox.test %>% relocate(group2, .after = group1)
        wilcox.test$marker <- marker
        wilcox.test <- wilcox.test %>% relocate(marker, .after = group2)
        
        # add y.position and Pr(>F) and remove rows from same samples
        y_increase <- max(samples_kruskal1$percentage)*0.2
        for (n in 1:nrow(wilcox.test)){
          wilcox.test$y.position[n] <- max(samples_kruskal1$percentage)+(y_increase+((y_increase/2)*(n-1)))
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
  
  dat_plot <- dat
  
  # re-order groups order
  dat_plot <- dat_plot %>% arrange(sapply(group, function(y) which(y == x)))
  
  # filter p-adj that are not significant !!!1 DO NOT RUN UNTIL YOU HAVE SIGNIFICANT MARKERS !!!
  #stat.test.signif <- filter(stat.test, p.adj <= 0.05)
  stat.test.signif <- stat.test # use this instead
  
  # ylim 
  ylim_plot <- max(stat.test.signif$y.position)
  
  # set as dataframe 
  dat_plot <- as.data.frame(dat_plot)
  
  # define colors for each group
  groups_label <- c(
    "AlzCV" = "#D02D27",
    "AlzTREM2" = "#D39898",
    "AlzR47H" = "#D39898",
    "AlzR62H" = "#D8C5C5",
    "CtrlCV" = "#314B93",
    "CtrlTREM2" = "#8297C4"
  )
  
  # # change y.positions (only if necessary)
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
  
  # Box plot facetted by "marker"
  p <- ggboxplot(dat_plot, x = "group", y = "percentage",
                 color = "group", 
                 palette = groups_label,
                 add = "jitter", 
                 ylab = '% of marker+ cells',
                 short.panel.labs = TRUE, 
                 legend = NULL) +
    facet_wrap(~ marker,
               ncol = length(markers_oi),
               #switch = "x",
               #labeller = layers_label
    ) +
    theme(legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          #strip.text.x = element_text(angle = 90), 
          axis.title.x = element_blank()) +
    #ylim(0,ylim_plot) +
    border()
  
  # Use only p.format as label. Remove method name.
  p <- p + stat_pvalue_manual(stat.test.signif,
                              label = "label",
                              tip.length = 0,
                              label.size = 7
  )
  
  p
  
  pdf(paste0(folder_plots, 
             "/cell_positive_for_", 
             paste(markers_oi, collapse = "_"), 
             "_in_", 
             paste(x, collapse = "_"), 
             ".pdf"), 
      width = 7, 
      height = 3)
  print(p)
  dev.off()
  
  write.xlsx(stat.test, paste0(folder_tables, 
                              "/stats_cell_positive_for_", 
                              paste(markers_oi, collapse = "_"), 
                              "_in_", 
                              paste(x, collapse = "_"), 
                              ".xlsx")) 
  
  # print plot done
  print(paste0("Plot for selected markers ", paste(markers_oi, collapse = "_"), " was generated"))
  
  
  print(paste0(" ----- END CONDITION #", cond, " ------"))
  
}
