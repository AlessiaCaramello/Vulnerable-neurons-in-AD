# Calculating number and density of cells per cluster

# input data: 
# - clustered_cells.csv (SIMPLI output)

####   CELL NUMBER PER CLUSTER ####
###### 0a. Load libraries ######

#install.packages("tidyverse")
library(tidyverse)
library(reshape)
library(data.table)
library(plyr)
library(dplyr)
library(tibble)
#install.packages("progress")
library(progress)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggsignif)
#install.packages("ggpubr")
library(ggpubr)
library(rstatix)
#install.packages("DirichletReg")
library(DirichletReg)
#install.packages("pheatmap")
library(pheatmap)
library(ggsci)
library(RColorBrewer)
library(viridis)
#install.packages("writexl")
library(writexl)
#install.packages("stringr")
library(stringr)

# set working directory
setwd("~/UK Dementia Research Institute Dropbox/Alessia Caramello/Alessia Lab Dropbox/Projects/UKDRI project/Image analysis/HistoCat:SIMPLI analysis/051022 TREM2 cohort/R code")

# create output folder
# dir.create("plots/clusters/dirichlet")
# dir.create("tables/clusters/dirichlet")

###### 0b. Import table ######

# import file
clustered_cells <- read.csv("input/clustered_cells.csv")

# pick which clustering resolution to analyse
res <- 'res_0.9_ids'

# remove clusters negative for all and positive for all and synapses
clustered_cells <- clustered_cells[!clustered_cells[,res] %in% c(1,16,23),] 

###### 1a. Count number of occurrences (cells per cluster) #####

# extract column with Metadata_sample_name and clusters count
counts <- clustered_cells[, colnames(clustered_cells) %in% c('Metadata_sample_name', res)]

# add column with sample_name
counts <- add_column(counts, "sample_name" = stringr::str_extract(counts$Metadata_sample_name, "[^_]*_[^_]*"))
counts <- add_column(counts, "sample_replicate" = stringr::str_extract(counts$Metadata_sample_name, "[^_]+$"))
counts$Metadata_sample_name <- NULL

# change column name
colnames(counts)[which(colnames(counts) == res)] <- 'clusters'

# convert counts to data.table
counts <- as.data.frame(counts)

# count number of rows (= number of cells)
counted <- aggregate(counts$clusters, 
                     by=list( counts$sample_name, counts$sample_replicate, counts$clusters), 
                     FUN = length)

# change colnames
colnames(counted) <- c('sample_name', 'sample_replicate', 'cluster', 'cell_number')

# reshape
countedw <- reshape(counted, idvar = c("sample_name", 'sample_replicate'), 
                    timevar = "cluster", 
                    direction = "wide")

# replace NA with 0
countedw <- countedw %>% replace(is.na(.), 0)

# change colnames
colnames(countedw)[3:length(colnames(countedw))] <- 
  paste0("cluster_", gsub("cell_number.","",
                          colnames(countedw)[3:length(colnames(countedw))]))

###### 1b. Count tot cells per sample ######
countedw$sum <- rowSums(countedw[,3:ncol(countedw)])

# check we have the right amount of cells
sum(countedw$sum) == nrow(counts)
nrow(counts) == sum(counted$cell_number)
sum(countedw$sum) == sum(counted$cell_number)

# total cells per cluster
tot_clust <- as.data.frame(colSums(countedw[,3:ncol(countedw)]))
tot_clust$clusters <- rownames(tot_clust)
rownames(tot_clust) <- NULL
colnames(tot_clust)[1] <- 'cell_number'

# remove sum col, we don't need it anymore
countedw$sum <- NULL

# remove sum row, we don't need it anymore
tot_clust <- tot_clust[!tot_clust$clusters == 'sum', ]

###### 1c. Plot clusters total cell number ######

plot <- ggplot(data=tot_clust, 
       aes(y=cell_number, x=reorder(clusters, -cell_number)))+
  geom_boxplot()+
  xlab('Cluster #')+
  ylab('Cell number per cluster')+
  scale_y_continuous(breaks = pretty(tot_clust$cell_number, n = 10))+
  scale_x_discrete(labels=gsub("cluster_","",tot_clust$clusters))

plot

pdf(paste0("plots/clusters/clusters_tot_cells_", res, ".pdf"))
print(plot)
dev.off()

# save cells per cluster file
write_csv(countedw, paste0("tables/clusters/cells_per_cluster_", res, ".csv"))

###### 1d. Plot clusters cell number in individual samples and ROIs ######

dir.create("plots/clusters/cells_per_ROI")

counted$sample_replicate <- paste0("ROI#", counted$sample_replicate)

for (i in unique(counted$sample_name)){
  
  counted_plot <- counted[counted$sample_name == i, ]
  counted_plot$cluster <- paste0("cluster_", counted_plot$cluster) 
  
  plot <- ggboxplot(counted_plot, x = "sample_replicate", y = "cell_number",
                    color = "sample_replicate", 
                    ylab = paste0("total cells in individual ROIs in sample ", i) ,
                    xlab = "",
                    short.panel.labs = TRUE) +
    facet_wrap(~cluster) +
    theme(
      legend.title = element_blank(), 
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1))
  
  plot
  
  pdf(paste0("plots/clusters/cells_per_ROI/cells_per_cluster_per_ROI_in_", i, ".pdf"))
  print(plot)
  dev.off()

}

####   DIRICHLET OF CELL NUMBER PER CLUSTER #####
###### 2a. Adapt tables to Dirichlet analysis ######

# remove sum col, we don't need it anymore
countedw$sum <- NULL

# turn countedw in long table
countedl <- countedw %>% pivot_longer(cols = starts_with("cluster_"),
                                      names_to = 'mark',
                                      values_to = 'value') 

# add group column in both countedl and countedw
countedw <- add_column(countedw, "group" = stringr::str_extract(countedw$sample_name, "[^_]*"))
countedl <- add_column(countedl, "group" = stringr::str_extract(countedl$sample_name, "[^_]*"))

# now we have both wide and long table

# get clusters/mark names
cell_cols <-
  colnames(countedw)[!colnames(countedw) %in% c('sample_name','sampleID',
                                      'group','disease', 'sample_replicate')]
# Save 
write_csv(countedw, paste0("tables/clusters/cells_per_cluster_wide.csv"))
write_csv(countedl, paste0("tables/clusters/cells_per_cluster_long.csv"))

###### 2b. Define comparison and create folders (OPEN -->) ########

# load input file
# countedl <- read_csv("tables/clusters/cells_per_cluster_long.csv")
# countedw <- read_csv("tables/clusters/cells_per_cluster_wide.csv")

# define comparisons
# 1st slot = group 1 name
# 2nd slot = group 2 name
# 3rd slot = column with into of group names

comp1 <- c('CtrlCV', 'AlzCV', 'group') 
comp2 <- c('CtrlCV', 'AlzTREM2', 'group') 
comp3 <- c('AlzCV', 'AlzTREM2', 'group') 

comp_all <- list(comp1, comp2, comp3)

# define paths and create folders
path_plots <- "plots/clusters/dirichlet/"
path_tables <- "tables/clusters/dirichlet/"

dir.create(path_plots)
dir.create(path_tables)

subfolder <- "cell_number"

dir.create(paste0(path_plots, subfolder))
dir.create(paste0(path_tables, subfolder)) 

path_plots <- paste0("plots/clusters/dirichlet/", subfolder, "/")
path_tables <- paste0("tables/clusters/dirichlet/", subfolder, "/")

###### 2c. Dirichlet and plot for loop for all comparisons defined ######

for (i in 1:length(comp_all)){
  #i=1
  ###### Prep tasks #####
  
  comp <- comp_all[[i]]

  #restore initial table
  dat <- countedw
  dat_long <- countedl
  
  # rename grouping column indicated in comp to "group"
  if (comp[3] != "group"){
    dat$group <- NULL
    dat_long$group <- NULL
    colnames(dat)[which(colnames(dat) == comp[3])] <- "group"
    colnames(dat_long)[which(colnames(dat_long) == comp[3])] <- "group"
  }
  
  dat <- filter(dat, dat$group %in% comp[1:2])
  dat_long <- filter(dat_long, dat_long$group %in% comp[1:2])
  
  # change comp
  #comp <- comp[1:2]

  ###### Dirichlet #######
  
  # keep only useful columns 
  dat <- dat %>% dplyr::select(starts_with(c('group',"cluster_")))
  dat_long <- dat_long %>% dplyr::select(c('mark', 'value', 'group'))
  dat <- data.table(dat)
  dat_long <- as.data.table(dat_long)
  
  # make all numbers class numeric 
  dat[,2:length(dat)] <- lapply(dat[,2:length(dat)],as.numeric)
  
  
  #now pass to DirichletReg
  
  #create DirichletReg matrix from cell type proportions
  #will normalise for you
  DR_dat <- DirichletReg::DR_data(dat[,colnames(dat) %in% cell_cols, with=FALSE])
  #run DR
  fit <- DirichletReg::DirichReg(DR_dat ~ group, dat)
  #get p-values
  #function taken from scFLOW
  u <- summary(fit)
  pvals <- u$coef.mat[grep("Intercept", rownames(u$coef.mat), invert = TRUE), 4]
  dir <- u$coef.mat[grep("Intercept", rownames(u$coef.mat), invert = TRUE), 1]
  dir <- ifelse(dir>0,'Down','Up') #done in relation to disease ### change > to < id the direction is wrong
  v <- names(pvals)
  pvals <- matrix(pvals, ncol = length(u$varnames))
  rownames(pvals) <- gsub(as.name(dat$group), "", v[1:nrow(pvals)])
  colnames(pvals) <- u$varnames
  #do same for dir
  v <- names(dir)
  dir <- matrix(dir, ncol = length(u$varnames))
  rownames(dir) <- gsub(as.name(dat$group), "", v[1:nrow(dir)])
  colnames(dir) <- u$varnames
  
  # correction for multiple testing
  pvals_adj <- t(matrix(p.adjust(pvals,method='fdr')))
  colnames(pvals_adj) <- colnames(pvals)
  rownames(pvals_adj) <- rownames(pvals)
  
  # add the stars
  # run pvals for seeing the table 
  pvals <- as.data.frame(t(pvals_adj)) %>% # for multiple comparison: pvals <- as.data.frame(t(pvals_adj))
    tibble::rownames_to_column(var = 'celltype') %>%
    tidyr::pivot_longer(
      cols = rownames(pvals),
      names_to = 'group',
      values_to = "pval"
    ) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        pval <= 0.001 ~ "***",
        pval <= 0.01 ~ "**",
        pval <= 0.05 ~ "*",
        pval > 0.05 ~ "",
        is.na(pval) ~ ""
      )
    ) %>%
    as.data.frame()
  pvals$direction <- t(dir)
  #rename cell types to just mark used
  pvals$celltype <- gsub(".*_","",pvals$celltype)
  #Up and down is in relation to AD so Up would be the number of cells increased
  
  ###### Plot and save #####
  
  #plot with significance mark
  library(ggplot2)
  library(cowplot)
  pal <- c("#F5CDB4","#78A2CC")#"#FDDDA0")
  #add maxcell type prop to pvals for plotting
  pvals <- as.data.table(pvals)
  pvals[,group:=comp[2]]
  mx_scores <- dat_long[,.(max_val=max(value)),by=mark]
  setkey(mx_scores,mark)
  setkey(pvals,celltype)
  pvals[mx_scores,max_val:=i.max_val]
  mx_scores <- mx_scores[order(mx_scores$max_val, decreasing=TRUE),] # reorder mx_scores based on max_val
  level_order <- mx_scores$mark # define order of plot based on mx_scores$mark (ordered for max_val)
  
  #keep label for just disease
  dat_long$value <- as.numeric(dat_long$value)
  
  plot <- ggplot(data=dat_long,
                 aes(fill = group,y = value, x=mark))+ 
    geom_boxplot(outlier.colour = 'white',outlier.size = 0)+
    #significance from DirichletReg
    geom_text(data=pvals,aes(x=celltype, y=as.numeric(max_val),label=label), 
              col='red', size=5, nudge_x = -.05)+ #remove angle = 90 if not flipping the plot
    #coord_flip()+ # for flipping the plot
    cowplot::theme_cowplot()+
    scale_fill_manual(values=pal)+
    theme(legend.title = element_blank(), axis.text.x = element_text(size = 9, angle = 45, hjust = 1))+
    xlab('')+ylab('# of cells per cluster')+
    scale_x_discrete(limits = level_order) 
  
  pdf(paste0(path_plots, "/dirichlet_cell_number_plot_", comp[2], "_vs_", comp[1], ".pdf"))
  print(plot)
  dev.off()
  
  # save p values file
  write_csv(pvals, paste0(path_tables, "/dirichlet_cell_number_stats_", comp[2], "_vs_", comp[1], ".csv"))
  
}

####   CELL DENSITY PER CLUSTER ####
###### 3a. Calculate total area and cell density #####

# find max X and Y for each ROI
groups <- unique(clustered_cells$Metadata_sample_name)

countedl <- read_csv("tables/clusters/cells_per_cluster_long.csv")
countedl$max_length <- NA
countedl$max_width <- NA
countedl$sample_replicate <- paste0(countedl$sample_name , "_", countedl$sample_replicate)

for (i in 1:length(groups)) {
  #i=1
  
  # subset rows of each ROI
  rows <- clustered_cells[clustered_cells$Metadata_sample_name == groups[i],]
  
  # identify max value for X and Y
  max_x <- max(rows$Location_Center_X, na.rm = TRUE)
  max_y <- max(rows$Location_Center_Y, na.rm = TRUE)
  
  countedl_id <- countedl$sample_replicate == groups[i]
  countedl$max_length[countedl_id] <- max_x
  countedl$max_width[countedl_id] <- max_y
  
}

# calculate total area in mm2
countedl$um2_area <- countedl$max_length * countedl$max_width # um2
countedl$mm2_area <- countedl$um2_area / 1000000 # mm2
countedl$um2_area <- NULL

# calculate cell density (cell number / mm2)
countedl$cell_density <- as.numeric(countedl$value) / as.numeric(countedl$mm2_area)
countedl <- countedl %>% replace(is.na(.), 0) # remove NA

#save
write_csv(countedl, paste0("tables/clusters/cell_density_per_cluster_long.csv"))

# select columns of interest
countedl <- countedl[,c("sample_name","sample_replicate", "mark", "cell_density")] 

# rename cell_density column
colnames(countedl)[colnames(countedl) == "cell_density"] <- "value" 

# generate the wide format
countedw <- countedl %>% pivot_wider(names_from = mark, values_from = value)

#save
write_csv(countedw, paste0("tables/clusters/cell_density_per_cluster_wide.csv"))

####   DIRICHLET WITH CELL DENSITY PER CLUSTER ####
###### 4a. Define comparison and create folders (OPEN -->) ########

# load input files
# countedl <- read_csv("tables/clusters/cell_density_per_cluster_long.csv")
# countedw <- read_csv("tables/clusters/cell_density_per_cluster_wide.csv")

# define comparisons
# 1st slot = group 1 name
# 2nd slot = group 2 name
# 3rd slot = column with into of group names

comp1 <- c('CtrlCV', 'AlzCV', 'group') 
comp2 <- c('CtrlCV', 'AlzTREM2', 'group') 
comp3 <- c('AlzCV', 'AlzTREM2', 'group') 

comp_all <- list(comp1, comp2, comp3)

# define paths and create folders
path_plots <- "plots/clusters/dirichlet/"
path_tables <- "tables/clusters/dirichlet/"

dir.create(path_plots)
dir.create(path_tables)

subfolder <- "cell_density"

dir.create(paste0(path_plots, subfolder))
dir.create(paste0(path_tables, subfolder)) 

path_plots <- paste0("plots/clusters/dirichlet/", subfolder, "/")
path_tables <- paste0("tables/clusters/dirichlet/", subfolder, "/")

# get clusters/mark names
cell_cols <-
  colnames(countedw)[!colnames(countedw) %in% c('sample_name','sampleID',
                                                'group','disease', 'sample_replicate')]

###### 4b. Dirichlet and plot for loop for all comparisons defined ######

for (i in 1:length(comp_all)){
  #i=1
  ###### Prep tasks #####
  
  comp <- comp_all[[i]]
  
  #restore initial table
  dat <- countedw
  dat_long <- countedl
  
  # rename grouping column indicated in comp to "group"
  if (comp[3] != "group"){
    dat$group <- NULL
    dat_long$group <- NULL
    colnames(dat)[which(colnames(dat) == comp[3])] <- "group"
    colnames(dat_long)[which(colnames(dat_long) == comp[3])] <- "group"
  }
  
  dat <- filter(dat, dat$group %in% comp[1:2])
  dat_long <- filter(dat_long, dat_long$group %in% comp[1:2])
  
  # change comp
  #comp <- comp[1:2]
  
  ###### Dirichlet #######
  
  # keep only useful columns 
  dat <- dat %>% dplyr::select(starts_with(c('group',"cluster_")))
  dat_long <- dat_long %>% dplyr::select(c('mark', 'value', 'group'))
  dat <- data.table(dat)
  dat_long <- as.data.table(dat_long)
  
  # make all numbers class numeric 
  dat[,2:length(dat)] <- lapply(dat[,2:length(dat)],as.numeric)
  
  
  #now pass to DirichletReg
  
  #create DirichletReg matrix from cell type proportions
  #will normalise for you
  DR_dat <- DirichletReg::DR_data(dat[,colnames(dat) %in% cell_cols, with=FALSE])
  #run DR
  fit <- DirichletReg::DirichReg(DR_dat ~ group, dat)
  #get p-values
  #function taken from scFLOW
  u <- summary(fit)
  pvals <- u$coef.mat[grep("Intercept", rownames(u$coef.mat), invert = TRUE), 4]
  dir <- u$coef.mat[grep("Intercept", rownames(u$coef.mat), invert = TRUE), 1]
  dir <- ifelse(dir>0,'Down','Up') #done in relation to disease ### change > to < id the direction is wrong
  v <- names(pvals)
  pvals <- matrix(pvals, ncol = length(u$varnames))
  rownames(pvals) <- gsub(as.name(dat$group), "", v[1:nrow(pvals)])
  colnames(pvals) <- u$varnames
  #do same for dir
  v <- names(dir)
  dir <- matrix(dir, ncol = length(u$varnames))
  rownames(dir) <- gsub(as.name(dat$group), "", v[1:nrow(dir)])
  colnames(dir) <- u$varnames
  
  # correction for multiple testing
  pvals_adj <- t(matrix(p.adjust(pvals,method='fdr')))
  colnames(pvals_adj) <- colnames(pvals)
  rownames(pvals_adj) <- rownames(pvals)
  
  # add the stars
  # run pvals for seeing the table 
  pvals <- as.data.frame(t(pvals_adj)) %>% # for multiple comparison: pvals <- as.data.frame(t(pvals_adj))
    tibble::rownames_to_column(var = 'celltype') %>%
    tidyr::pivot_longer(
      cols = rownames(pvals),
      names_to = 'group',
      values_to = "pval"
    ) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        pval <= 0.001 ~ "***",
        pval <= 0.01 ~ "**",
        pval <= 0.05 ~ "*",
        pval > 0.05 ~ "",
        is.na(pval) ~ ""
      )
    ) %>%
    as.data.frame()
  pvals$direction <- t(dir)
  #rename cell types to just mark used
  pvals$celltype <- gsub(".*_","",pvals$celltype)
  #Up and down is in relation to AD so Up would be the number of cells increased
  
  ###### Plot and save #####
  
  #plot with significance mark
  library(ggplot2)
  library(cowplot)
  pal <- c("#F5CDB4","#78A2CC")#"#FDDDA0")
  #add maxcell type prop to pvals for plotting
  pvals <- as.data.table(pvals)
  pvals[,group:=comp[2]]
  mx_scores <- dat_long[,.(max_val=max(value)),by=mark]
  setkey(mx_scores,mark)
  setkey(pvals,celltype)
  pvals[mx_scores,max_val:=i.max_val]
  mx_scores <- mx_scores[order(mx_scores$max_val, decreasing=TRUE),] # reorder mx_scores based on max_val
  level_order <- mx_scores$mark # define order of plot based on mx_scores$mark (ordered for max_val)
  
  #keep label for just disease
  dat_long$value <- as.numeric(dat_long$value)
  
  plot <- ggplot(data=dat_long,
                 aes(fill = group,y = value, x=mark))+ 
    geom_boxplot(outlier.colour = 'white',outlier.size = 0)+
    #significance from DirichletReg
    geom_text(data=pvals,aes(x=celltype, y=as.numeric(max_val),label=label), 
              col='red', size=5, nudge_x = -.05)+ #remove angle = 90 if not flipping the plot
    #coord_flip()+ # for flipping the plot
    cowplot::theme_cowplot()+
    scale_fill_manual(values=pal)+
    theme(legend.title = element_blank(), axis.text.x = element_text(size = 9, angle = 45, hjust = 1))+
    xlab('')+ylab('# of cells per cluster')+
    scale_x_discrete(limits = level_order) 
  
  pdf(paste0(path_plots, "/dirichlet_cell_density_plot_", comp[2], "_vs_", comp[1], ".pdf"))
  print(plot)
  dev.off()
  
  # save p values file
  write_csv(pvals, paste0(path_tables, "/dirichlet_cell_density_stats_", comp[2], "_vs_", comp[1], ".csv"))
  
}

###### PLOT CLUSTERS SEPERATLY - STATISTIC BASED ON DIRICHLET CELL NUMBER PER CLUSTERS ######

# load input file
countedw <- read_csv("tables/clusters/cells_per_cluster_wide.csv")

# Create table with conditions for for loop 
conditions <- data.frame(matrix(ncol = 3, nrow = 3))
conditions <- data.frame(
  cond_name = c(1, 2, 3),
  order = I(list(c("CtrlCV", "AlzCV"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzTREM2"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzR62H", "Alz472H"))),
  group_column = c("group", "group", "TREM2var")
)

# define paths and create folders
path_plots <- "plots/clusters/dirichlet/"
path_tables <- "tables/clusters/dirichlet/"

dir.create(path_plots)
dir.create(path_tables)

subfolder <- "cell_number"

dir.create(paste0(path_plots, subfolder))
dir.create(paste0(path_tables, subfolder)) 

path_plots <- paste0("plots/clusters/dirichlet/", subfolder, "/")
path_tables <- paste0("tables/clusters/dirichlet/", subfolder, "/")

# Start for loop 

for (cond in 1:nrow(conditions)){
  #cond=1
  
  ##### Prep tasks #####
  
  # define sample name, column_oi and folder name
  column_oi <- conditions$group_column[cond] 
  x <- conditions$order[cond]
  sample_folder <- paste(x, collapse = "_") 
  
  # adjust dataframe
  dat <- as.data.frame(countedw)
  
  # rename grouping column indicated in comp to "group"
  if (column_oi != "group"){
    dat$group <- NULL
    colnames(dat)[which(colnames(dat) == column_oi)] <- 'group'
  }

  # create folder for samples
  dir.create(paste0(path_tables, sample_folder, "/"))
  dir.create(paste0(path_plots, sample_folder, "/"))
  
  # define which clusters we want to analyse
  clusters_oi <- colnames(dat %>% select(starts_with("cluster_"))) # using this as we're analysing all clusters
  # clusters_oi <- c("cluster_1", "cluster_2", "cluster_5", "cluster_8", "cluster_11", "cluster_15", "cluster_18")
  
  ##### Prepare table for plotting #####
  
  # filter only groups needed for comparison
  dat <- dat[dat$group %in% x,]
  
  # re-order groups order
  dat <- dat %>% arrange(sapply(group, function(y) which(y == x)))
  
  # select columns of with clusters of interest
  dat <- dat %>% dplyr::select(c("group", clusters_oi))
  
  # change colname
  #colnames(dat)[which(colnames(dat) == "name")] <- "clusters"
  
  # add unique identifier
  dat$ID <- 1:nrow(dat)
  dat$ID_group <- paste0(dat$ID, "_", dat$group)
  dat$ID <- NULL
  
  # turn long 
  dat_sel <- dat %>% pivot_longer(cols = names(dat)[grepl("cluster", names(dat))],
                                  names_to = 'clusters',
                                  values_to = 'cell_number')

  # Reordering group factor levels
  dat_sel$clusters <- factor(dat_sel$clusters,
                             levels = unique(dat_sel$clusters))
  
  ##### Save cell number file for paper #####
  
  # generate wide dat table
  dat_sel_wide <- dat_sel %>% pivot_wider(names_from = clusters, 
                                          values_from = cell_number)
  
  dat_sel_wide$ID_group <- NULL
  
  # save
  write_csv(dat_sel_wide, paste0(path_tables, sample_folder, "/cell_number_per_cluster.csv"))
  
  # calculate average and sd
  dat_sel$ID_group <- NULL 
  
  dat_sel_avg <- dat_sel %>% 
    dplyr::group_by(group, clusters) %>%
    dplyr::summarise(cells = paste(round(mean(cell_number), 2), 
                                   round(sd(cell_number), 2), 
                                   sep = ' ± ')
    )
  
  dat_sel_avg <- dat_sel_avg %>% pivot_wider(names_from = clusters, values_from = cells)
  
  write_xlsx(dat_sel_avg, paste0(path_tables, sample_folder, "/avg_cell_number_per_cluster.xlsx"))
  
  
  ##### Generate and save stat.test file (for ggplot) from Dirichlet analysis #####
  
  dat_sel$ID_group <- NULL
  
  # get list of files in folder of selected analysis
  files <- list.files(path = path_tables, 
                      pattern = "dirichlet_*")
  
  # create data frame to store all tables
  stat.test <- data.frame(matrix(ncol = 5))
  colnames(stat.test) <- c("clusters", "group1", "group2", "p.adj", "label")
  
  # loop through the files to import them, rbind and add comparison
  for (i in 1:length(files)){
    # read one file at time
    temp <- read.csv(paste0(path_tables, files[i]))
    # change colnames
    temp$direction <- NULL
    temp$max_val <- NULL
    colnames(temp)[1:3] <- c("clusters", "group1", "p.adj")
    # extract group2 comparison and insert into table
    group2 <- sub(paste0(".*", temp$group1[1], "_vs_"), "", files[i])
    group2 <- sub(paste0(".csv.*"), "", group2)
    temp <- add_column(temp,  group2 = group2, .after = "group1")
    # store into final dataframe
    stat.test <- rbind(stat.test, temp)
  }
  
  # remove first row
  if (is.na(stat.test[1,1])){
    stat.test <- stat.test[-1,]
  }
  rm(temp)
  
  # replace NA values with nothing, if any
  stat.test$label[is.na(stat.test$label)] <- ""
  
  # # filter only groups needed for comparison
  stat.test <- stat.test[stat.test$group1 %in% x,]
  stat.test <- stat.test[stat.test$group2 %in% x,]
  
  # keep only clusters of interest
  stat.test$clusters <- paste0("cluster_", stat.test$clusters)
  stat.test <- stat.test[stat.test$clusters %in% clusters_oi,]
  stat.test <- unique(stat.test)
  
  # save for paper
  # stat.test.paper <- stat.test[!stat.test$group1 %in% c("CtrlTREM2", "AlzTREM2"), ]
  # stat.test.paper <- stat.test.paper[!stat.test.paper$group2 %in% c("CtrlTREM2", "AlzTREM2"), ]

  write_xlsx(stat.test, paste0(path_tables, sample_folder, "/stats_cell_number_per_cluster.xlsx"))
  
  # filter p-adj that are not significant
  if (any(stat.test$p.adj <= 0.05)){
    stat.test.signif <- filter(stat.test, p.adj <= 0.05)
  } else {
    stat.test.signif <- stat.test[!is.na(stat.test$label),]
  }
  
  # add column for y.position
  stat.test.signif$y.position <- 0
  
  # replace NA labels with nothing
  #stat.test.signif$label[is.na(stat.test.signif$label)] <- ""
  
  # adjust y.position in stat.test.signif
  # get markers with significant difference
  markers <- unique(stat.test.signif$clusters)
  
  # keep only clusters with a statistical significance
  dat_sel <- dat_sel[dat_sel$clusters %in% markers, ]
  
  # create temporary stat.test.signif_temp file to store new y.position
  stat.test.signif_temp <- data.frame(matrix(ncol = ncol(stat.test.signif)))
  colnames(stat.test.signif_temp) <- colnames(stat.test.signif)
  
  # calculate new y.position
  for (marker in 1:length(markers)) {
    #marker <- 1
    temp <- dat_sel[str_detect(dat_sel$clusters, markers[marker]),]
    temp2 <- stat.test.signif[stat.test.signif$clusters == markers[marker],]
    for (i in 1:nrow(temp2)){
      temp2$y.position[i] <- max(temp$cell_number)+((max(dat_sel$cell_number)/10)*i)
    }
    stat.test.signif_temp <- rbind(stat.test.signif_temp, temp2)
  }
  
  # remove first row
  if (is.na(stat.test.signif_temp[1,1])){
    stat.test.signif_temp <- stat.test.signif_temp[-1,]
  }
  
  # re-assign to stat.test.signif
  stat.test.signif <- stat.test.signif_temp
  rm(stat.test.signif_temp)
  
  # ylim
  ylim_plot <- max(as.numeric(stat.test.signif$y.position))
  
  # change clusters labels
  labels <- read_csv("input/clusters_labels.csv")
  # labels$full_label <- paste0(labels$cluster, " - ", labels$label)
  # write_csv(labels, "input/clusters_labels.csv")
  cluster_labels <- labels[,1:3]
  cluster_labels <- cluster_labels[!cluster_labels$label == "positive_for_all",]
  cluster_labels <- na.omit(cluster_labels)
  
  # change cluster names (for plotting reasons)
  dat_sel$clusters <- cluster_labels$full_label[match(dat_sel$clusters, cluster_labels$cluster)]
  stat.test.signif$clusters <- cluster_labels$full_label[match(stat.test.signif$clusters, cluster_labels$cluster)]
  
  # remove duplicates
  stat.test.signif <- 
    stat.test.signif[!duplicated(apply(stat.test.signif[,1:3], 1, function(row) paste(sort(row),collapse=""))),]
  
  # re-order groups order
  stat.test.signif <- stat.test.signif %>% arrange(sapply(clusters, function(y) which(y == unique(dat_sel$clusters))))
  dat_sel <- dat_sel %>% arrange(sapply(group, function(y) which(y == x)))
  
  # define order for facet_wrap
  order <- unique(stat.test.signif$clusters)
  
  # round stat.test.signif
  stat.test.signif$p.adj <- round(stat.test.signif$p.adj, 4)
  
  ##### Plot ####
  #ylim_plot <- 500
  
  #SKIP FOR NOW
  # pick palette (when CtrlTREM2 samples are analysed separately)
  # if (any(x == "AlzTREM2")){
  #   print("AlzTREM2")
  #   sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898")
  # } else if (any(x %in% "AlzR47H")){
  #   print("AlzR47H")
  #   sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898", "#D8C5C5")
  # } else {
  #   print("AlzCV")
  #   sel_palette = c("#314B93", "#D02D27")
  # }
  # 
  
  
  # plot
  plot <- ggboxplot(dat_sel, x = "group", y = "cell_number",
                    color = "group", 
                    #palette = sel_palette,
                    add = "jitter", 
                    ylab = paste0("cell number per cluster") ,
                    xlab = "",
                    short.panel.labs = TRUE) +
    facet_wrap(~factor(clusters, order)) +
    theme(legend.title = element_blank(), axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) +
    ylim(NA,ylim_plot+5)+
    rremove("legend")
  
  # Add stats
  plot <- plot + stat_pvalue_manual(stat.test.signif,
                                    label = "label", 
                                    label.size = 5,
                                    tip.length = 0)
  plot
  
  # define plot size
  if (length(markers) == 1){
    width_plot = 2.5
    height_plot = 3
  } else if (length(markers) == 2) {
    width_plot = 5
    height_plot = 3
  } else if (length(markers) == 3) {
    width_plot = 7
    height_plot = 4
  }else if (length(markers) > 3 & length(markers) <= 6) {
    width_plot = 7
    height_plot = 5
  } else if (length(markers) > 6 & length(markers) <= 9) {
    width_plot = 7
    height_plot = 7
  } else if (length(markers) > 9 ) {
    width_plot = 10
    height_plot = 10
  } 
  
  pdf(paste0(path_plots, sample_folder, "/", 
             paste(parse_number(markers), 
                   collapse = "_"),
             "_clusters_cell_number_in_", paste(x, collapse = "_"), 
             ".pdf"), 
      width = width_plot, 
      height = height_plot)
  print(plot)
  dev.off()
  
}


###### PLOT CLUSTERS SEPERATLY - STATISTIC BASED ON DIRICHLET CELL DENSITY PER CLUSTER ######

# load input file
countedw <- read_csv("tables/clusters/cell_density_per_cluster_wide.csv")

# Create table with conditions for for loop 
conditions <- data.frame(matrix(ncol = 3, nrow = 3))
conditions <- data.frame(
  cond_name = c(1, 2, 3),
  order = I(list(c("CtrlCV", "AlzCV"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzTREM2"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzR62H", "Alz472H"))),
  group_column = c("group", "group", "TREM2var")
)

# define paths and create folders
path_plots <- "plots/clusters/dirichlet/"
path_tables <- "tables/clusters/dirichlet/"

dir.create(path_plots)
dir.create(path_tables)

subfolder <- "cell_density"

dir.create(paste0(path_plots, subfolder))
dir.create(paste0(path_tables, subfolder)) 

path_plots <- paste0("plots/clusters/dirichlet/", subfolder, "/")
path_tables <- paste0("tables/clusters/dirichlet/", subfolder, "/")

# Start for loop 

for (cond in 1:nrow(conditions)){
  #cond=1
  
  ##### Prep tasks #####
  
  # define sample name, column_oi and folder name
  column_oi <- conditions$group_column[cond] 
  x <- conditions$order[cond]
  sample_folder <- paste(x, collapse = "_") 
  
  # adjust dataframe
  dat <- as.data.frame(countedw)
  
  # rename grouping column indicated in comp to "group"
  if (column_oi != "group"){
    dat$group <- NULL
    colnames(dat)[which(colnames(dat) == column_oi)] <- 'group'
  }
  
  # create folder for samples
  dir.create(paste0(path_tables, sample_folder, "/"))
  dir.create(paste0(path_plots, sample_folder, "/"))
  
  # define which clusters we want to analyse
  clusters_oi <- colnames(dat %>% select(starts_with("cluster_"))) # using this as we're analysing all clusters
  # clusters_oi <- c("cluster_1", "cluster_2", "cluster_5", "cluster_8", "cluster_11", "cluster_15", "cluster_18")
  
  ##### Prepare table for plotting #####
  
  # filter only groups needed for comparison
  dat <- dat[dat$group %in% x,]
  
  # re-order groups order
  dat <- dat %>% arrange(sapply(group, function(y) which(y == x)))
  
  # select columns of with clusters of interest
  dat <- dat %>% dplyr::select(c("group", clusters_oi))
  
  # change colname
  #colnames(dat)[which(colnames(dat) == "name")] <- "clusters"
  
  # add unique identifier
  dat$ID <- 1:nrow(dat)
  dat$ID_group <- paste0(dat$ID, "_", dat$group)
  dat$ID <- NULL
  
  # turn long 
  dat_sel <- dat %>% pivot_longer(cols = names(dat)[grepl("cluster", names(dat))],
                                  names_to = 'clusters',
                                  values_to = 'cell_density')
  
  # Reordering group factor levels
  dat_sel$clusters <- factor(dat_sel$clusters,
                             levels = unique(dat_sel$clusters))
  
  ##### Save cell number file for paper #####
  
  # generate wide dat table
  dat_sel_wide <- dat_sel %>% pivot_wider(names_from = clusters, 
                                          values_from = cell_density)
  
  dat_sel_wide$ID_group <- NULL
  
  # save
  write_csv(dat_sel_wide, paste0(path_tables, sample_folder, "/cell_density_per_cluster.csv"))
  
  # calculate average and sd
  dat_sel$ID_group <- NULL 
  
  dat_sel_avg <- dat_sel %>% 
    dplyr::group_by(group, clusters) %>%
    dplyr::summarise(cells = paste(round(mean(cell_density), 2), 
                                   round(sd(cell_density), 2), 
                                   sep = ' ± ')
    )
  
  dat_sel_avg <- dat_sel_avg %>% pivot_wider(names_from = cluster, values_from = cells)
  
  write_xlsx(dat_sel_avg, paste0(path_tables, sample_folder, "/avg_cell_density_per_cluster.xlsx"))
  
  
  ##### Generate and save stat.test file (for ggplot) from Dirichlet analysis #####
  
  dat_sel$ID_group <- NULL
  
  # get list of files in folder of selected analysis
  files <- list.files(path = path_tables, 
                      pattern = "dirichlet_*")
  
  # create data frame to store all tables
  stat.test <- data.frame(matrix(ncol = 5))
  colnames(stat.test) <- c("clusters", "group1", "group2", "p.adj", "label")
  
  # loop through the files to import them, rbind and add comparison
  for (i in 1:length(files)){
    # read one file at time
    temp <- read.csv(paste0(path_tables, files[i]))
    # change colnames
    temp$direction <- NULL
    temp$max_val <- NULL
    colnames(temp)[1:3] <- c("clusters", "group1", "p.adj")
    # extract group2 comparison and insert into table
    group2 <- sub(paste0(".*", temp$group1[1], "_vs_"), "", files[i])
    group2 <- sub(paste0(".csv.*"), "", group2)
    temp <- add_column(temp,  group2 = group2, .after = "group1")
    # store into final dataframe
    stat.test <- rbind(stat.test, temp)
  }
  
  # remove first row
  if (is.na(stat.test[1,1])){
    stat.test <- stat.test[-1,]
  }
  rm(temp)
  
  # replace NA values with nothing, if any
  stat.test$label[is.na(stat.test$label)] <- ""
  
  # # filter only groups needed for comparison
  stat.test <- stat.test[stat.test$group1 %in% x,]
  stat.test <- stat.test[stat.test$group2 %in% x,]
  
  # keep only clusters of interest
  stat.test$clusters <- paste0("cluster_", stat.test$clusters)
  stat.test <- stat.test[stat.test$clusters %in% clusters_oi,]
  stat.test <- unique(stat.test)
  
  # save for paper
  # stat.test.paper <- stat.test[!stat.test$group1 %in% c("CtrlTREM2", "AlzTREM2"), ]
  # stat.test.paper <- stat.test.paper[!stat.test.paper$group2 %in% c("CtrlTREM2", "AlzTREM2"), ]
  
  write_xlsx(stat.test, paste0(path_tables, sample_folder, "/stats_cell_density_per_cluster.xlsx"))
  
  # filter p-adj that are not significant
  if (any(stat.test$p.adj <= 0.05)){
    stat.test.signif <- filter(stat.test, p.adj <= 0.05)
  } else {
    stat.test.signif <- stat.test[!is.na(stat.test$label),]
  }
  
  # add column for y.position
  stat.test.signif$y.position <- 0
  
  # replace NA labels with nothing
  #stat.test.signif$label[is.na(stat.test.signif$label)] <- ""
  
  # adjust y.position in stat.test.signif
  # get markers with significant difference
  markers <- unique(stat.test.signif$clusters)
  
  # keep only clusters with a statistical significance
  dat_sel <- dat_sel[dat_sel$clusters %in% markers, ]
  
  # create temporary stat.test.signif_temp file to store new y.position
  stat.test.signif_temp <- data.frame(matrix(ncol = ncol(stat.test.signif)))
  colnames(stat.test.signif_temp) <- colnames(stat.test.signif)
  
  # calculate new y.position
  for (marker in 1:length(markers)) {
    #marker <- 1
    temp <- dat_sel[str_detect(dat_sel$clusters, markers[marker]),]
    temp2 <- stat.test.signif[stat.test.signif$clusters == markers[marker],]
    for (i in 1:nrow(temp2)){
      temp2$y.position[i] <- max(temp$cell_density)+((max(dat_sel$cell_density)/10)*i)
    }
    stat.test.signif_temp <- rbind(stat.test.signif_temp, temp2)
  }
  
  # remove first row
  if (is.na(stat.test.signif_temp[1,1])){
    stat.test.signif_temp <- stat.test.signif_temp[-1,]
  }
  
  # re-assign to stat.test.signif
  stat.test.signif <- stat.test.signif_temp
  rm(stat.test.signif_temp)
  
  # ylim
  ylim_plot <- max(as.numeric(stat.test.signif$y.position))
  
  # change clusters labels
  labels <- read_csv("input/clusters_labels.csv")
  # labels$full_label <- paste0(labels$cluster, " - ", labels$label)
  # write_csv(labels, "input/clusters_labels.csv")
  cluster_labels <- labels[,1:3]
  cluster_labels <- cluster_labels[!cluster_labels$label == "positive_for_all",]
  cluster_labels <- na.omit(cluster_labels)
  
  # change cluster names (for plotting reasons)
  dat_sel$clusters <- cluster_labels$full_label[match(dat_sel$clusters, cluster_labels$cluster)]
  stat.test.signif$clusters <- cluster_labels$full_label[match(stat.test.signif$clusters, cluster_labels$cluster)]
  
  # remove duplicates
  stat.test.signif <- 
    stat.test.signif[!duplicated(apply(stat.test.signif[,1:3], 1, function(row) paste(sort(row),collapse=""))),]
  
  # re-order groups order
  stat.test.signif <- stat.test.signif %>% arrange(sapply(clusters, function(y) which(y == unique(dat_sel$clusters))))
  dat_sel <- dat_sel %>% arrange(sapply(group, function(y) which(y == x)))
  
  # define order for facet_wrap
  order <- unique(stat.test.signif$clusters)
  
  # round stat.test.signif
  stat.test.signif$p.adj <- round(stat.test.signif$p.adj, 4)
  
  ##### Plot ####
  #ylim_plot <- 500
  
  #SKIP FOR NOW
  # pick palette (when CtrlTREM2 samples are analysed separately)
  # if (any(x == "AlzTREM2")){
  #   print("AlzTREM2")
  #   sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898")
  # } else if (any(x %in% "AlzR47H")){
  #   print("AlzR47H")
  #   sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898", "#D8C5C5")
  # } else {
  #   print("AlzCV")
  #   sel_palette = c("#314B93", "#D02D27")
  # }
  # 
  
  
  # plot
  plot <- ggboxplot(dat_sel, x = "group", y = "cell_density",
                    color = "group", 
                    #palette = sel_palette,
                    add = "jitter", 
                    ylab = paste0("cell density per cluster") ,
                    xlab = "",
                    short.panel.labs = TRUE) +
    facet_wrap(~factor(clusters, order)) +
    theme(legend.title = element_blank(), axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) +
    ylim(NA,ylim_plot+5)+
    rremove("legend")
  
  # Add stats
  plot <- plot + stat_pvalue_manual(stat.test.signif,
                                    label = "label", 
                                    label.size = 5,
                                    tip.length = 0)
  plot
  
  # define plot size
  if (length(markers) == 1){
    width_plot = 2.5
    height_plot = 3
  } else if (length(markers) == 2) {
    width_plot = 5
    height_plot = 3
  } else if (length(markers) == 3) {
    width_plot = 7
    height_plot = 4
  }else if (length(markers) > 3 & length(markers) <= 6) {
    width_plot = 7
    height_plot = 5
  } else if (length(markers) > 6 & length(markers) <= 9) {
    width_plot = 7
    height_plot = 7
  } else if (length(markers) > 9 ) {
    width_plot = 10
    height_plot = 10
  } 
  
  pdf(paste0(path_plots, sample_folder, "/", 
             paste(parse_number(markers), 
                   collapse = "_"),
             "_clusters_cell_density_in_", paste(x, collapse = "_"), 
             ".pdf"), 
      width = width_plot, 
      height = height_plot)
  print(plot)
  dev.off()
  
}

