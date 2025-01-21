# Heatmap of markers expression in each cluster

# input data: 
# - clustered_cells.csv (SIMPLI output)

###### 0. Load libraries ######

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
setwd("~/R code")

# create output folder already done
# dir.create("plots/clusters")
# dir.create("plots/clusters/heatmap")
# dir.create("tables/clusters")

###### 0. Import table ######

# import file
clustered_cells <- read.csv("input/clustered_cells.csv")

# pick clustering resolution to analyse
res <- 'res_0.9_ids'

###### 1. Plot heatmap of clusters and markers expression ######

# select only columns of interest
intensity <- clustered_cells[, colnames(clustered_cells) %in% 
                                   c('Metadata_sample_name', 
                                     'CellName',
                                     res, 
                                     colnames(clustered_cells)[stringr::str_detect(colnames(clustered_cells),"Intensity_MeanIntensity_")]
                                   )]

# rename markers columns
colnames(intensity)[stringr::str_detect(colnames(intensity),"Intensity_MeanIntensity_")] <-
  str_remove(colnames(intensity)[stringr::str_detect(colnames(intensity),"Intensity_MeanIntensity_")], "Intensity_MeanIntensity_")

# calculate average intensity by cluster
intensity_mean <- intensity[, -(1:2)] # remove useless columns
colnames(intensity_mean)[1] <- "clusters" # rename cluster column
intensity_mean <- aggregate(intensity_mean[, 2:length(intensity_mean)], # calculate mean
                            list(intensity_mean$clusters), mean)
colnames(intensity_mean)[1] <- "clusters" # rename cluster column
intensity_mean$clusters <- paste0("cluster_", intensity_mean$clusters) # rename clusters

# remove markers that we don't need (DNA, AD markers, myelin and synapses)
intensity_mean <- intensity_mean[, -(grep("DNA1", colnames(intensity_mean)))] 
intensity_mean <- intensity_mean[, -(grep("pTau", colnames(intensity_mean)))] 
intensity_mean <- intensity_mean[, -(grep("APP", colnames(intensity_mean)))] 
intensity_mean <- intensity_mean[, -(grep("Ab", colnames(intensity_mean)))] 
intensity_mean <- intensity_mean[, -(grep("PLP1", colnames(intensity_mean)))] 
intensity_mean <- intensity_mean[, -(grep("Synapto", colnames(intensity_mean)))] 
intensity_mean <- intensity_mean[, -(grep("NTNG2", colnames(intensity_mean)))] 
#intensity_mean <- intensity_mean[, -(grep("Ab4G8", colnames(intensity_mean)))] 
#intensity_mean <- intensity_mean[, -(grep("H31L21", colnames(intensity_mean)))] 
#intensity_mean <- intensity_mean[, -(grep("LAMP1", colnames(intensity_mean)))] 
#intensity_mean <- intensity_mean[, -(grep("LC3AB", colnames(intensity_mean)))] 
#intensity_mean <- intensity_mean[, -(grep("MOAB2", colnames(intensity_mean)))] 
#intensity_mean <- intensity_mean[, -(grep("p62", colnames(intensity_mean)))] 
#intensity_mean <- intensity_mean[, -(grep("ubFK2", colnames(intensity_mean)))] 
#intensity_mean <- intensity_mean[, -(grep("ubK48", colnames(intensity_mean)))] 

# remove clusters positive or negative for all or expressing synapses markers only
intensity_mean <- intensity_mean[-which(intensity_mean$clusters %in% c("cluster_16")), ] 

# save for paper
write_csv(intensity_mean, "tables/clusters/intensity_mean.csv")

# define order of clusters (rows) 
clusters_order <- c("cluster_5", # excitatory
                    "cluster_10",
                    "cluster_12",
                    "cluster_13",
                    "cluster_18",
                    
                    "cluster_3", # inhibitory
                    "cluster_14",
                    "cluster_17",
                    "cluster_21",
                    "cluster_22",
                    
                    "cluster_4", # unclassified
                    "cluster_7",
                    "cluster_8",
                    
                    "cluster_9", # astro
                    "cluster_11", 
                    
                    "cluster_6", # microglia
                    "cluster_19",
                    
                    "cluster_0", # oligos
                    "cluster_2",
                    "cluster_15",
                    
                    "cluster_1", # positive_for_all
                    "cluster_23"
                    )

intensity_mean <- as.data.frame(intensity_mean[match(clusters_order, 
                                                     intensity_mean$clusters),]) # reorder dataframe

# define rownames of heatmap
rownames(intensity_mean) <- intensity_mean$clusters
intensity_mean$clusters <- NULL

# define order of markers (columns)
markers_order <- c("RORB", 
                   "FOXP2", "GPC5", "PCP4", "LMO3", "CUX2", "Calretinin", "LHX6", 
                   "ADARB1", "NPY", "SST", "PVALB", "VIP", "GAD1", "CCK", "CALB1",
                   "NeuN", "MAP2all", "MAP2","S100B", "GFAP",  "Iba1", "CD68", "OLIG2"
                   )

intensity_mean <- intensity_mean[,match(markers_order, colnames(intensity_mean))] # reorder dataframe

# make into matrix for heatmap plotting
m <- as.matrix(intensity_mean)

# read in labels and modify 
labels <- read_csv("input/clusters_labels.csv")

labels$full_label <- paste0(labels$cluster, " - ", labels$label)

cluster_labels <- as.data.frame(labels[, c(1,4)]) # for annotation_row
cluster_labels <- na.omit(cluster_labels)
rownames(cluster_labels) <- cluster_labels$cluster
cluster_labels[, 1] <- NULL
colnames(cluster_labels) <- "assigned_cluster"

markers_labels <- as.data.frame(labels[, 5:6]) # for annotation_col
#markers_labels <- markers_labels[1:17,] 
rownames(markers_labels) <- markers_labels$marker
markers_labels$marker <- NULL
colnames(markers_labels) <- "cell_type_marker"

# make as factor to define order
markers_labels$cell_type_marker <- factor(markers_labels$cell_type_marker,
                                          levels = c("inh_neurons", "ex_neurons", "microglia", 
                                                     "astrocytes", "oligos", "uncl_neurons"))

cluster_labels$assigned_cluster <- factor(cluster_labels$assigned_cluster,
                                          levels = c("positive_for_all", "microglia", "astrocytes", "oligos", 
                                                     "uncl_neurons", "ex_neurons", "inh_neurons"))
# define colors for legend
colors = list(
  cell_type_marker = c(
    ex_neurons = "#f66d9b", 
    inh_neurons = "#9561e2",
    uncl_neurons = "#3490dc", 
    astrocytes = "#38c172",
    microglia = "#ffed4a", 
    oligos = "#4dc0b5"), 
  assigned_cluster = c(
    ex_neurons = "#f66d9b", 
    inh_neurons = "#9561e2",
    uncl_neurons = "#3490dc", 
    astrocytes = "#38c172",
    microglia = "#ffed4a", 
    oligos = "#4dc0b5",
    positive_for_all = "#d5d5d5")
  )

# Plot

plot <- pheatmap(m,
                 scale = "column",
                 color = colorRampPalette(c("navy", "white", "red"))(50), 
                 annotation_col = markers_labels, 
                 annotation_row = cluster_labels,
                 annotation_colors = colors,
                 border_color = NA,
                 cluster_rows = FALSE, # remove dendograms
                 cluster_cols = FALSE, # remove dendograms
                 main = "markers expression in clusters",
                 gaps_col = c(6, 16, 19, 21, 23),
                 gaps_row = c(5, 10, 13, 15, 17, 20)
)
plot


# save heatmap
pdf(paste0("plots/clusters/heatmap/all_markers_expression_in_selected_clusters_w_all_positives.pdf"), 
    height = 10, 
    width =10)
print(plot)
dev.off()





