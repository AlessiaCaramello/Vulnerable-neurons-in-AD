# Spatial distirbution between clusters

# input data:
# - all_cells-clusters.csv (SIMPLI output)
# - UMAPs.rds (generated in 2b. UMAPs of clusters)
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
#BiocManager::install("SpatialExperiment")
library(SpatialExperiment)
#install.packages("Seurat")
#install.packages("devtools")
#devtools::install_github("mojaveazure/seurat-disk")
#devtools::install_github('satijalab/seurat-data')
#devtools::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(SeuratData) # for loading seurat object
library(SeuratDisk) # for saving seurat object
#BiocManager::install("imcRtools")
library(imcRtools)
library(scales)

# set working directory
setwd("~/R code")

# create output folder
# dir.create("plots/SIMPLI_analyses/clusters/interactions")
# dir.create("tables/SIMPLI_analyses/clusters/interactions")

##### 1. Start analysis ######

for (i in 1:4){
  #i=1
  if (i==1){
    # filter for excitatory neurons and microglia
    clusters_oi <- c("5", "8", "10", "11", "13", "16", "18", "23", "27",
                     "12", "24", "15")
    exclude <- c("12", "24", "15")
    exp <- "excitatory_microglia"
  }
  
  if (i==2){
  # filter for inhibitory neurons and microglia
  clusters_oi <- c("9", "14", "19", "21", "22", "25", "26", "28", "29", "31", "32", "33",
                   "12", "24", "15")
  exclude <- c("12", "24", "15")
  exp <- "inhibitory_microglia"
  }
  
  if (i==3){
  # filter for excitatory neurons and microglia
  clusters_oi <- c("5", "8", "10", "11", "13", "16", "18", "23", "27",
                   "2", "3", "34", "35")
  exclude <- c("2", "3", "34", "35")
  exp <- "excitatory_astrocytes"
  }
  
  if (i==4){
  # filter for inhibitory neurons and microglia
  clusters_oi <- c("9", "14", "19", "21", "22", "25", "26", "28", "29", "31", "32", "33",
                   "2", "3", "34", "35")
  exclude <- c("2", "3", "34", "35")
  exp <- "inhibitory_astrocytes"
  }
  
  print(paste0("START processing analysis: ", exp))
  
  ##### 1. import clusters and UMAP file and turn into Seurat object #####
  
  # read original all_cells-clusters.csv file
  dt <- read.csv("input/all_cells-clusters.csv")
  
  # read UMAPS file
  UMAPS <- readRDS("tables/clusters/UMAPs.rds")
  
  # filter for clusters of interest (excitatory neurons and microglia)
  dt_sel <- dt %>% filter(res_0.9_ids %in% clusters_oi)
  
  # generate seurat object
  mat <- dt_sel %>%
    select(c("CellName", starts_with("Intensity_MeanIntensity")))
  
  rownames(mat) <- mat$CellName
  mat$CellName <- NULL
  
  colnames(mat) <- gsub("Intensity_MeanIntensity_", "", colnames(mat))
  
  mat <- t(mat)
  metadata <- dt_sel %>%
    select(c("Metadata_sample_name", "CellName", "res_0.9_ids", 
             "Location_Center_X", "Location_Center_Y", "AreaShape_Area", 
             "AreaShape_Eccentricity", "AreaShape_MinorAxisLength", 
             "AreaShape_MaxFeretDiameter"
    ))
  rownames(metadata) <- metadata$CellName
  all(rownames(metadata) == colnames(mat))
  
  # filter for clusters of interest (excitatory neurons and microglia)
  UMAPS_sel <- UMAPS$res_0.9_ids
  UMAPS_sel <- UMAPS_sel %>% filter(res_0.9_ids %in% clusters_oi)
  umap_dt <- UMAPS_sel[, 1:2]
  
  colnames(umap_dt) <-  paste0("UMAP_", 1:2)
  umap_dt <- as.matrix(umap_dt)
  dimnames(umap_dt)[[1]] <- metadata$CellName
  
  seu <- CreateSeuratObject(counts = mat,
                            meta.data = metadata)
  seu <- NormalizeData(seu)
  seu[["umap"]] <- CreateDimReducObject(embeddings = umap_dt, key = "UMAP_")
  
  qs::qsave(seu, file = paste0("tables/clusters/interactions/IMC_seu_", exp, ".qs"))
  
  print(paste0("1/4 SEU object generated for: ", exp))
  
  ##### 2. Convert seu to spe #####
  
  # load seurat object
  #seu <- qs::qread(paste0("tables/SIMPLI_analyses/clusters/IMC_seu_", exp, ".qs"))
  
  # Extract spatial coordinates
  spatialCoords <- data.frame(as.matrix(
    seu@meta.data[, c("Location_Center_X", "Location_Center_Y")]))
  
  names(spatialCoords) <- c("Pos_X", "Pos_Y")
  
  # change sample_id
  seu$sample_id <- seu$Metadata_sample_name
  seu$Metadata_sample_name <- NULL
  
  # extract counts and metadata
  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")  
  metadata <- seu@meta.data
  
  spe <- SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = as.matrix(spatialCoords),
    colData = metadata
  )
  
  # test if it's a valid spatial object
  validObject(spe)
  
  # change some names
  names(spe@colData)[names(spe@colData) == "res_0.9_ids"] <- "cluster"
  names(spe@colData)[names(spe@colData) == "CellName"] <- "ObjectNumber"
  names(spe@colData)[names(spe@colData) == "AreaShape_Area"] <- "area"
  names(spe@colData)[names(spe@colData) == "AreaShape_Eccentricity"] <- "eccentricity"
  names(spe@colData)[names(spe@colData) == "AreaShape_MinorAxisLength"] <- "axis_minor_length"
  names(spe@colData)[names(spe@colData) == "AreaShape_MaxFeretDiameter"] <- "axis_major_length"
  spe$nCount_RNA <- NULL
  spe$nFeature_RNA <- NULL
  rownames(spe@colData) <- NULL
  spe$Location_Center_X <- NULL
  spe$Location_Center_Y <- NULL
  
  qs::qsave(spe, file = paste0("tables/clusters/interactions/IMC_spe_", exp, ".qs"))
  
  print(paste0("2/4 SPE object generated for: ", exp))
  
  ##### 3. calculate distances ######
  
  #spe <- qs::qread(paste0("tables/SIMPLI_analyses/clusters/IMC_spe_", exp, ".qs"))
  
  # from here
  # https://bodenmillergroup.github.io/imcRtools/reference/buildSpatialGraph.html
  
  spe <- buildSpatialGraph(spe, 
                           img_id = "sample_id", 
                           type = "knn", 
                           k = 20
  )
  
  # spe <- buildSpatialGraph(spe, 
  #                          img_id = "sample_id", 
  #                          type = "expansion", 
  #                          threshold = 20)
  # 
  # spe <- buildSpatialGraph(spe, 
  #                          img_id = "sample_id", 
  #                          type = "delaunay", 
  #                          max_dist = 20)
  # 
  # colPairNames(spe)
  
  # colPair(spe, "knn_interaction_graph")
  # colPair(spe, "expansion_interaction_graph")
  # colPair(spe, "delaunay_interaction_graph")
  
  qs::qsave(spe, file = paste0("tables/clusters/interactions/IMC_spe_", exp, ".qs"))
  
  print(paste0("3/4 knn distances calculated for: ", exp))
  
  ##### 4. perform spatial analyses ######
  library(scales)
  library(BiocParallel)
  register(SerialParam())
  # from here:
  # https://bodenmillergroup.github.io/IMCDataAnalysis/performing-spatial-analysis.html#interaction-analysis
  
  spe <- qs::qread(paste0("tables/clusters/interactions/IMC_spe_", exp, ".qs"))
  
  out <- testInteractions(spe, 
                          group_by = "sample_id",
                          label = "cluster", 
                          colPairName = "knn_interaction_graph",
                          BPPARAM = SerialParam(RNGseed = 123))
  
  #head(out)
  
  # prepare for plot
  out_plot <- out %>% as_tibble() %>%
    dplyr::group_by(from_label, to_label) %>%
    dplyr::summarize(sum_sigval = sum(sigval, na.rm = TRUE))
  
  # filter
  out_plot <- out_plot[out_plot$from_label %in% exclude,]
  out_plot <- out_plot[!(out_plot$to_label %in% exclude),]
  
  # rescale signal
  out_plot$sum_sigval <- rescale(out_plot$sum_sigval, to = c(-5, 5))
  
  #out_plot$facet <- ifelse(out_plot$from_label %in% c("12", "15", "24", "2", "3", "34", "35"), "glia", "neurons")
  
  # change cluster name
  out_plot$from_label <- paste0("cluster_" , out_plot$from_label)
  out_plot$to_label <- paste0("cluster_" , out_plot$to_label)
  
  # save
  write_csv(out_plot, paste0("tables/clusters/interactions/out_plot_interactions_", exp, ".csv"))
  
  clusters_labels <- read_csv("input/clusters_labels.csv")
  out_plot$from_label <- clusters_labels$full_label[match(out_plot$from_label, clusters_labels$cluster)]
  out_plot$to_label <- clusters_labels$full_label[match(out_plot$to_label, clusters_labels$cluster)]
  
  # order 
  order <- setdiff(clusters_oi, exclude)
  order <- paste0("cluster_", order)
  order <- clusters_labels$full_label[match(order, clusters_labels$cluster)]
  
  #out_plot$from_label <- factor(out_plot$from_label, levels = order)
  out_plot$to_label <- factor(out_plot$to_label, levels = order)
  
  # plot
  p <- ggplot(out_plot) +
    geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
    scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          panel.background = element_blank())
    #facet_grid(.~facet, space = "free", scales = "free_x")
  
  #p
  
  # save heatmap
  pdf(paste0("plots/clusters/interactions/interaction_", exp, "_knn.pdf"), 
      height = 4, 
      width = 3.8)
  print(p)
  dev.off()
  
  print(paste0("4/4 knn interactions calculated for: ", exp))
  
  print(paste0("FINSHED processing analysis: ", exp))
  
}

