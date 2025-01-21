# UMAP plotting for samples, markers and clusters 
# adapted from SIMPLI plotting functions

# Input data:
# - clustered_cells.csv (SIMPLI output)

####### SETUP from SIMPLI #######
library(data.table)
library(uwot)
library(ggplot2)
#install.packages("qs", type = "binary")
library(qs)

#source("scripts/Plot_Functions.R")

# set working directory
setwd("~/R code")

#define variables
cell_file_name <- "input/clustered_cells.csv"                        # $clustered_cell_file \\                                                                                         
#sample_metadata_file_name <- "input/sample_metadata_TREM2.csv"       # $sample_metadata_file \\
#cell_population <- "all_cells"                                       # $cell_type of clustered_cell \\                                                                                        
markers <- strsplit("Intensity_MeanIntensity_CCK@Intensity_MeanIntensity_GAD1@Intensity_MeanIntensity_CALB1@Intensity_MeanIntensity_VIP@Intensity_MeanIntensity_PVALB@Intensity_MeanIntensity_SST@Intensity_MeanIntensity_NPY@Intensity_MeanIntensity_ADARB1@Intensity_MeanIntensity_LHX6@Intensity_MeanIntensity_Calretinin@Intensity_MeanIntensity_CUX2@Intensity_MeanIntensity_LMO3@Intensity_MeanIntensity_PCP4@Intensity_MeanIntensity_GPC5@Intensity_MeanIntensity_FOXP2@Intensity_MeanIntensity_RORB@Intensity_MeanIntensity_MAP2@Intensity_MeanIntensity_MAP2all@Intensity_MeanIntensity_NeuN@Intensity_MeanIntensity_OLIG2@Intensity_MeanIntensity_GFAP@Intensity_MeanIntensity_S100B@Intensity_MeanIntensity_Iba1@Intensity_MeanIntensity_CD68@Intensity_MeanIntensity_Synapto@Intensity_MeanIntensity_NTNG2@Intensity_MeanIntensity_pTau@Intensity_MeanIntensity_Ab@Intensity_MeanIntensity_APP", "@")[[1]]     # $markers \\                                                        
resolutions <- as.numeric("0.9")                                # $resolutions \\                                                   
high_color <- "#0000FF"                                         # $params.high_color \\                                                       
mid_color <- "#FFFFFF"                                          # $params.mid_color \\                                                           
low_color <- "#DCDCDC"                                          # $params.low_color \\                                                        
output_folder <- "plots/clusters/UMAPs"  
dir.create(output_folder)
dir.create("tables/clusters/UMAPs")

resolutions <- paste0("res_", resolutions, "_ids")

######## Cell data #######
Cells <- fread(cell_file_name)
# Samples <- fread(sample_metadata_file_name)
# Samples[is.na(comparison), comparison := "NA" ]
# suppressWarnings(Cells[, color := NULL])
# suppressWarnings(Cells[, comparison := NULL])

# Cells <- merge(Cells, Samples, by.x = "Metadata_sample_name", by.y = "sample_name")
# Cells <- Cells[comparison != "NA" & cell_type == cell_population]

######## Heatmaps #######
# heatmaps <- lapply(resolutions, function(res){
#   heatmapper(Cells[, c(res, markers), with = F], res_column = res, cols = markers, high_color,
#              mid_color, low_color)
# })
# names(heatmaps) <- resolutions

################## UMAPS #####################
set.seed(666)

checkna_warn <- function(X){  ### Override the uwot:::checkna which stops with an error!
  if (!is.null(X) && any(is.na(X))) {
    warning("Missing values found in 'X'")
  }
}	
rlang::env_unlock(env = asNamespace('uwot'))
rlang::env_binding_unlock(env = asNamespace('uwot'))
assign('checkna', checkna_warn, envir = asNamespace('uwot'))
rlang::env_binding_lock(env = asNamespace('uwot'))
rlang::env_lock(asNamespace('uwot'))

UMAPS <- lapply(resolutions, function(res){
  umap_coords <- umap(Cells, n_neighbors = 40, min_dist = 0.9, learning_rate = 0.5, init = "random",
                      metric = list("euclidean" = markers, "categorical" = res), 
                      n_sgd_threads = 1, n_threads = 1)
  umap_table <- cbind(umap_coords, Cells[, c("Metadata_sample_name", markers, res), with = F])
  setnames(umap_table, c("V1", "V2"), c("umap_x", "umap_y"))
  return(umap_table)
})
names(UMAPS) <- resolutions

# save
#qsave(UMAPS, "tables/clusters/UMAPs/UMAPs.qs")
saveRDS(UMAPS, "tables/clusters/UMAPs/UMAPs_new.rds")

###################### UMAP Plots by Cluster ##########################

label_dot_plotter <- function(data, x_name, y_name, label_name, plot_title = NULL)
{
  plot_data <- copy(data)
  data[[label_name]] <- as.character(data[[label_name]])	
  my_plot <- ggplot(data) +
    geom_point(mapping = aes_string(x = x_name, y = y_name, color = label_name), size = 0.5) +
    labs(x = x_name, y = y_name) +
    theme_bw(base_size = 15, base_family = "sans") +
    labs(title = plot_title) +
    theme(plot.title = element_blank(), 
          legend.position = "right", 
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black", size = 0.75),
          strip.background = element_rect(colour = NA, fill = NA),
          legend.key.size = unit(1, 'cm'))
  return(my_plot)
}

umap_cluster_plots <- lapply(resolutions, function(res){
  label_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", res), with = F], "umap_x", "umap_y", res, "Cluster")
})
names(umap_cluster_plots) <- resolutions

png(file=paste0("plots/clusters/UMAPs/UMAP_clusters_", resolutions, ".png"), # file with smaller size
    width=3200, height=2500, res=400)
umap_cluster_plots[[resolutions]]
dev.off()

pdf(paste0("plots/clusters/UMAPs/UMAP_clusters_", resolutions, ".pdf"))
umap_cluster_plots[[resolutions]]
dev.off()


###################### UMAP Plots by Sample ##########################
umap_sample_plots <- lapply(resolutions, function(res){
  label_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", "Metadata_sample_name"), with = F], "umap_x", "umap_y",
                    "Metadata_sample_name", "Sample")
})
names(umap_sample_plots) <- resolutions

png(file=paste0("plots/clusters/UMAPs/UMAP_samples_", resolutions, ".png"), # file with smaller size
    width=3500, height=2500, res=400)
umap_sample_plots[[resolutions]]
dev.off()

pdf(paste0("plots/clusters/UMAPs/UMAP_samples_", resolutions, ".pdf"))
umap_sample_plots[[resolutions]]
dev.off()

###################### UMAP Plots by Marker ##########################

gradient_dot_plotter <- function(data, x_name, y_name, marker, highcol, lowcol)
{
  plot_data <- copy(data)
  setnames(plot_data, c(x_name, y_name, marker), c("x_name", "y_name", "marker"))
  plot_data[, marker := marker / max(marker)]
  my_plot <- ggplot(plot_data) +
    geom_point(mapping = aes(x = x_name, y = y_name, color = marker), size = 0.5) +
    scale_color_gradient(low = lowcol, high = highcol) +
    theme_bw(base_size = 15, base_family = "sans") +
    labs(x = x_name, y = y_name, title = gsub("Intensity_MeanIntensity_", "", marker)) +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.title = element_blank(),
          legend.position = "right", 
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black", size = 0.75),
          strip.background = element_rect(colour = NA, fill = NA),
          legend.key.size = unit(1, 'cm'))
  return(my_plot)
}

umap_marker_plots <- lapply(resolutions, function(res){
  plots <- lapply(markers, function(marker){
    gradient_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", res, marker), with = F],
                         "umap_x", "umap_y", marker, high_color, low_color)  
  })
  names(plots) <- markers
  return(plots)
})
names(umap_marker_plots) <- resolutions


for (i in 1:length(names(umap_marker_plots[[resolutions]]))){
  plot_name <- gsub("Intensity_MeanIntensity_", "", names(umap_marker_plots[[resolutions]])[i])
  png(file=paste0("plots/clusters/UMAPs/UMAP_", plot_name, "_", resolutions, ".png"),
      width=2900, height=2600, res=400)
  print(umap_marker_plots[[resolutions]][i])
  dev.off()
  print(paste0(plot_name, " plot printed"))
}
