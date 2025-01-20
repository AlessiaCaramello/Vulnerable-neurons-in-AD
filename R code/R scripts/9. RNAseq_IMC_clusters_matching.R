# matching IMC clusters to snRNASeq clusters + markers weighted (by intensity of marker expression)

# input 
# - IMC_markers_expression_in_RNAseq.csv (from snRNAseq)
# - clustered_cells.csv (input files)
# - clusters_labels.csv (metadata)
# - clusters_quartile_matching_score.csv (generated in 8.)


##### 0a. Load libraries ####

library(stringr) 
library(tidyverse)
library(plyr)
library(readr)
library(stringr)
library("rqdatatable")
library(dplyr)
#library(ggExtra)
library(gridExtra)
library(ggsignif)
#install.packages("ggpubr")
library(ggpubr)
library(rstatix)
library(plyr)
library(reshape2)
#install.packages("pheatmap")
library(pheatmap)
#install.packages("scales")                              
library("scales")
library(grid)
#install.packages("gridExtra")  
library(gridBase)
library(gridExtra)
library(ComplexHeatmap)
library('circlize')

setwd("~/UK Dementia Research Institute Dropbox/Alessia Caramello/Alessia Lab Dropbox/Projects/UKDRI project/Image analysis/HistoCat:SIMPLI analysis/051022 TREM2 cohort/R code")

##### 0b. pick method >>OPEN<< ####

# define range for rescaling markers expression in both IMC and RNAseq datasets
start <- as.numeric(-2) 
end <- as.numeric(5)

# define number of top expressed genes to take into consideration for calculating matching score
markers_n <- as.numeric(7)

##### 1. import snRNAseq file and find top genes expressed in each cluster #####

RNAseq_clusters <- read_csv("input/IMC_markers_expression_in_RNAseq.csv")

# remove APP, MAP2, RBFOX3, SYP (expressed in all cells) and all non neuronal markers
RNAseq_clusters <- RNAseq_clusters[!(RNAseq_clusters$features.plot %in% 
                                       #c("APP", "MAP2", "RBFOX3", "SYP", "AIF1", "CD68", "GFAP", "OLIG2", "PLP1", "S100B", "ADARB2", "NTNG1")
                                       c("APP", "RBFOX3", "SYP", "AIF1", "CD68", "GFAP", "OLIG2", "PLP1", "S100B", "ADARB2", "NTNG1")
),]

# re-scale avg.exp per marker in range -5 to +5
RNAseq_clusters$avg.exp.scaled2 <- 0

RNAmarkers <- unique(RNAseq_clusters$features.plot)

for (m in 1:length(RNAmarkers)){
  temp <- subset(RNAseq_clusters, RNAseq_clusters$features.plot == RNAmarkers[m])
  RNAseq_clusters$avg.exp.scaled2[which(RNAseq_clusters$features.plot == RNAmarkers[m])] <- rescale(temp$avg.exp, to = c(start, end))
}

rm(temp)

# pick highest between GAD1/2

temp_GAD <- RNAseq_clusters %>% 
  subset(features.plot %in% c("GAD1", "GAD2")) %>%
  group_by(id) %>%
  slice(which.max(avg.exp.scaled2)) %>%
  ungroup()

# remove rows with GAD1/2 and add back those generated above

RNAseq_clusters <- RNAseq_clusters %>% 
  filter(!grepl(c("GAD"), features.plot))


RNAseq_clusters <- rbind(RNAseq_clusters, temp_GAD)
rm(temp_GAD)

#  change GAD2 --> GAD1 in RNAseq clusters
RNAseq_clusters$features.plot[which(RNAseq_clusters$features.plot == "GAD2")] <- "GAD1"

# rename CALB2 to Calretinin
RNAseq_clusters$features.plot[which(RNAseq_clusters$features.plot == "CALB2")] <- "Calretinin"

# extract top 7 markers
require(data.table)
RNAseq_clusters <- data.table(RNAseq_clusters, key="avg.exp.scaled2")
RNAseq_clusters_top <- RNAseq_clusters[, tail(.SD, markers_n), by=id]

##### 2. import IMC file and find top genes expressed in each cluster #####

IMC_clusters <- read_csv("input/clustered_cells_original.csv")

# select only columns of interest
IMC_clusters_mean <- IMC_clusters[, colnames(IMC_clusters) %in% 
                                    c('res_0.9_ids', 
                                      colnames(IMC_clusters)[stringr::str_detect(colnames(IMC_clusters),"Intensity_MeanIntensity_")]
                                    )]

# rename markers columns
colnames(IMC_clusters_mean)[stringr::str_detect(colnames(IMC_clusters_mean),"Intensity_MeanIntensity_")] <-
  str_remove(colnames(IMC_clusters_mean)[stringr::str_detect(colnames(IMC_clusters_mean),"Intensity_MeanIntensity_")], "Intensity_MeanIntensity_")

# calculate average intensity by cluster
colnames(IMC_clusters_mean)[1] <- "clusters" # rename cluster column
IMC_clusters_mean <- aggregate(IMC_clusters_mean[, 2:length(IMC_clusters_mean)], # calculate mean
                               list(IMC_clusters_mean$clusters), mean)
colnames(IMC_clusters_mean)[1] <- "clusters" # rename cluster column
IMC_clusters_mean$clusters <- paste0("cluster_", IMC_clusters_mean$clusters) # rename clusters

# remove APP, MAP2, RBFOX3 (NeuN), SYP (Synapto) columns (expressed in all cells in snRNAseq) and non-neuronal
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("DNA1", colnames(IMC_clusters_mean)))]
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("APP", colnames(IMC_clusters_mean)))]
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("NeuN", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("Synapto", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("pTau", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("MAP2all", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("Ab", colnames(IMC_clusters_mean)))] 
#IMC_clusters_mean <- IMC_clusters_mean[, -(grep("MAP2", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("PLP1", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("OLIG2", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("GFAP", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("S100B", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("CD68", colnames(IMC_clusters_mean)))] 
IMC_clusters_mean <- IMC_clusters_mean[, -(grep("Iba1", colnames(IMC_clusters_mean)))]

# turn table to long
IMC_clusters_long <- IMC_clusters_mean %>% pivot_longer(cols=!c('clusters'),
                                                        names_to='marker',
                                                        values_to='avg.exp')

# rename LMO3 to LMO4
IMC_clusters_long$marker[which(IMC_clusters_long$marker == "LMO3")] <- "LMO4"

# remove non-neuronal clusters
IMC_clusters_long <- IMC_clusters_long[!(IMC_clusters_long$clusters %in% 
                                           c("cluster_2", "cluster_3", "cluster_34", "cluster_35", "cluster_0", 
                                             "cluster_30", "cluster_15", "cluster_24", "cluster_12", "cluster_4", 
                                             "cluster_6", "cluster_20")
),]

# re-scale avg.exp per marker in range defined above
IMC_clusters_long$avg.exp.scaled2 <- 0

IMCmarkers <- unique(IMC_clusters_long$marker)

for (m in 1:length(IMCmarkers)){
  temp <- subset(IMC_clusters_long, IMC_clusters_long$marker == IMCmarkers[m])
  IMC_clusters_long$avg.exp.scaled2[which(IMC_clusters_long$marker == IMCmarkers[m])] <- rescale(temp$avg.exp, to = c(start, end))
}

# adjust outlier in CCK expression
IMC_clusters_long$avg.exp.scaled2[which(IMC_clusters_long$marker == "CCK" & IMC_clusters_long$avg.exp.scaled2 == 5)] <- 0.1

# extract top 7 markers
require(data.table)
IMC_clusters_top <- data.table(IMC_clusters_long, key="avg.exp.scaled2")
IMC_clusters_top <- IMC_clusters_top[, tail(.SD, markers_n), by=clusters]

##### 3. calculate score (markers only) #####

# create datatable for storing info
markers_matrix <- data.frame(matrix(ncol = 7, 
                                    nrow = (length(unique(IMC_clusters_top$clusters))*length(unique(RNAseq_clusters_top$id)))))

colnames(markers_matrix) <- c('IMC_cluster', 'RNAseq_cluster', 'common_markers', 'score_IMC_w',
                              'score_RNAseq_w', 'shared_top', 'combined_marker_score')

# calculate correlation
for (i in 1:length(unique(IMC_clusters_top$clusters))){
  #i=26
  # subset IMC clusters
  IMCclust <- subset(IMC_clusters_top, IMC_clusters_top$clusters == unique(IMC_clusters_top$clusters)[i])
  IMCmark <- IMCclust$marker
  
  for (r in 1:length(unique(RNAseq_clusters_top$id))){
    #r=2
    # subset RNAseq clusters
    RNAclust <- subset(RNAseq_clusters_top, RNAseq_clusters_top$id == unique(RNAseq_clusters_top$id)[r])
    RNAmark <- RNAclust$features.plot
    
    # find matches
    common <- intersect(IMCmark, RNAmark)
    #uncommon <- dplyr::setdiff(IMCmark, RNAmark)
    
    # starting score
    scoreIMC <- 0
    scoreRNA <- 0
    
    # calculate slot
    slot <- (length(unique(RNAseq_clusters_top$id))*(i-1))+r
    
    # calculate RNA and IMC scores
    if (!is_empty(common)){
      for (m in 1:length(common)){
        #m=2
        # calculate weighted IMC score
        scoreIMC <<- scoreIMC + IMCclust$avg.exp.scaled2[grep(common[m], IMCclust$marker)]
        # calculate weighted RNAseq score
        scoreRNA <<- scoreRNA + RNAclust$avg.exp.scaled2[grep(common[m], RNAclust$features.plot)]
      }
    }
    
    # check if top marker is the same
    top_RNA <- RNAclust[avg.exp.scaled2 == max(avg.exp.scaled2), features.plot]
    top_IMC <- IMCclust[avg.exp.scaled2 == max(avg.exp.scaled2), marker]
    
    # add extra points if top marker is shared
    top_shared <- (top_RNA %in% top_IMC)
    
    if (length(top_shared[top_shared == TRUE]) > 0){
      markers_matrix$shared_top[slot] <- as.numeric(2)
    } else {
      markers_matrix$shared_top[slot] <- as.numeric(1)
    }
    
    # fill up matrix
    markers_matrix$IMC_cluster[slot] <- unique(IMC_clusters_top$clusters)[i]
    markers_matrix$RNAseq_cluster[slot] <- unique(RNAseq_clusters_top$id)[r]
    markers_matrix$common_markers[slot] <- paste0(common, collapse = ", ")
    markers_matrix$score_IMC_w[slot] <- scoreIMC
    markers_matrix$score_RNAseq_w[slot] <- scoreRNA
    markers_matrix$combined_marker_score[slot] <- scoreIMC*scoreRNA*markers_matrix$shared_top[slot]
    
    # change to minus sign if both score are negative
    if ((scoreIMC < 0) & (scoreRNA < 0)){
      markers_matrix$combined_marker_score[slot] <- (markers_matrix$combined_marker_score[slot])*-1
    }
    
  }
  
}

##### 4. import quartile scores ####

quartile_scores <- read_csv("tables/RNAseq/clusters_matching/clusters_quartile_matching_score.csv")
markers_scores <- markers_matrix

quartile_scores[,1] <- NULL

colnames(markers_scores)[which(colnames(markers_scores) == "RNAseq_cluster")] <- "RNA_cluster"

##### 5. calculate combined score (markers + quartiles) ####

# add quartile scores to markers score
combined_scores <- left_join(markers_scores, quartile_scores, by = c("RNA_cluster", "IMC_cluster"))
colnames(combined_scores)[which(colnames(combined_scores) == "score")] <- "quant_score"

# calculate combined score (mult + top marker) (S)
combined_scores$combined_marker_quant_score <- combined_scores$combined_marker_score * combined_scores$quant_score

# save 
write_csv(combined_scores, paste0("tables/RNAseq/clusters_matching/clusters_matches_combined_score.csv"))

##### 6. make heatmap #####

# select only one type of score 
markers_matrix_plot <- combined_scores[, names(combined_scores) %in% c("IMC_cluster", "RNA_cluster", "combined_marker_quant_score")]

# turn to wide table 
markers_matrix_plot <- markers_matrix_plot %>% pivot_wider(names_from = RNA_cluster, values_from = "combined_marker_quant_score")

# define order of clusters (rows)
clusters_order <- c("cluster_5",
                    "cluster_8",
                    "cluster_10",
                    "cluster_11",
                    "cluster_13",
                    "cluster_16",
                    "cluster_18",
                    "cluster_23",
                    "cluster_27",
                    "cluster_9",
                    "cluster_14",
                    "cluster_19",
                    "cluster_21",
                    "cluster_22",
                    "cluster_25",
                    "cluster_26",
                    "cluster_28",
                    "cluster_29",
                    "cluster_31",
                    "cluster_32",
                    "cluster_33",
                    "cluster_1",
                    "cluster_7",
                    "cluster_17")

markers_matrix_plot <- markers_matrix_plot[match(clusters_order, 
                                                 markers_matrix_plot$IMC_cluster),] # reorder dataframe
# turn column to rownames
markers_matrix_plot <- column_to_rownames(markers_matrix_plot, var = "IMC_cluster")

# define order of markers (columns)
markers_order <- c("Exc-L2-3-LINC00507-RPL9P17", "Exc-L3-5-FEZF2-ONECUT1", "Exc-L3-5-RORB-HSPB3", "Exc-L3-LINC00507-CTXN3", 
                   "Exc-L3-RORB-CARTPT", "Exc-L4-6-RORB-LCN15", "Exc-L5-6-FEZF2-RSAD2", "Exc-L5-6-THEMIS-GPR21",
                   "Exc-L5-RORB-LINC01202", "Exc-L6-FEZF2", "Exc-L6-THEMIS-LINC00343", "Inh-L1-4-VIP-CHRNA2", "Inh-L1-5-VIP-KCNJ2", 
                   "Inh-L3-4-PVALB-HOMER3", "Inh-L6-SST-NPY", "Inh-LAMP5"
)

markers_matrix_plot <- markers_matrix_plot[,match(markers_order, colnames(markers_matrix_plot))] # reorder dataframe

# read in labels and modify 
labels <- read.csv("input/clusters_labels_original.csv")
cluster_labels <- labels[, c(1,4)] # for annotation_row

# remove non-neuronal clusters
cluster_labels <- cluster_labels[!(cluster_labels$cluster %in% 
                                     c("cluster_2", "cluster_3", "cluster_4", "cluster_34", "cluster_35", "cluster_0", 
                                       "cluster_30", "cluster_15", "cluster_24", "cluster_12", "cluster_6", "cluster_20")
),]

# turn column to rownames
rownames(cluster_labels) <- cluster_labels$cluster
cluster_labels[, 1] <- NULL
colnames(cluster_labels) <- "assigned_cluster"

# make as factor to define order
cluster_labels$assigned_cluster <- factor(cluster_labels$assigned_cluster,
                                          levels = c("ex_neurons", "inh_neurons", "pan_neurons"))

# define colors for legend
colors = list(
  assigned_cluster = c(
    ex_neurons = "#f66d9b", 
    inh_neurons = "#9561e2",
    pan_neurons = "#3490dc"))

# make into matrix for heatmap plotting
m <- as.matrix(markers_matrix_plot)

# Plot
plot <- pheatmap(m,
                 scale = "row",
                 #color = colorRampPalette(c("navy", "white", "red"))(50), 
                 #annotation_col = markers_labels, 
                 annotation_row = cluster_labels,
                 annotation_colors = colors,
                 #display_numbers = TRUE,
                 #number_color = "black", 
                 #fontsize_number = 4, 
                 border_color = NA,
                 cluster_rows = FALSE, # remove dendograms
                 cluster_cols = FALSE, # remove dendograms
                 #main = "markers expression in clusters", 
                 gaps_col = c(11),
                 gaps_row = c(9, 21)
) 
plot

# save heatmap
pdf(paste0("plots/RNAseq/clusters_matching/combined_scores_matching_IMC_and_RNA_clusters.pdf"), 
    height = 10, 
    width = 8)
print(plot)
dev.off()

##### 7. select top match (IMC to RNA) #####

# select only one type of score 
markers_matrix_table <- combined_scores[, names(combined_scores) %in% c("IMC_cluster", "RNA_cluster", "combined_marker_quant_score")]

# get top matching RNA cluster for each IMC cluster
top_match <- markers_matrix_table %>% group_by(IMC_cluster) %>% slice_max(combined_marker_quant_score, n = 1)

colnames(top_match)[3] <- "score"

# add column for common markers
top_match$common_markers <- 0

# add common markers
for (i in 1:nrow(top_match)){
  top_match$common_markers[i] <- combined_scores$common_markers[which(combined_scores$IMC_cluster == top_match$IMC_cluster[i] 
                                                                     & combined_scores$RNA_cluster == top_match$RNA_cluster[i])]
}

# read in labels and modify 
labels <- read.csv("input/clusters_labels_original.csv")
cluster_labels <- labels[, c(1,4)] # for annotation_row

# remove non-neuronal clusters
cluster_labels <- cluster_labels[!(cluster_labels$cluster %in% 
                                     c("cluster_2", "cluster_3", "cluster_4", "cluster_34", "cluster_35", "cluster_0", 
                                       "cluster_30", "cluster_15", "cluster_24", "cluster_12", "cluster_6", "cluster_20")
),]

# add column for IMC clusters labels
top_match$IMC_clusters_label <- 0

# add clusters labels
for (i in 1:nrow(top_match)){
  top_match$IMC_clusters_label[i] <- labels$label[which(labels$cluster == top_match$IMC_cluster[i])]
}

# save
write_csv(top_match, "tables/RNAseq/clusters_matching/top_match_between_IMC_to_RNA_clusters.csv")

##### 8. make chord plot #####

# Import top_match
#top_match <- read_csv("tables/RNAseq/clusters_matching/top_match_between_IMC_to_RNA_clusters.csv")

# select columns of interest
top1_circle <- top_match %>% select(c("IMC_cluster", "RNA_cluster", "score"))

# normalise score
top1_circle_norm <- top1_circle
top1_circle_norm$score <- scales::rescale(top1_circle$score, to = c(1, 10))

# define order of markers (columns)
markers_order <- c("Inh-L1-4-VIP-CHRNA2", "Inh-L1-5-VIP-KCNJ2", "Inh-L3-4-PVALB-HOMER3", "Inh-L6-SST-NPY", "Inh-LAMP5", 
                   "Exc-L2-3-LINC00507-RPL9P17", "Exc-L3-5-FEZF2-ONECUT1", "Exc-L3-5-RORB-HSPB3", "Exc-L3-LINC00507-CTXN3", 
                   "Exc-L3-RORB-CARTPT", "Exc-L4-6-RORB-LCN15", "Exc-L5-6-FEZF2-RSAD2", "Exc-L5-6-THEMIS-GPR21",
                   "Exc-L5-RORB-LINC01202", "Exc-L6-FEZF2", "Exc-L6-THEMIS-LINC00343"
)

# define order of clusters (rows)
clusters_order <- c("cluster_5",
                    "cluster_8",
                    "cluster_10",
                    "cluster_11",
                    "cluster_13",
                    "cluster_16",
                    "cluster_18",
                    "cluster_23",
                    "cluster_27",
                    "cluster_1",
                    "cluster_7",
                    "cluster_17",
                    "cluster_9",
                    "cluster_14",
                    "cluster_19",
                    "cluster_21",
                    "cluster_22",
                    "cluster_25",
                    "cluster_26",
                    "cluster_28",
                    "cluster_29",
                    "cluster_31",
                    "cluster_32",
                    "cluster_33"
)


# turn into wide table
#m <- top1_circle_norm %>% pivot_wider(names_from = RNAseq_cluster, values_from = score, values_fill = 0)

# add missing columns
# diff <- setdiff(markers_order, colnames(m))
# for (i in length(diff)){
#   m[, diff] <- NA
# }

# keep 
same <- intersect(markers_order, colnames(m))
markers_order <- markers_order[markers_order %in% same]

# order columns
#m <- m[, c("IMC_cluster", markers_order)] # reorder dataframe

# order rows
#m <- m[order(factor(m$IMC_cluster, levels = clusters_order)), ] # reorder dataframe

order_plot <- c(clusters_order, markers_order)

# Define colors
colors <- c(cluster_5 = "#f66d9b", 
            cluster_8 = "#f66d9b", 
            cluster_10 = "#f66d9b", 
            cluster_11 = "#f66d9b", 
            cluster_13 = "#f66d9b",
            cluster_16 = "#f66d9b",
            cluster_18 = "#f66d9b",
            cluster_23 = "#f66d9b",
            cluster_27 = "#f66d9b",
            cluster_1 = "#3490dc",
            cluster_7 = "#3490dc",
            cluster_17 = "#3490dc",
            cluster_9 = "#9561e2",
            cluster_14 = "#9561e2",
            cluster_19 = "#9561e2",
            cluster_21 = "#9561e2",
            cluster_22 = "#9561e2",
            cluster_25 = "#9561e2",
            cluster_26 = "#9561e2",
            cluster_28 = "#9561e2",
            cluster_29 = "#9561e2",
            cluster_31 = "#9561e2",
            cluster_32 = "#9561e2",
            cluster_33 = "#9561e2", 
            "Inh-L1-4-VIP-CHRNA2" = "#9561e2",
            "Inh-L3-4-PVALB-HOMER3" = "#9561e2",
            "Inh-L6-SST-NPY" = "#9561e2",
            "Inh-LAMP5" = "#9561e2",
            "Exc-L3-5-FEZF2-ONECUT1" = "#f66d9b",
            "Exc-L2-3-LINC00507-RPL9P17" = "#f66d9b",
            "Exc-L3-RORB-CARTPT" = "#f66d9b",
            "Exc-L3-5-RORB-HSPB3" = "#f66d9b",
            "Exc-L4-6-RORB-LCN15" = "#f66d9b",
            "Exc-L5-6-FEZF2-RSAD2" = "#f66d9b",
            "Exc-L5-6-THEMIS-GPR21" = "#f66d9b",
            "Exc-L5-RORB-LINC01202" = "#f66d9b",
            "Exc-L6-FEZF2" = "#f66d9b",
            "Exc-L6-THEMIS-LINC00343" = "#f66d9b"
)

pdf(paste0("plots/RNAseq/clusters_matching/clusters_matching_chord_plot.pdf"), 
    width=8, 
    height=8)

# edit initialising parameters
circlize::circos.par(canvas.ylim=c(-1.5,1.5), # edit  canvas size 
           #gap.after = gaps, # adjust gaps between regions
           track.margin = c(0.01, 0), # adjust bottom and top margin
           # track.margin = c(0.01, 0.1)
           track.height = 0.03)

# plot
chordDiagram(top1_circle_norm, 
             order = order_plot, 
             grid.col = colors,
             annotationTrack = "grid", 
             preAllocateTracks = 2) # <-- define number of tracks (sectors highlight and labels)

# highlight sectors
highlight.sector(top1_circle_norm$IMC_cluster, track.index = 2, col = "#DD726F")
highlight.sector(top1_circle_norm$RNA_cluster, track.index = 2, col = "#37ABB7")

# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

# Add legend
lgd_list = Legend(
  labels = c("IMC", "RNA"),
  title = "clusters_from",
  #background = c("#DD726F","#37ABB7"))
  legend_gp = gpar(fill = c("#DD726F","#37ABB7")))

lgd_col = Legend(
  labels = c("ex_neurons", "pan_neurons", "inh_neurons"),                 
  title = "assigned_clusters",    
  #background = unique(colors))
  legend_gp = gpar(fill = unique(colors)))

lgd_list = packLegend(lgd_list, lgd_col, direction = "vertical")
circle_size = unit(1, "snpc") 
# snpc unit gives you a square region

draw(lgd_list, 
     x = unit(20, "mm"), 
     y = unit(25, "mm"),
     #x = circle_size, 
     just = "left")

# Restart circular layout parameters
circos.clear()

dev.off()


