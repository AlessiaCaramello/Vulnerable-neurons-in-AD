## Calculate IMC to RNAseq clusters matching score based on size (cell number)

# input files:
# - RNAseq_sample_metadata.tsv (metadata)
# - RNAseq_cluster_size.tsv (RNAseq data)
# - cells_per_cluster_wide.csv (generated in 3. cells number or density per cluster.R)

##### 0. Load libraries ####

library(stringr) 
library(tidyverse)
library(plyr)
library(readr)
library("rqdatatable")
library(dplyr)
#library(ggExtra)
library(gridExtra)
library(ggsignif)
library(ggpubr)
library(rstatix)
library(reshape2)
library(pheatmap)
library("scales")
#install.packages("ggthemes")
library(ggthemes)
library(viridis)   

setwd("~/Dropbox (UK Dementia Research Institute)/Alessia Lab Dropbox/Projects/UKDRI project/Image analysis/HistoCat:SIMPLI analysis/051022 TREM2 cohort/R code")

dir.create("tables/RNAseq")
dir.create("tables/RNAseq/clusters_matching")

dir.create("plots/RNAseq")
dir.create("plots/RNAseq/clusters_matching")

##### 1. Load input datasets ####

RNA_metadata <- read.table(file = 'input/RNAseq_sample_metadata.tsv', header = TRUE)
RNA_size <- read.table(file = 'input/RNAseq_cluster_size.tsv', sep = '\t', header = TRUE)
colnames(RNA_size) <- gsub("[.]", "-", colnames(RNA_size))

IMC_size <- read.csv("tables/clusters/cells_per_cluster_wide.csv")

##### 2. calculate total cells in RNA and IMC clusters #####

############### RNA 
# add column with group (diagnosis + trem2var)
RNA_metadata <- RNA_metadata %>% mutate(group =
                                          case_when(diagnosis == "AD" & trem2_group == "CV" ~ "AlzCV", 
                                                    diagnosis == "AD" & trem2_group == "TREM2var" ~ "AlzTREM2", 
                                                    diagnosis == "Control" & trem2_group == "CV" ~ "CtrlCV", 
                                                    diagnosis == "Control" & trem2_group == "TREM2var" ~ "CtrlTREM2", )
)

# add group column to RNA_size dataset
for (i in 1:nrow(RNA_size)){
  RNA_size$group[i] <- 
    RNA_metadata$group[grep(RNA_size$manifest[i], RNA_metadata$manifest)] 
}   

# turn into long table
RNA_size <- RNA_size %>% pivot_longer(cols = -c("manifest", "group"), 
                                      names_to ='cluster',
                                      values_to ='cell_number')
RNA_size_tot <- RNA_size %>% 
  group_by(cluster) %>% 
  dplyr::summarise(tot_RNA = sum(cell_number))

colnames(RNA_size_tot)[which(colnames(RNA_size_tot) == "cluster")] <- "RNA_cluster"

############### IMC
IMC_size_tot <- IMC_size %>% 
  select(starts_with("cluster_")) %>%
  colSums() %>%
  as.data.frame() 

names(IMC_size_tot) <- "tot_IMC"
IMC_size_tot$IMC_cluster <- rownames(IMC_size_tot)
rownames(IMC_size_tot) <- NULL

# remove non neuronal cells
remove <- paste0("cluster_", c("2", "3", "34", "35", "12", "15", "24", "30", "0", "4"))
IMC_size_tot <- IMC_size_tot[!(IMC_size_tot$IMC_cluster %in% remove), ]

tot_IMC_neuro <- sum(IMC_size_tot$tot_IMC)

##### 3. calculate quartile #####

IMC_quant <- quantile(IMC_size_tot$tot_IMC)
RNA_quant <- quantile(RNA_size_tot$tot_RNA)

# add quartile
# IMC
IMC_size_tot$quartile[IMC_size_tot$tot_IMC < IMC_quant[2]] <- as.numeric(1)
IMC_size_tot$quartile[IMC_size_tot$tot_IMC <= IMC_quant[3] & IMC_size_tot$tot_IMC > IMC_quant[2]] <- as.numeric(2)
IMC_size_tot$quartile[IMC_size_tot$tot_IMC <= IMC_quant[4] & IMC_size_tot$tot_IMC > IMC_quant[3]] <- as.numeric(3)
IMC_size_tot$quartile[IMC_size_tot$tot_IMC > IMC_quant[4]] <- as.numeric(4)
# RNA
RNA_size_tot$quartile[RNA_size_tot$tot_RNA < RNA_quant[2]] <- as.numeric(1)
RNA_size_tot$quartile[RNA_size_tot$tot_RNA <= RNA_quant[3] & RNA_size_tot$tot_RNA > RNA_quant[2]] <- as.numeric(2)
RNA_size_tot$quartile[RNA_size_tot$tot_RNA <= RNA_quant[4] & RNA_size_tot$tot_RNA > RNA_quant[3]] <- as.numeric(3)
RNA_size_tot$quartile[RNA_size_tot$tot_RNA > RNA_quant[4]] <- as.numeric(4)

##### 4. calculate similarity based on quartile ######

quartile_matches <- data.frame(matrix(ncol = 5))
colnames(quartile_matches) <- c("RNA_cluster", "RNA_quant", "IMC_cluster", "IMC_quant", "score")

temp <- data.frame(matrix(nrow = nrow(RNA_size_tot)))
colnames(temp) <- "RNA_cluster"
temp$RNA_cluster <- RNA_size_tot$RNA_cluster
temp$RNA_quant <- RNA_size_tot$quartile

for (clu in IMC_size_tot$IMC_cluster){
  #clu = "cluster_1"
  
  # add IMC quartile info
  temp2 <- temp
  temp2$IMC_cluster <- clu
  temp2$IMC_quant <- IMC_size_tot$quartile[IMC_size_tot$IMC_cluster == clu]
  
  for (row in 1:nrow(temp2)){
    
    if (temp2$RNA_quant[row] == temp2$IMC_quant[row]){ # identical (x2)
      temp2$score[row] <- as.numeric(2)
    } else if (temp2$RNA_quant[row] == (temp2$IMC_quant[row]+1) | temp2$RNA_quant[row] == (temp2$IMC_quant[row]-1)){ # similar ±1 (x1.5)
      temp2$score[row] <- as.numeric(1.75)
    } else if (temp2$RNA_quant[row] == (temp2$IMC_quant[row]+2) | temp2$RNA_quant[row] == (temp2$IMC_quant[row]-2)){ # far ±2 (x1)
      temp2$score[row] <- as.numeric(1)
    } else if (temp2$RNA_quant[row] == (temp2$IMC_quant[row]+3) | temp2$RNA_quant[row] == (temp2$IMC_quant[row]-3)){ # opposite (x0.5)
      temp2$score[row] <- as.numeric(0.5)
    }
    
  }
  
  quartile_matches <- rbind(quartile_matches, temp2)
  
}

if (is.na(quartile_matches[1,1])){
  quartile_matches <- quartile_matches[-1,]
}

rm(temp, temp2)

# save
write.csv(quartile_matches, "tables/RNAseq/clusters_matching/clusters_quartile_matching_score.csv")


