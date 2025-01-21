# Flip images to correct orientation by changing XY coordinates 

# input data: 
#  - clustered_cells.csv (SIMPLI output)
#  - sample_rotation.csv (metadata - degrees of rotation for each image)
  
##### 0. Load libraries #####

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
#install.packages('openxlsx')
library(openxlsx)
#install.packages("writexl")
library(writexl)
#install.packages("readxl)
library(readxl)
#install.packages("openxlsx") 
library(openxlsx)

##### 1. Set working directory and import files #####

# set working directory
setwd("~/R code")

# load input files
clustered_cells <- read.csv("input/clustered_cells.csv")

sample_names <- unique(clustered_cells$Metadata_sample_name)

##### 2. Change X and Y for rotating images - LONG #####

# create sample_rotation file (already done)
#sample_rotation <- data.frame(sample_name = c(sample_names), rotation = 0) # only clockwise rotation - insert the degree of rotation (90-180-270) to end up with layer 1 on top
#write.xlsx(x = sample_rotation, 'input/sample_rotation.xlsx')

#insert rotational angle and load file again
sample_rotation <- read.xlsx('input/sample_rotation.xlsx')

# add new columns for the new X and Y coordinates
colnames(clustered_cells)[which(colnames(clustered_cells) == 'Location_Center_X')] <- 'X'
colnames(clustered_cells)[which(colnames(clustered_cells) == 'Location_Center_Y')] <- 'Y'

clustered_cells <- add_column(clustered_cells, "new_X" = 0 , .after = "X")
clustered_cells <- add_column(clustered_cells, "new_Y" = 0 , .after = "Y")

# add new column for degrees of rotation 
clustered_cells <- add_column(clustered_cells, "rotation" = 0)

# fill in with degrees of rotation
for (i in 1:nrow(clustered_cells)){
  degrees <- grep(clustered_cells$Metadata_sample_name[i], 
                  sample_rotation$sample_name)
  degrees <- as.numeric(sample_rotation$rotation[degrees])
  clustered_cells$rotation[i] <- degrees
}

# change X and Y - this may take a long time!!
for (i in 1:nrow(clustered_cells)){
  
  print(paste0("number of rows done:", i, "/", nrow(clustered_cells)))
  
  # extract rows of sample replicate to process
  sample_rows <- clustered_cells %>% filter(
    clustered_cells$Metadata_sample_name == 
      clustered_cells$Metadata_sample_name[i]) 
  
  if (clustered_cells$rotation[i] == 0) {
    
    # new Y = old Y
    clustered_cells$new_Y[i] <- clustered_cells$Y[i]
    
    # new X = old X
    clustered_cells$new_X[i] <- clustered_cells$X[i]
    
  }
  
  if (clustered_cells$rotation[i] == 90) {
    
    # new Y = old X
    clustered_cells$new_Y[i] <- clustered_cells$X[i]
    
    # new X = max old Y - old Y
    clustered_cells$new_X[i] <- max(sample_rows$Y) - clustered_cells$Y[i]
    
  }
  if (clustered_cells$rotation[i] == 180) {
    
    # new X = max old X - old X
    clustered_cells$new_X[i] <- max(sample_rows$X) - clustered_cells$X[i]
    
    # new Y = max old Y - old Y
    clustered_cells$new_Y[i] <- max(sample_rows$Y) - clustered_cells$Y[i]
    
  }
  if (clustered_cells$rotation[i] == 270) {
    
    # new X = old Y
    clustered_cells$new_X[i] <- clustered_cells$Y[i]
    
    # new Y = max old X - old X
    clustered_cells$new_Y[i] <- max(sample_rows$X) - clustered_cells$X[i]
  }
}

# remove columns of old X and Y
clustered_cells <- subset(clustered_cells, select = -c(X, Y, rotation))

# rename X and Y columns
colnames(clustered_cells)[which(names(clustered_cells) == "new_X" )] <- "X"
colnames(clustered_cells)[which(names(clustered_cells) == "new_Y" )] <- "Y"

# save new table
write.csv(clustered_cells, "input/clustered_cells.csv")

