# Extract pixel coordinates from masks PNG files

# input files:
# - PNG files of plaques and rings masks from ImageJ macro


##### 0. Install packages and set working directory #####

library(dplyr)
library('png')

# set working directory
setwd("~/UK Dementia Research Institute Dropbox/Alessia Caramello/Alessia Lab Dropbox/Projects/UKDRI project/Image analysis/HistoCat:SIMPLI analysis/051022 TREM2 cohort/R code NatComm revision")

dir.create("tables/masks")

##### 1. Load all PNG files and extract X,Y and intensity #####

files <- list.files("ImageJ_macro/enlarged_output_files/", pattern = "*.png")

for (f in files){
  
  print(paste0("processing ", f))
  
  # load image
  img <- readPNG(paste0 ("ImageJ_macro/enlarged_output_files/", f))
  
  # Get image dimensions
  dim_x <- dim(img)[1]  # Number of rows (height)
  dim_y <- dim(img)[2]  # Number of columns (width)
  
  # Create a data frame with X, Y, and intensity
  pixel_data <- data.frame(
    X = rep(1:dim_y, each = dim_x),
    Y = rep(1:dim_x, times = dim_y),
    value = as.vector(img) # Flatten the matrix into a vector
  )
  
  write_csv(pixel_data, paste0("tables/masks/coordinates/", gsub(".png", "", f), ".csv"))
  
}









