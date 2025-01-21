# Calculate number of cells in plaques and rings

# input files:
# - clustered_cells.csv (SIMPLI output)
# - csv. files of plaques and rings masks (generated in 0b. pixel coordinates extraction from masks.R)
# - clusters_labels.csv (metadata)

  
##### 0. Install packages and set working directory #####

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
library(ggplot2)
library(dplyr)
library(writexl)

# set working directory
setwd("~/R code")

# read in clustered_cells
clustered_cells <- read_csv("input/clustered_cells.csv")

# define resolution
res <- "res_0.9_ids"

##### 1. Merge csv coordinate files and filter for white pixels only #####

# Get list of CSV files in directory
file_list <- list.files(path = "./tables/masks/coordinates/",
                        recursive = TRUE,
                        pattern = ".csv")

replicates <- unique(clustered_cells$Metadata_sample_name)
samples <- unique(stringr::str_extract(replicates, "[^_]*_[^_]*"))

# Create an empty list to store data frames
df_list <- list()

# Loop over CSV files and read them into R
for (file in file_list) {
  #file <- file_list[1]
  # Read CSV file into R
  df <- read.csv(paste0("./tables/masks/coordinates/", file))
  
  if (nrow(df) > 1){

    if (str_extract(file, "[^-]+") %in% replicates){
      
      df$sample_replicate <- str_extract(file, "[^-]+")
      df$region <- gsub(".csv", "", strsplit(file, "_(?!.*_)", perl=TRUE)[[1]][2])
      df <- df[!df$Value == 0, ]
      df_list[[file]] <- df
    }
  }
  print(paste0(file, " uploaded"))
}

# Combine all data frames into a single data frame
combined_df <- do.call(rbind, df_list)

# remove colnames
rownames(combined_df) <- NULL
combined_df$Value <- NULL
rm(df_list)
rm(df)

# save
write_csv(combined_df, "tables/masks/white_pixels_coordinates.csv")

replicates_pixels <- unique(combined_df$sample_replicate)
length(replicates_pixels)
setdiff(replicates, replicates_pixels)

##### 2. Find cells with same coordinates as 255 pixels #####

#combined_df <- read_csv("tables/masks/white_pixels_coordinates.csv")

# select astrocytes cluster 
astrocytes <- clustered_cells[clustered_cells[[res]] == 2, ]
astrocytes <- astrocytes[, colnames(astrocytes) %in% c("Metadata_sample_name", "Location_Center_X", "Location_Center_Y")]
colnames(astrocytes) <- c("sample_replicate", "X", "Y")

# select microglia cluster
microglia <- clustered_cells[clustered_cells[[res]] == 15, ]
microglia <- microglia[, colnames(microglia) %in% c("Metadata_sample_name", "Location_Center_X", "Location_Center_Y")]
colnames(microglia) <- c("sample_replicate", "X", "Y")

# round X and Y coordinates to nearest integer 
astrocytes$X <- round(astrocytes$X)
astrocytes$Y <- round(astrocytes$Y)
microglia$X <- round(microglia$X)
microglia$Y <- round(microglia$Y)

# concatenate X and Y coordinates into a single string/column
astrocytes$XY <- paste(astrocytes$X, astrocytes$Y, sep = "_")
microglia$XY <- paste(microglia$X, microglia$Y, sep = "_")
combined_df$XY <- paste(combined_df$X, combined_df$Y, sep = "_")

# find rows that match in xy coordinates and sample_replicate name 
matches_astrocytes <- merge(combined_df, astrocytes[, c("XY", "sample_replicate")], by = c("XY", "sample_replicate"), all = FALSE)
matches_microglia <- merge(combined_df, microglia[, c("XY", "sample_replicate")], by = c("XY", "sample_replicate"), all = FALSE)

write.csv(astrocytes, 'tables/masks/astrocytes_XY_merged.csv')
write.csv(microglia, 'tables/masks/microglia_XY_merged.csv')
write.csv(matches_astrocytes, 'tables/masks/microglia_in_plaques_rings.csv')
write.csv(matches_microglia, 'tables/masks/microglia_in_plaques_rings.csv')

##### 3. Count cells #####

# add info to table
matches_astrocytes <- matches_astrocytes[,c(2,5)]
matches_microglia <- matches_microglia[,c(2,5)]

matches_astrocytes$sample_name <- stringr::str_extract(matches_astrocytes$sample_replicate, "[^_]*_[^_]*")
matches_microglia$sample_name <- stringr::str_extract(matches_microglia$sample_replicate, "[^_]*_[^_]*")

# count and average
count_astrocytes <- matches_astrocytes %>% 
  group_by(sample_name, region) %>%
  dplyr::summarise(cell_number = n()
            )

count_microglia <- matches_microglia %>% 
  group_by(sample_name, region) %>%
  dplyr::summarise(cell_number = n()
  )

# add missing samples
temp <- data.frame(matrix(ncol = 3, nrow = length(samples)))
colnames(temp) <- colnames(count_astrocytes)
temp$sample_name <- samples
temp$cell_number <- as.numeric(0)
regions <- c("ring", "plaques")

for (r in regions){
  #r = "ring"
  
  # subset for region
  temp2 <- count_astrocytes[count_astrocytes$region == r,]
  temp$region <- r
  
  # find missing samples
  temp3 <- anti_join(temp, temp2, by = "sample_name")
  
  # join dataframes
  count_astrocytes <- rbind(count_astrocytes, temp3)
  
}

for (r in regions){
  #r = "ring"
  
  # subset for region
  temp2 <- count_microglia[count_microglia$region == r,]
  temp$region <- r
  
  # find missing samples
  temp3 <- anti_join(temp, temp2, by = "sample_name")
  
  # join dataframes
  count_microglia <- rbind(count_microglia, temp3)
  
}

rm(temp, temp2, temp3)

count_microglia$group <- lapply(str_split(count_microglia$sample_name, "_"), "[", 1)
count_astrocytes$group <- lapply(str_split(count_astrocytes$sample_name, "_"), "[", 1)

count_microglia$group <- as.character(count_microglia$group)
count_astrocytes$group <- as.character(count_astrocytes$group)

# prepare for paper (total cells)
paper_astrocytes <- count_astrocytes %>% pivot_wider(names_from = region, 
                                                   values_from = cell_number)

paper_microglia <- count_microglia %>% pivot_wider(names_from = region, 
                                                     values_from = cell_number)

paper_astrocytes[is.na(paper_astrocytes)] <- as.numeric(0)
paper_microglia[is.na(paper_microglia)] <- as.numeric(0)

paper_astrocytes$group <- as.character(paper_astrocytes$group)
paper_microglia$group <- as.character(paper_microglia$group)

paper_astrocytes <- paper_astrocytes[order(paper_astrocytes$group), ]
paper_microglia <- paper_microglia[order(paper_microglia$group), ]

# prepare for paper (average + sd)
paper_astrocytes_avg <- count_astrocytes %>% 
  dplyr::group_by(group, region) %>%
  dplyr::summarise(cells = paste(round(mean(cell_number), 2), 
                                 round(sd(cell_number), 2), 
                                 sep = ' ± ')
  )

paper_microglia_avg <- count_microglia %>% 
  dplyr::group_by(group, region) %>%
  dplyr::summarise(cells = paste(round(mean(cell_number), 2), 
                                 round(sd(cell_number), 2), 
                                 sep = ' ± ')
  )

paper_astrocytes_avg$group <- as.character(paper_astrocytes_avg$group)
paper_microglia_avg$group <- as.character(paper_microglia_avg$group)

# save
write_xlsx(paper_astrocytes, "tables/masks/astrocytes_in_plaques_rings_total_paper.xlsx")
write_xlsx(paper_microglia, "tables/masks/microglia_in_plaques_rings_total_paper.xlsx")
write_xlsx(paper_astrocytes_avg, "tables/masks/astrocytes_in_plaques_rings_avg_paper.xlsx")
write_xlsx(paper_microglia_avg, "tables/maskss/microglia_in_plaques_rings_avg_paper.xlsx")
write_csv(count_astrocytes, "tables/masks/astrocytes_in_plaques_rings.csv")
write_csv(count_microglia, "tables/masks/microglia_in_plaques_rings.csv")

##### 4. Statistical analyses and plot ######

# define analyses 
analyses <- c("astrocytes", "microglia")

dats <- list(count_astrocytes, count_microglia)
names(dats) <-  c("astrocytes", "microglia")

samples_names <- unique(count_microglia$sample_name)

dir.create("tables/masks/stats")
dir.create("plots/masks/")

for (a in 1:length(analyses)){
  #a=1
  analysis <- analyses[a]
  dat <- dats[[analyses[a]]]
  dat$sample_name <- NULL

##### Test normality ########
# testing normality with a Shapiro-Wilk Test.
# source: https://www.statology.org/test-for-normality-in-r/

# rename colnames
colnames(dat)[which(colnames(dat) == 'region')] <- "cell_type"
colnames(dat)[which(colnames(dat) == 'cell_number')] <- "count"

cell_types <- unique(dat$cell_type)
groups <- unique(dat$group)

#prepare data table for result of normality test
tests_todo <- data.frame(matrix(ncol = 2))
colnames(tests_todo) <- c('cell_type', 'dist')
tests_todo_temp <- tests_todo

# test normality 
# group_by: group, layer and cell_type 
for (group in 1:length(groups)){
  skip_to_next <- FALSE
  tryCatch({
    #group = 1
    samples1 <- dat[dat$group == groups[group],]
    
    for (cell in 1:length(cell_types)){
      skip_to_next <- FALSE
      tryCatch({
        #cell = 1
        samples2 <- samples1[samples1$cell_type == cell_types[cell],]
        
        if ((shapiro.test(samples2$count)$p.value) < 0.05) {
          print(paste0(cell_types[cell], " in ", groups[group], " **NOT** normally distributed"))
          
          tests_todo_temp$cell_type <- cell_types[cell]
          tests_todo_temp$dist <- 'kruskal.wilcox'
          
          tests_todo <- rbind(tests_todo, tests_todo_temp)
          
        } 
        
        if ((shapiro.test(samples2$count)$p.value) > 0.05) {
          print(paste0(cell_types[cell], " in ", groups[group], " normally distributed"))
          
          tests_todo_temp$cell_type <- cell_types[cell]
          tests_todo_temp$dist <- 'anova.tukey'
          
          tests_todo <- rbind(tests_todo, tests_todo_temp)
        }
      }, error = function(e) {skip_to_next <<- TRUE})
      if(skip_to_next) { next } 
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
colnames(tests_todo_final) <- c('cell_type', 'test')
tests_todo_final_temp <- tests_todo_final 

# check which cell_types and groups are normally distributed
# and define stat test to use for analysis
for (cell in 1:length(cell_types)){
  samples <- filter(tests_todo, cell_type == cell_types[cell])
  if (nrow(samples) > 1) {
    if (all(samples$dist == "anova.tukey")){
      tests_todo_final_temp$cell_type <- cell_types[cell]
      tests_todo_final_temp$test <- "anova.tukey"
    } else {
      tests_todo_final_temp$cell_type <- cell_types[cell]
      tests_todo_final_temp$test <- "kruskal.wilcox"
    }
  } else {
    tests_todo_final_temp$cell_type <- cell_types[cell]
    tests_todo_final_temp$test <- "kruskal.wilcox"
  }
  tests_todo_final <- rbind(tests_todo_final, tests_todo_final_temp)
}

# remove first column and unnecessary variables
if (is.na(tests_todo_final[1,1])){
  tests_todo_final <- tests_todo_final[-1,] 
}
rm(tests_todo_final_temp)
rm(tests_todo_temp)

# only some data are normally distributed:
# --> ANOVA + Tukey on normally distributed data
# --> Kruskal-wallis + Wilcoxon test on NOT normally distributed data

##### ANOVA + Tukey (for normally distributed data) #####
# source: https://www.scribbr.com/statistics/anova-in-r/#:~:text=We%20can%20perform%20an%20ANOVA,levels%20of%20the%20independent%20variable.

# filter samples for ANOVA
temp <- filter(tests_todo_final, test == 'anova.tukey')

cells_anova <- temp$cell_type

samples_anova <- dat %>% filter(cell_type %in% cells_anova)

# generate stat.test file to add tukey p.adj of count, for plotting
stat.test <- data.frame(matrix(ncol = 7))
colnames(stat.test) <- c('cell_type', 'group1', 'group2', 'p.adj', 'method', 'anova_krus', 'y.position')

if (nrow(samples_anova) > 0){
  # check we have all cells
  if (nrow(samples_anova) == length(cells_anova)*length(sample_names)){
    print("All cells present for ANOVA + Tukey analyses")
  }
  
  # ANOVA + Tukey test - grouping by cell_type
  for (cell in 1:length(cells_anova)){
    # filter cell type
    samples <- filter(samples_anova, samples_anova$cell_type == !!cells_anova[cell])
    
    # save ANOVA results
    one.way <- aov(count ~ group, data = samples)
    
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
    tukey.one.way.df$cell_type <- cells_anova[cell]
    tukey.one.way.df$rowname <- NULL
    tukey.one.way.df$group.diff <- NULL
    tukey.one.way.df$group.lwr <- NULL
    tukey.one.way.df$group.upr <- NULL
    
    tukey.one.way.df <- tukey.one.way.df %>% relocate(cell_type, .before = group1)
    tukey.one.way.df <- tukey.one.way.df %>% relocate(p.adj, .after = group2)
    tukey.one.way.df$method <- 'anova.tukey'
    
    # add y.position and Pr(>F)
    y_increase <- max(samples$count)*0.2
    for (n in 1:nrow(tukey.one.way.df)){
      tukey.one.way.df$y.position[n] <- max(samples$count)+(y_increase+((y_increase/2)*(n-1)))
      tukey.one.way.df$anova_krus[n] <- summary(one.way)[[1]][["Pr(>F)"]][1]
    }
    
    # add to previous table
    stat.test <- rbind(stat.test, tukey.one.way.df)
    
  }
  
  # remove first row
  if (is.na(stat.test[1,1])){
    stat.test <- stat.test[-1,] 
  }
  
} else {
  print("No samples for ANOVA")
}

##### Kruskal + Wilcoxon (for not normally distributed data) #####

# filter samples for Kruskal
temp <- filter(tests_todo_final, test == 'kruskal.wilcox')

cells_kruskal <- temp$cell_type

samples_kruskal <- dat %>% filter(cell_type %in% cells_kruskal)

if (nrow(samples_kruskal) > 0){
  # check we have all cells
  
  if (nrow(samples_kruskal) == length(cells_kruskal)*length(samples_names)){
    print("All cells present for Kruskal + Wilcoxon analysis")
  }
  
  # Wilcox test - grouping by cell_type
  for (cell in 1:length(cells_kruskal)){
    #cell = 1
    # filter for cell type
    samples1 <- filter(samples_kruskal, samples_kruskal$cell_type == !!cells_kruskal[cell])
    
    if (nrow(samples1) > 2) {
      # perform Kruskal test and save
      kruskal.test <- kruskal.test(count ~ group, data = samples1)
      
      # perform Wilcox test and save
      wilcox.test <- pairwise.wilcox.test(samples1$count, samples1$group,
                                          p.adjust.method = "BH", 
                                          exact = FALSE)
      
      # modify wilcox.test to fit stat.test
      wilcox.test <- rownames_to_column(as.data.frame(wilcox.test$p.value))
      wilcox.test <- melt(wilcox.test)
      colnames(wilcox.test) <- c('group2', 'group1', 'p.adj')
      wilcox.test$method <- 'kruskal.wilcox'
      wilcox.test$cell_type <- cells_kruskal[cell]
      
      wilcox.test <- wilcox.test %>% relocate(cell_type, .before = group1)
      wilcox.test <- wilcox.test %>% relocate(p.adj, .after = group2)
      wilcox.test <- wilcox.test %>% relocate(group2, .after = group1)
      
      # add y.position and Pr(>F) and remove rows from same samples
      y_increase <- max(samples1$count)*0.2
      for (n in 1:nrow(wilcox.test)){
        wilcox.test$y.position[n] <- max(samples1$count)+(y_increase+((y_increase/2)*(n-1)))
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
  
  # remove first row
  if (is.na(stat.test[1,1])){
    stat.test <- stat.test[-1,] 
  }
  
} else {
  print("No samples for Kruskal")
}

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

#check we have all cells
length(unique(stat.test$cell_type)) == length(cell_types)

# cell types missing
setdiff(cell_types, unique(stat.test$cell_type))

# remove unnecessary comparisons
stat.test <- stat.test[!(stat.test$group1 == "CtrlTREM2" & stat.test$group2 == "AlzCV"), ]
stat.test <- stat.test[!(stat.test$group1 == "AlzCV" & stat.test$group2 == "CtrlTREM2"), ]

# save file
write_csv(stat.test, 
          paste0("tables/masks/stats/stats_", analyses[a] ,"_in_plaques_rings.csv"))

##### Plot ######

dat_sel <- dat

# filter p-adj that are not significant
if (any(stat.test$p.adj <= 0.05)){
  stat.test.signif <- filter(stat.test, p.adj <= 0.05)
} else {
  stat.test.signif <- stat.test[!is.na(stat.test$label),]
}

# add column for y.position
stat.test.signif$y.position <- 0

# adjust y.position in stat.test.signif
# get markers with significant difference
markers <- unique(stat.test.signif$cell_type)

# create temporary stat.test.signif_temp file to store new y.position
stat.test.signif_temp <- data.frame(matrix(ncol = ncol(stat.test.signif)))
colnames(stat.test.signif_temp) <- colnames(stat.test.signif)

# calculate new y.position
for (marker in 1:length(markers)) {
  temp <- dat_sel[str_detect(dat_sel$cell_type, markers[marker]),]
  temp2 <- stat.test.signif[stat.test.signif$cell_type == markers[marker],]
  for (i in 1:nrow(temp2)){
    temp2$y.position[i] <- max(temp$count)+((max(dat_sel$count)/10)*i)
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

# read in clusters label
cluster_labels <- read_csv("input/clusters_labels.csv")

# re-order groups order
x <- c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzTREM2")
dat_sel <- dat_sel %>% arrange(sapply(group, function(y) which(y == x)))
#stat.test.signif <- stat.test.signif %>% arrange(sapply(group, function(y) which(y == x)))


# pick palette
if (any(x == "AlzTREM2")){
  print("AlzTREM2")
  sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898")
} else if (any(x %in% "AlzR47H")){
  print("AlzR47H")
  sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898", "#D8C5C5")
} else {
  print("AlzCV")
  sel_palette = c("#314B93", "#D02D27")
}

# plot
plot <- ggboxplot(dat_sel, x = "group", y = "count",
                  color = "group", 
                  palette = sel_palette,
                  add = "jitter", 
                  ylab = paste0("total ", analyses[a], " number") ,
                  xlab = "",
                  short.panel.labs = TRUE) +
  facet_wrap(~factor(cell_type)) +
  theme(legend.title = element_blank(), axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) +
  ylim(NA,ylim_plot+5)+
  rremove("legend")

# Add stats
plot <- plot + stat_pvalue_manual(stat.test.signif,
                                  label = "label", 
                                  label.size = 5,
                                  tip.length = 0)
plot

pdf(paste0("plots/masks/", 
           analyses[a], 
           "_in_plaques_rings.pdf"), 
    width = 3.5, 
    height = 4)
print(plot)
dev.off()



}


##### 5. Calculate cell density #####

# calculate total plaque and ring area
count_pixels <- combined_df %>% 
  group_by(sample_replicate, region) %>%
  dplyr::summarise(pixels = n()
  )

count_pixels$sample_name <- stringr::str_extract(count_pixels$sample_replicate, "[^_]*_[^_]*")

count_pixels <- count_pixels %>% 
  group_by(sample_name, region) %>%
  dplyr::summarise(pixels = sum(pixels)
  )

# add missing samples
samples <- unique(stringr::str_extract(replicates, "[^_]*_[^_]*"))
temp <- data.frame(matrix(ncol = 3, nrow = length(samples)))
colnames(temp) <- colnames(count_pixels)
temp$sample_name <- samples
temp$pixels <- as.numeric(0)
regions <- c("ring", "plaques")

for (r in regions){
  #r = "ring"
  
  # subset for region
  temp2 <- count_pixels[count_pixels$region == r,]
  temp$region <- r
  
  # find missing samples
  temp3 <- anti_join(temp, temp2, by = "sample_name")
  
  # join dataframes
  count_pixels <- rbind(count_pixels, temp3)
  
}

# add mm2 pixels area
count_pixels$pixels_mm2 <- count_pixels$pixels / 1000000 # mm2

# add pixel mm2 area to microglia and astrocytes
count_astrocytes <- left_join(count_astrocytes, count_pixels, by = c("sample_name", "region"))
count_microglia <- left_join(count_microglia, count_pixels, by = c("sample_name", "region"))

# calculate cell density
count_astrocytes$cell_density <- count_astrocytes$cell_number / count_astrocytes$pixels_mm2
count_microglia$cell_density <- count_microglia$cell_number / count_microglia$pixels_mm2

count_astrocytes$cell_density[is.nan(count_astrocytes$cell_density)] <- 0
count_microglia$cell_density[is.nan(count_microglia$cell_density)] <- 0

write_csv(count_astrocytes, "tables/masks/astrocytes_in_plaques_rings.csv")
write_csv(count_microglia, "tables/masks/microglia_in_plaques_rings.csv")

##### 6. Statistical analyses and plot ######

# define analyses 
analyses <- c("astrocytes", "microglia")

dats <- list(count_astrocytes, count_microglia)
names(dats) <-  c("astrocytes", "microglia")

samples_names <- unique(count_microglia$sample_name)

for (a in 1:length(analyses)){
  #a=1
  analysis <- analyses[a]
  dat <- dats[[analyses[a]]]
  dat$sample_name <- NULL
  
  ##### Test normality ########
  # testing normality with a Shapiro-Wilk Test.
  # source: https://www.statology.org/test-for-normality-in-r/
  
  # rename colnames
  colnames(dat)[which(colnames(dat) == 'region')] <- "cell_type"
  colnames(dat)[which(colnames(dat) == 'cell_density')] <- "count"
  
  dat <- dat[colnames(dat) %in% c("cell_type", "count", "group")]
  
  cell_types <- unique(dat$cell_type)
  groups <- unique(dat$group)
  
  #prepare data table for result of normality test
  tests_todo <- data.frame(matrix(ncol = 2))
  colnames(tests_todo) <- c('cell_type', 'dist')
  tests_todo_temp <- tests_todo
  
  # test normality 
  # group_by: group, layer and cell_type 
  for (group in 1:length(groups)){
    skip_to_next <- FALSE
    tryCatch({
      #group = 1
      samples1 <- dat[dat$group == groups[group],]
      
      for (cell in 1:length(cell_types)){
        skip_to_next <- FALSE
        tryCatch({
          #cell = 1
          samples2 <- samples1[samples1$cell_type == cell_types[cell],]
          
          if ((shapiro.test(samples2$count)$p.value) < 0.05) {
            print(paste0(cell_types[cell], " in ", groups[group], " **NOT** normally distributed"))
            
            tests_todo_temp$cell_type <- cell_types[cell]
            tests_todo_temp$dist <- 'kruskal.wilcox'
            
            tests_todo <- rbind(tests_todo, tests_todo_temp)
            
          } 
          
          if ((shapiro.test(samples2$count)$p.value) > 0.05) {
            print(paste0(cell_types[cell], " in ", groups[group], " normally distributed"))
            
            tests_todo_temp$cell_type <- cell_types[cell]
            tests_todo_temp$dist <- 'anova.tukey'
            
            tests_todo <- rbind(tests_todo, tests_todo_temp)
          }
        }, error = function(e) {skip_to_next <<- TRUE})
        if(skip_to_next) { next } 
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
  colnames(tests_todo_final) <- c('cell_type', 'test')
  tests_todo_final_temp <- tests_todo_final 
  
  # check which cell_types and groups are normally distributed
  # and define stat test to use for analysis
  for (cell in 1:length(cell_types)){
    samples <- filter(tests_todo, cell_type == cell_types[cell])
    if (nrow(samples) > 1) {
      if (all(samples$dist == "anova.tukey")){
        tests_todo_final_temp$cell_type <- cell_types[cell]
        tests_todo_final_temp$test <- "anova.tukey"
      } else {
        tests_todo_final_temp$cell_type <- cell_types[cell]
        tests_todo_final_temp$test <- "kruskal.wilcox"
      }
    } else {
      tests_todo_final_temp$cell_type <- cell_types[cell]
      tests_todo_final_temp$test <- "kruskal.wilcox"
    }
    tests_todo_final <- rbind(tests_todo_final, tests_todo_final_temp)
  }
  
  # remove first column and unnecessary variables
  if (is.na(tests_todo_final[1,1])){
    tests_todo_final <- tests_todo_final[-1,] 
  }
  rm(tests_todo_final_temp)
  rm(tests_todo_temp)
  
  # only some data are normally distributed:
  # --> ANOVA + Tukey on normally distributed data
  # --> Kruskal-wallis + Wilcoxon test on NOT normally distributed data
  
  ##### ANOVA + Tukey (for normally distributed data) #####
  # source: https://www.scribbr.com/statistics/anova-in-r/#:~:text=We%20can%20perform%20an%20ANOVA,levels%20of%20the%20independent%20variable.
  
  # filter samples for ANOVA
  temp <- filter(tests_todo_final, test == 'anova.tukey')
  
  cells_anova <- temp$cell_type
  
  samples_anova <- dat %>% filter(cell_type %in% cells_anova)
  
  # generate stat.test file to add tukey p.adj of count, for plotting
  stat.test <- data.frame(matrix(ncol = 7))
  colnames(stat.test) <- c('cell_type', 'group1', 'group2', 'p.adj', 'method', 'anova_krus', 'y.position')
  
  if (nrow(samples_anova) > 0){
    # check we have all cells
    if (nrow(samples_anova) == length(cells_anova)*length(sample_names)){
      print("All cells present for ANOVA + Tukey analyses")
    }
    
    # ANOVA + Tukey test - grouping by cell_type
    for (cell in 1:length(cells_anova)){
      # filter cell type
      samples <- filter(samples_anova, samples_anova$cell_type == !!cells_anova[cell])
      
      # save ANOVA results
      one.way <- aov(count ~ group, data = samples)
      
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
      tukey.one.way.df$cell_type <- cells_anova[cell]
      tukey.one.way.df$rowname <- NULL
      tukey.one.way.df$group.diff <- NULL
      tukey.one.way.df$group.lwr <- NULL
      tukey.one.way.df$group.upr <- NULL
      
      tukey.one.way.df <- tukey.one.way.df %>% relocate(cell_type, .before = group1)
      tukey.one.way.df <- tukey.one.way.df %>% relocate(p.adj, .after = group2)
      tukey.one.way.df$method <- 'anova.tukey'
      
      # add y.position and Pr(>F)
      y_increase <- max(samples$count)*0.2
      for (n in 1:nrow(tukey.one.way.df)){
        tukey.one.way.df$y.position[n] <- max(samples$count)+(y_increase+((y_increase/2)*(n-1)))
        tukey.one.way.df$anova_krus[n] <- summary(one.way)[[1]][["Pr(>F)"]][1]
      }
      
      # add to previous table
      stat.test <- rbind(stat.test, tukey.one.way.df)
      
    }
    
    # remove first row
    if (is.na(stat.test[1,1])){
      stat.test <- stat.test[-1,] 
    }
    
  } else {
    print("No samples for ANOVA")
  }
  
  ##### Kruskal + Wilcoxon (for not normally distributed data) #####
  
  # filter samples for Kruskal
  temp <- filter(tests_todo_final, test == 'kruskal.wilcox')
  
  cells_kruskal <- temp$cell_type
  
  samples_kruskal <- dat %>% filter(cell_type %in% cells_kruskal)
  
  if (nrow(samples_kruskal) > 0){
    # check we have all cells
    
    if (nrow(samples_kruskal) == length(cells_kruskal)*length(samples_names)){
      print("All cells present for Kruskal + Wilcoxon analysis")
    }
    
    # Wilcox test - grouping by cell_type
    for (cell in 1:length(cells_kruskal)){
      #cell = 1
      # filter for cell type
      samples1 <- filter(samples_kruskal, samples_kruskal$cell_type == !!cells_kruskal[cell])
      
      if (nrow(samples1) > 2) {
        # perform Kruskal test and save
        kruskal.test <- kruskal.test(count ~ group, data = samples1)
        
        # perform Wilcox test and save
        wilcox.test <- pairwise.wilcox.test(samples1$count, samples1$group,
                                            p.adjust.method = "BH", 
                                            exact = FALSE)
        
        # modify wilcox.test to fit stat.test
        wilcox.test <- rownames_to_column(as.data.frame(wilcox.test$p.value))
        wilcox.test <- melt(wilcox.test)
        colnames(wilcox.test) <- c('group1', 'group2', 'p.adj')
        wilcox.test$method <- 'kruskal.wilcox'
        wilcox.test$cell_type <- cells_kruskal[cell]
        
        wilcox.test <- wilcox.test %>% relocate(cell_type, .before = group1)
        wilcox.test <- wilcox.test %>% relocate(p.adj, .after = group2)
        wilcox.test <- wilcox.test %>% relocate(group2, .after = group1)
        
        # add y.position and Pr(>F) and remove rows from same samples
        y_increase <- max(samples1$count)*0.2
        for (n in 1:nrow(wilcox.test)){
          wilcox.test$y.position[n] <- max(samples1$count)+(y_increase+((y_increase/2)*(n-1)))
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
    
    # remove first row
    if (is.na(stat.test[1,1])){
      stat.test <- stat.test[-1,] 
    }
    
  } else {
    print("No samples for Kruskal")
  }
  
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
  
  #check we have all cells
  length(unique(stat.test$cell_type)) == length(cell_types)
  
  # cell types missing
  #setdiff(cell_types, unique(stat.test$cell_type))
  
  # remove unnecessary comparisons
  stat.test <- stat.test[!(stat.test$group1 == "CtrlTREM2" & stat.test$group2 == "AlzCV"), ]
  stat.test <- stat.test[!(stat.test$group1 == "AlzCV" & stat.test$group2 == "CtrlTREM2"), ]
  
  # save file
  write_csv(stat.test, 
            paste0("tables/masks/stats/stats_", analyses[a] ,"_density_in_plaques_rings.csv"))
  
  ##### Plot ######
  
  dat_sel <- dat
  
  # filter p-adj that are not significant
  if (any(stat.test$p.adj <= 0.05)){
    stat.test.signif <- filter(stat.test, p.adj <= 0.05)
  } else {
    stat.test.signif <- stat.test[!is.na(stat.test$label),]
  }
  
  # add column for y.position
  stat.test.signif$y.position <- 0
  
  # adjust y.position in stat.test.signif
  # get markers with significant difference
  markers <- unique(stat.test.signif$cell_type)
  
  # create temporary stat.test.signif_temp file to store new y.position
  stat.test.signif_temp <- data.frame(matrix(ncol = ncol(stat.test.signif)))
  colnames(stat.test.signif_temp) <- colnames(stat.test.signif)
  
  # calculate new y.position
  for (marker in 1:length(markers)) {
    temp <- dat_sel[str_detect(dat_sel$cell_type, markers[marker]),]
    temp2 <- stat.test.signif[stat.test.signif$cell_type == markers[marker],]
    for (i in 1:nrow(temp2)){
      temp2$y.position[i] <- max(temp$count)+((max(dat_sel$count)/10)*i)
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
  
  # re-order groups order
  x <- c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzTREM2")
  dat_sel <- dat_sel %>% arrange(sapply(group, function(y) which(y == x)))
  #stat.test.signif <- stat.test.signif %>% arrange(sapply(group, function(y) which(y == x)))
  
  # pick palette
  if (any(x == "AlzTREM2")){
    print("AlzTREM2")
    sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898")
  } else if (any(x %in% "AlzR47H")){
    print("AlzR47H")
    sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898", "#D8C5C5")
  } else {
    print("AlzCV")
    sel_palette = c("#314B93", "#D02D27")
  }
  
  # plot
  plot <- ggboxplot(dat_sel, x = "group", y = "count",
                    color = "group", 
                    palette = sel_palette,
                    add = "jitter", 
                    ylab = paste0("density of ", analyses[a]) ,
                    xlab = "",
                    short.panel.labs = TRUE) +
    facet_wrap(~factor(cell_type)) +
    theme(legend.title = element_blank(), axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) +
    ylim(NA,ylim_plot+5)+
    rremove("legend")
  
  # Add stats
  plot <- plot + stat_pvalue_manual(stat.test.signif,
                                    label = "label", 
                                    label.size = 5,
                                    tip.length = 0)
  plot
  
  pdf(paste0("plots/masks/", analyses[a], 
             "_density_in_plaques_rings.pdf"), 
      width = 3.5, 
      height = 4)
  print(plot)
  dev.off()
  
}



##### 7. Calculate cell proportion #####

# read files
count_astrocytes <- read_csv("tables/masks/astrocytes_in_plaques_rings.csv")
count_microglia <- read_csv("tables/masks/microglia_in_plaques_rings.csv")
cells_cluster <- read_csv("tables/clusters/cells_per_cluster_long.csv")

# calculate sum of cells per cluster
cells_cluster <- cells_cluster %>% 
  dplyr::group_by(sample_name, mark) %>%
  dplyr::summarise(tot_cells = sum(value))

# turn wide
cells_cluster <- cells_cluster %>% pivot_wider(names_from = mark, 
                                               values_from = tot_cells)

# extract total microglia and astrocytes
all_astro <- cells_cluster[,colnames(cells_cluster) %in% c("sample_name", "cluster_2")]
all_micro <- cells_cluster[,colnames(cells_cluster) %in% c("sample_name", "cluster_15")]

# add regions 
all_astro$region <- "plaques"
all_astro2 <- all_astro
all_astro2$region <- "ring"
all_astro <- rbind(all_astro, all_astro2)

all_micro$region <- "plaques"
all_micro2 <- all_micro
all_micro2$region <- "ring"
all_micro <- rbind(all_micro, all_micro2)

rm(all_astro2, all_micro2)

# add all astrocytes and microglia to datasets
count_astrocytes <- left_join(count_astrocytes, all_astro, by = c("sample_name", "region"))
colnames(count_astrocytes)[colnames(count_astrocytes) == "cluster_2"] <- "tot_cells"

count_microglia <- left_join(count_microglia, all_micro, by = c("sample_name", "region"))
colnames(count_microglia)[colnames(count_microglia) == "cluster_15"] <- "tot_cells"

# calculate cell proportion
count_astrocytes$cell_prop <- (count_astrocytes$cell_number / count_astrocytes$tot_cells) * 100
count_microglia$cell_prop <- (count_microglia$cell_number / count_microglia$tot_cells) * 100

write_csv(count_astrocytes, "tables/masks/astrocytes_in_plaques_rings.csv")
write_csv(count_microglia, "tables/masks/microglia_in_plaques_rings.csv")

###### save table for paper

# prepare for paper (total cells)
paper_astrocytes <- count_astrocytes[, colnames(count_astrocytes) %in% c("region", "cell_prop", "sample_name")]

paper_astrocytes <- paper_astrocytes %>% pivot_wider(names_from = region, 
                                                     values_from = cell_prop)

paper_microglia <- count_microglia[, colnames(count_microglia) %in% c("region", "cell_prop", "sample_name")]

paper_microglia <- paper_microglia %>% pivot_wider(names_from = region, 
                                                   values_from = cell_prop)

paper_microglia$group <- lapply(str_split(paper_microglia$sample_name, "_"), "[", 1)
paper_astrocytes$group <- lapply(str_split(paper_astrocytes$sample_name, "_"), "[", 1)

paper_astrocytes$group <- as.character(paper_astrocytes$group)
paper_microglia$group <- as.character(paper_microglia$group)

paper_astrocytes <- paper_astrocytes[order(paper_astrocytes$group), ]
paper_microglia <- paper_microglia[order(paper_microglia$group), ]

# prepare for paper (average + sd)
paper_astrocytes_avg <- count_astrocytes %>% 
  dplyr::group_by(group, region) %>%
  dplyr::summarise(cells = paste(round(mean(cell_prop), 2), 
                                 round(sd(cell_prop), 2), 
                                 sep = ' ± ')
  )

paper_microglia_avg <- count_microglia %>% 
  dplyr::group_by(group, region) %>%
  dplyr::summarise(cells = paste(round(mean(cell_prop), 2), 
                                 round(sd(cell_prop), 2), 
                                 sep = ' ± ')
  )

paper_astrocytes_avg$group <- as.character(paper_astrocytes_avg$group)
paper_microglia_avg$group <- as.character(paper_microglia_avg$group)

write_xlsx(paper_astrocytes, "tables/masks/astrocytes_in_plaques_rings_prop_paper.xlsx")
write_xlsx(paper_microglia, "tables/masks/microglia_in_plaques_rings_prop_paper.xlsx")
write_xlsx(paper_astrocytes_avg, "tables/masks/astrocytes_in_plaques_rings_avg_paper.xlsx")
write_xlsx(paper_microglia_avg, "tables/masks/microglia_in_plaques_rings_avg_paper.xlsx")

##### 8. Statistical analyses and plot ######

# define analyses 
analyses <- c("astrocytes", "microglia")

dats <- list(count_astrocytes, count_microglia)
names(dats) <-  c("astrocytes", "microglia")

samples_names <- unique(count_microglia$sample_name)

for (a in 1:length(analyses)){
  #a=1
  analysis <- analyses[a]
  dat <- dats[[analyses[a]]]
  
  ##### Transform percentage ########
  # explanation arcsin: https://www.programmingr.com/tutorial/arcsine-transformation/
  # explanation arcsin and logit with nice plots: http://strata.uga.edu/8370/rtips/proportions.html
  
  # save dat for plotting
  dat <- dat[colnames(dat) %in% c("region", "cell_prop", "group")]
  colnames(dat)[which(colnames(dat) == 'region')] <- "cell_type"
  colnames(dat)[which(colnames(dat) == 'cell_prop')] <- "count"
  dat_sel <- dat

  # add column with arcsin, log and radarcsin
  dat$arcsin_perc <- asin(sqrt(dat$count/100))
  
  # keep only columns of interest
  dat <- dat[colnames(dat) %in% c("cell_type", "arcsin_perc", "group")]
  colnames(dat)[which(colnames(dat) == 'arcsin_perc')] <- "count"
  
  # dat$log_perc <- log((((dat$count)+0.00001)/100)/
  #                       (1-((dat$count+0.00001)/100)))
  
  # dat$radarcsin_perc <- rad2deg(asin(sqrt(((dat$count)+0.00001)/100)))
  
  ##### Test normality ########
  # testing normality with a Shapiro-Wilk Test.
  # source: https://www.statology.org/test-for-normality-in-r/
  
  cell_types <- unique(dat$cell_type)
  groups <- unique(dat$group)
  
  #prepare data table for result of normality test
  tests_todo <- data.frame(matrix(ncol = 2))
  colnames(tests_todo) <- c('cell_type', 'dist')
  tests_todo_temp <- tests_todo
  
  # test normality 
  # group_by: group, layer and cell_type 
  for (group in 1:length(groups)){
    skip_to_next <- FALSE
    tryCatch({
      #group = 1
      samples1 <- dat[dat$group == groups[group],]
      
      for (cell in 1:length(cell_types)){
        skip_to_next <- FALSE
        tryCatch({
          #cell = 1
          samples2 <- samples1[samples1$cell_type == cell_types[cell],]
          
          if ((shapiro.test(samples2$count)$p.value) < 0.05) {
            print(paste0(cell_types[cell], " in ", groups[group], " **NOT** normally distributed"))
            
            tests_todo_temp$cell_type <- cell_types[cell]
            tests_todo_temp$dist <- 'kruskal.wilcox'
            
            tests_todo <- rbind(tests_todo, tests_todo_temp)
            
          } 
          
          if ((shapiro.test(samples2$count)$p.value) > 0.05) {
            print(paste0(cell_types[cell], " in ", groups[group], " normally distributed"))
            
            tests_todo_temp$cell_type <- cell_types[cell]
            tests_todo_temp$dist <- 'anova.tukey'
            
            tests_todo <- rbind(tests_todo, tests_todo_temp)
          }
        }, error = function(e) {skip_to_next <<- TRUE})
        if(skip_to_next) { next } 
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
  colnames(tests_todo_final) <- c('cell_type', 'test')
  tests_todo_final_temp <- tests_todo_final 
  
  # check which cell_types and groups are normally distributed
  # and define stat test to use for analysis
  for (cell in 1:length(cell_types)){
    samples <- filter(tests_todo, cell_type == cell_types[cell])
    if (nrow(samples) > 1) {
      if (all(samples$dist == "anova.tukey")){
        tests_todo_final_temp$cell_type <- cell_types[cell]
        tests_todo_final_temp$test <- "anova.tukey"
      } else {
        tests_todo_final_temp$cell_type <- cell_types[cell]
        tests_todo_final_temp$test <- "kruskal.wilcox"
      }
    } else {
      tests_todo_final_temp$cell_type <- cell_types[cell]
      tests_todo_final_temp$test <- "kruskal.wilcox"
    }
    tests_todo_final <- rbind(tests_todo_final, tests_todo_final_temp)
  }
  
  # remove first column and unnecessary variables
  if (is.na(tests_todo_final[1,1])){
    tests_todo_final <- tests_todo_final[-1,] 
  }
  rm(tests_todo_final_temp)
  rm(tests_todo_temp)
  
  # only some data are normally distributed:
  # --> ANOVA + Tukey on normally distributed data
  # --> Kruskal-wallis + Wilcoxon test on NOT normally distributed data
  
  ##### ANOVA + Tukey (for normally distributed data) #####
  # source: https://www.scribbr.com/statistics/anova-in-r/#:~:text=We%20can%20perform%20an%20ANOVA,levels%20of%20the%20independent%20variable.
  
  # filter samples for ANOVA
  temp <- filter(tests_todo_final, test == 'anova.tukey')
  
  cells_anova <- temp$cell_type
  
  samples_anova <- dat %>% filter(cell_type %in% cells_anova)
  
  # generate stat.test file to add tukey p.adj of count, for plotting
  stat.test <- data.frame(matrix(ncol = 7))
  colnames(stat.test) <- c('cell_type', 'group1', 'group2', 'p.adj', 'method', 'anova_krus', 'y.position')
  
  if (nrow(samples_anova) > 0){
    # check we have all cells
    if (nrow(samples_anova) == length(cells_anova)*length(sample_names)){
      print("All cells present for ANOVA + Tukey analyses")
    }
    
    # ANOVA + Tukey test - grouping by cell_type
    for (cell in 1:length(cells_anova)){
      # filter cell type
      samples <- filter(samples_anova, samples_anova$cell_type == !!cells_anova[cell])
      
      # save ANOVA results
      one.way <- aov(count ~ group, data = samples)
      
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
      tukey.one.way.df$cell_type <- cells_anova[cell]
      tukey.one.way.df$rowname <- NULL
      tukey.one.way.df$group.diff <- NULL
      tukey.one.way.df$group.lwr <- NULL
      tukey.one.way.df$group.upr <- NULL
      
      tukey.one.way.df <- tukey.one.way.df %>% relocate(cell_type, .before = group1)
      tukey.one.way.df <- tukey.one.way.df %>% relocate(p.adj, .after = group2)
      tukey.one.way.df$method <- 'anova.tukey'
      
      # add y.position and Pr(>F)
      y_increase <- max(samples$count)*0.2
      for (n in 1:nrow(tukey.one.way.df)){
        tukey.one.way.df$y.position[n] <- max(samples$count)+(y_increase+((y_increase/2)*(n-1)))
        tukey.one.way.df$anova_krus[n] <- summary(one.way)[[1]][["Pr(>F)"]][1]
      }
      
      # add to previous table
      stat.test <- rbind(stat.test, tukey.one.way.df)
      
    }
    
    # remove first row
    if (is.na(stat.test[1,1])){
      stat.test <- stat.test[-1,] 
    }
    
  } else {
    print("No samples for ANOVA")
  }
  
  ##### Kruskal + Wilcoxon (for not normally distributed data) #####
  
  # filter samples for Kruskal
  temp <- filter(tests_todo_final, test == 'kruskal.wilcox')
  
  cells_kruskal <- temp$cell_type
  
  samples_kruskal <- dat %>% filter(cell_type %in% cells_kruskal)
  
  if (nrow(samples_kruskal) > 0){
    # check we have all cells
    
    if (nrow(samples_kruskal) == length(cells_kruskal)*length(samples_names)){
      print("All cells present for Kruskal + Wilcoxon analysis")
    }
    
    # Wilcox test - grouping by cell_type
    for (cell in 1:length(cells_kruskal)){
      #cell = 1
      # filter for cell type
      samples1 <- filter(samples_kruskal, samples_kruskal$cell_type == !!cells_kruskal[cell])
      
      if (nrow(samples1) > 2) {
        # perform Kruskal test and save
        kruskal.test <- kruskal.test(count ~ group, data = samples1)
        
        # perform Wilcox test and save
        wilcox.test <- pairwise.wilcox.test(samples1$count, samples1$group,
                                            p.adjust.method = "BH", 
                                            exact = FALSE)
        
        # modify wilcox.test to fit stat.test
        wilcox.test <- rownames_to_column(as.data.frame(wilcox.test$p.value))
        wilcox.test <- melt(wilcox.test)
        colnames(wilcox.test) <- c('group1', 'group2', 'p.adj')
        wilcox.test$method <- 'kruskal.wilcox'
        wilcox.test$cell_type <- cells_kruskal[cell]
        
        wilcox.test <- wilcox.test %>% relocate(cell_type, .before = group1)
        wilcox.test <- wilcox.test %>% relocate(p.adj, .after = group2)
        wilcox.test <- wilcox.test %>% relocate(group2, .after = group1)
        
        # add y.position and Pr(>F) and remove rows from same samples
        y_increase <- max(samples1$count)*0.2
        for (n in 1:nrow(wilcox.test)){
          wilcox.test$y.position[n] <- max(samples1$count)+(y_increase+((y_increase/2)*(n-1)))
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
    
    # remove first row
    if (is.na(stat.test[1,1])){
      stat.test <- stat.test[-1,] 
    }
    
  } else {
    print("No samples for Kruskal")
  }
  
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
  
  #check we have all cells
  length(unique(stat.test$cell_type)) == length(cell_types)
  
  # cell types missing
  #setdiff(cell_types, unique(stat.test$cell_type))
  
  # remove unnecessary comparisons
  stat.test <- stat.test[!(stat.test$group1 == "CtrlTREM2" & stat.test$group2 == "AlzCV"), ]
  stat.test <- stat.test[!(stat.test$group1 == "AlzCV" & stat.test$group2 == "CtrlTREM2"), ]
  
  # save file
  write_csv(stat.test, 
            paste0("tables/masks/stats/stats_", analyses[a] ,"_prop_in_plaques_rings.csv"))
  
  ##### Plot ######
  
  # filter p-adj that are not significant
  if (any(stat.test$p.adj <= 0.05)){
    stat.test.signif <- filter(stat.test, p.adj <= 0.05)
  } else {
    stat.test.signif <- stat.test[!is.na(stat.test$label),]
  }
  
  # add column for y.position
  stat.test.signif$y.position <- 0
  
  # adjust y.position in stat.test.signif
  # get markers with significant difference
  markers <- unique(stat.test.signif$cell_type)
  
  # create temporary stat.test.signif_temp file to store new y.position
  stat.test.signif_temp <- data.frame(matrix(ncol = ncol(stat.test.signif)))
  colnames(stat.test.signif_temp) <- colnames(stat.test.signif)
  
  # calculate new y.position
  for (marker in 1:length(markers)) {
    temp <- dat_sel[str_detect(dat_sel$cell_type, markers[marker]),]
    temp2 <- stat.test.signif[stat.test.signif$cell_type == markers[marker],]
    for (i in 1:nrow(temp2)){
      temp2$y.position[i] <- max(temp$count)+((max(dat_sel$count)/10)*i)
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
  
  # re-order groups order
  x <- c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzTREM2")
  dat_sel <- dat_sel %>% arrange(sapply(group, function(y) which(y == x)))
  #stat.test.signif <- stat.test.signif %>% arrange(sapply(group, function(y) which(y == x)))
  
  # pick palette
  if (any(x == "AlzTREM2")){
    print("AlzTREM2")
    sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898")
  } else if (any(x %in% "AlzR47H")){
    print("AlzR47H")
    sel_palette = c("#314B93", "#8297C4","#D02D27", "#D39898", "#D8C5C5")
  } else {
    print("AlzCV")
    sel_palette = c("#314B93", "#D02D27")
  }
  
  # plot
  plot <- ggboxplot(dat_sel, x = "group", y = "count",
                    color = "group", 
                    palette = sel_palette,
                    add = "jitter", 
                    ylab = paste0("% of ", analyses[a]) ,
                    xlab = "",
                    short.panel.labs = TRUE) +
    facet_wrap(~factor(cell_type)) +
    theme(legend.title = element_blank(), axis.text.x = element_text(size = 9, angle = 45, hjust = 1)) +
    ylim(NA,ylim_plot+5)+
    rremove("legend")
  
  # Add stats
  plot <- plot + stat_pvalue_manual(stat.test.signif,
                                    label = "label", 
                                    label.size = 5,
                                    tip.length = 0)
  plot
  
  pdf(paste0("plots/masks/", analyses[a], 
             "_prop_in_plaques_rings.pdf"), 
      width = 3.5, 
      height = 4)
  print(plot)
  dev.off()
  
}

