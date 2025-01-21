# Calculate average channel intensity per sample

# Input data: 
# - area_measurement.csv (SIMPLI output)

#### 0. Load libraries ####

library(stringr) 
library(tidyverse)
library(plyr)
library(readr)
library(stringr)
#install.packages("rqdatatable")
library(rqdatatable)
library(dplyr)
#library(ggExtra)
library(gridExtra)
library(ggsignif)
#install.packages("ggpubr")
library(ggpubr)
library(rstatix)
library(plyr)
#install.packages("reshape2")
library(reshape2)
#install.packages("writexl")
library(writexl)
library(readxl)
library(openxlsx)

setwd("~/R code")

dir.create("plots/pixel_area")
dir.create("tables/pixel_area")

#### 1. Load data from previous analysis ####
simpli_output <- read.csv('input/area_measurements.csv')

#### 2. Prepare table for analysis ####

# add rows with sample_name and group
colnames(simpli_output)[1] <- 'sample_replicate'

# extract sample name and group
simpli_output$sample_name = str_extract(simpli_output$sample_replicate, "[^_]*_[^_]*")
simpli_output$group = sub("_.*", "", simpli_output$sample_replicate)

# exception: run when you have only 1 sample per group
simpli_output$sample_name <- paste0(simpli_output$group, '_', sapply(strsplit(simpli_output$sample_replicate, "_"), "[[", 3))

#### 3. Calculate average percentage of marker+ area per total ROI area ####

# calculate average per sample -> avg of roi sample
avg_areal <- aggregate(simpli_output$percentage, 
                           list(simpli_output$sample_name,
                                simpli_output$marker 
                           ), 
                           mean)
# rename columns
colnames(avg_areal) <- c("sample_name", "marker", "perc_area")

# add column with group
avg_areal <- 
  add_column(avg_areal, 
             "group" = stringr::str_extract(avg_areal$sample_name, 
                                            "[^_]*"), .before = "sample_name")

write_csv(avg_areal,'tables/pixel_area/avg_perc_area_per_sample.csv')


#### 4. Calculate average total marker+ area and save ##### 

# select columns of interest
total_area_long <- subset(total_ROI_area, select = c('sample_replicate', 'sample_name', 'group', 'marker', 'area'))

# turn into wide table for calculations
total_area_wide <- reshape(total_area_long,
                          idvar = c("sample_name", "sample_replicate"),
                          v.names = "area",
                          timevar = "marker",
                          direction = "wide")

# extract markers
markers <- unique(total_ROI_area$marker)


# change colnames
colnames(total_area_wide) <- sub("area.","", colnames(total_area_wide))


# calculate average of column intensity (per sample)
avg_tot_area <- aggregate(total_area_wide[, 4:ncol(total_area_wide)],
                     list(total_area_wide$sample_name), mean)

names(avg_tot_area)[names(avg_tot_area) == 'Group.1'] <- 'sample_name'

# save csv
write_csv(avg_tot_area,'tables/pixel_area/avg_total_marker_area.csv')

# calculate average count of cells per sample -> avg of roi sample
avg_tot_areal <- aggregate(total_area_long$area,
                                    list(total_area_long$sample_name,
                                         total_area_long$marker
                                    ),
                                    mean)
# rename columns
colnames(avg_tot_areal) <- c("sample_name", "marker", "area")

# add column with group
avg_tot_areal <-
  add_column(avg_tot_areal,
             "group" = stringr::str_extract(avg_tot_areal$sample_name,
                                            "[^_]*"), .before = "sample_name")

# save csv
write_csv(avg_tot_areal,'tables/avg_total_marker_area_long.csv')

#### 5. Add variants #####

# load input file
#avg_areal <- read_csv('tables/pixel_area/avg_perc_area_per_sample.csv')

# load file with variants
variants <- read_csv('input/vars_per_sample.csv')

# add columns for hosting variants
avg_tot_areal$TREM2var <- "CV"
avg_tot_areal$ApoE4pos <- "E3/E3"
avg_tot_areal$AD_stage <- "a"
avg_tot_areal$TREM2var_ADstage <- "a"
avg_tot_areal$ApoE4pos_TREM2 <- "a"
avg_tot_areal$ApoE4pos_TREM2var <- "a"

# Add TREM2 variants
for (i in 1:nrow(avg_tot_areal)){
  avg_tot_areal$TREM2var[i] <- variants$TREM2var[grep(avg_tot_areal$sample_name[i], variants$sample_name)]
  avg_tot_areal$AD_stage[i] <- variants$AD_stage[grep(avg_tot_areal$sample_name[i], variants$sample_name)]
  avg_tot_areal$ApoE4pos[i] <- variants$ApoE4pos[grep(avg_tot_areal$sample_name[i], variants$sample_name)]
  avg_tot_areal$TREM2var_ADstage[i] <- variants$TREM2var_ADstage[grep(avg_tot_areal$sample_name[i], variants$sample_name)]
  avg_tot_areal$ApoE4pos_TREM2var[i] <- variants$ApoE4pos_TREM2var[grep(avg_tot_areal$sample_name[i], variants$sample_name)]
}

#### 6. Plot all markers ####

# calculate stats
stat.test <- avg_areal %>% 
  group_by(marker) %>%
  t_test(perc_area ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_significance("p")

stat.test <- stat.test %>%
  add_xy_position(fun = "max", x = "marker", dodge = 0.8)

# plot

plot <- ggbarplot(
  avg_areal, 
  x = "marker", 
  y = "perc_area", 
  xlab = "", 
  ylab = "Average pixel area",
  title = "Pixel area per marker",
  add = c("mean_sd", "jitter"), 
  color = "group", 
  palette = c("#D02D27", "#D39898", "#314B93", "#8297C4"), 
  position = position_dodge(0.8)
) + 
  scale_y_continuous(labels = scales::comma,  breaks = scales::extended_breaks(n = 10)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  # stat_pvalue_manual(
  #   stat.test,  
  #   label = "p.signif", 
  #   tip.length = 0.01
  #)

pdf("plots/pixel_area/average_marker_area.pdf", width = 12)
print(plot)
dev.off()


#### 7. Plot selected markers (OPEN -->) ####

# input markers to plot
selected_markers = "Ab|pTau"


#### 8. Filter table ####

# select markers to plot
avg_areal_sel <- avg_areal %>% 
  filter(str_detect(avg_areal$marker, selected_markers))

#### 9. Define comparisons and folders (OPEN -->) #####

# Create table with conditions for for loop 
conditions <- data.frame(matrix(ncol = 3, nrow = 3))
conditions <- data.frame(
  cond_name = c(1, 2, 3),
  order = I(list(c("CtrlCV", "AlzCV"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzTREM2"), c("CtrlCV", "CtrlTREM2", "AlzCV", "AlzR62H", "Alz472H"))),
  group_column = c("group", "group", "TREM2var")
)

# define column with value to test
column_oi <- 'perc_area'

# define if value is a percentage
is_percentage <- TRUE

# define paths to folders
folder <- 'pixel_area'
subfolder <- "normal_comp"
#subsubfolder <- 'add_if_necessary'

#### 10. Statistical analysis and plotting #####

# make directories
dir.create(paste0("plots/", folder))
dir.create(paste0("tables/", folder))

dir.create(paste0("plots/", folder, "/", subfolder))
dir.create(paste0("tables/", folder, "/", subfolder))

# define tables and plots path
folder_tables <- paste0("tables/", folder, "/", subfolder) 
folder_plots <- paste0("plots/", folder, "/", subfolder)

# run if you have a subsubfolder

# dir.create(paste0("plots/", folder, "/", subfolder, "/", subsubfolder))
# dir.create(paste0("tables/", folder, "/", subfolder, "/", subsubfolder))

# folder_tables <- paste0("tables/", folder, "/", subfolder, "/", subsubfolder, "/") 
# folder_plots <- paste0("plots/", folder, "/", subfolder, "/", subsubfolder, "/")


##### Start for loop 

for (cond in 1:nrow(conditions)) {
  #cond = 1
  
  print(paste0(" ----- START CONDITION #", cond, " ------"))
  
  dat <- as.data.frame(avg_areal_sel)
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
                     ylab = '% of pixel positive area',
                     short.panel.labs = TRUE, 
                     legend = NULL) +
        facet_wrap(~ marker,
                   ncol = 6,
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
                 "/pixel_area_", 
                 selected_markers, 
                 "_in_", 
                 paste(x, collapse = "_"), 
                 ".pdf"), 
          width = 4, 
          height = 3)
      print(p)
      dev.off()
      
      write_csv(stat.test, paste0(folder_tables, 
                             "/stats_pixel_area_", 
                             selected_markers, 
                             "_in_", 
                             paste(x, collapse = "_"), 
                             ".csv")) 
      
      # print plot done
      print(paste0("Plot for selected markers ", selected_markers, " was generated"))
  
  
  print(paste0(" ----- END CONDITION #", cond, " ------"))
  
}
