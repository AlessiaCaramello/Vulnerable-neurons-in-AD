### Purpose of repository

This repository contains all R scripts used for performing data analysis for our project on neuronal vulnerability in Alzheimer's Disease (AD). A preprint of our study is available:

_Intra-cellular accumulation of amyloid is a marker of selective neuronal vulnerability in Alzheimer’s disease_
Alessia Caramello, Nurun Fancy, Clotilde Tournerie, Maxine Eklund, Vicky Chau, Emily Adair, Marianna Papageorgopoulou, Johanna Jackson, John Hardy, Paul M. Matthews
_medRxiv_ 2023.11.23.23298911; doi: https://doi.org/10.1101/2023.11.23.23298911

### Overview of the study

We used imaging mass cytometry (IMC) to identify the first neurons to be selectively lost in the human middle temporal gyrus in AD (vulnerable neurons). IMC images were processed using the software [SIMPLI](https://github.com/ciccalab/SIMPLI) 
(the version used is available in my GitHub page). Raw and processed images, as well as SIMPLI output files are available on my
[figshare](https://figshare.com/account/home#/projects/228897/edit) repository (still private while paper is under review). 

We then performed pathway analysis on vulnerable neuronal populations, using an [RNAseq dataset]([url](https://doi.org/10.1101/2022.07.12.22277509)) generated internally from the same sample cohort. Clusters from IMC and RNAseq were first matched based on size and markers expression.


### How to set up

1. **Gather all input files**:

   IMC data (available on Figshare):
   - area_measurements.csv
   - clustered_cells.csv
   - all_cells-clusters.csv
  
   snRNAseq data (provided):
   - RNAseq_cluster_size.tsv (cell number in each cluster)
   - IMC_markers_expression_in_RNAseq.csv (scaled expression per cluster of markers used in IMC)
   
   Metadata files (provided): 
   - clusters_labels.csv (IMC)
   - sample_rotation.csv (IMC)
   - vars_per_sample.csv (IMC)
   - markers_threshold.csv (IMC)
   - RNAseq_sample_metadata.tsv (RNAseq)

2. **Prepare the "R code" folder**:
   This should should include:
   - "input" -> all input files listed above
   - "R scripts" -> all scripts available in this repository
   - "ImageJ_macro" -> containing the "ImageJ macro/Plaques mask and rings.ijm" script and the "mask_channels" folder (which should contain all amyloid-b ROIs - not provided)
   - "plots" -> [empty] output for plots
   - "tables" -> [empty] output for tables
     
3. **Adapt metadata**:
   - clusters_labels.csv --> labels to assign to each cluster
   - sample_rotation.csv --> degree of clockwise rotation for turning each ROI in the correct orientation
   - vars_per_sample.csv --> additional variables that can be used to group samples for comparisons (sex, age, etc..)
   - markers_threshold.csv --> threshold at which a cells should be considered positive for a given marker
   - RNAseq_sample_metadata.tsv --> info for grouping RNAseq samples

4. **Run ImageJ macro**:
   This will generate masks of plaques and surrounding area (rings) masks from the amyloid-b channel (rings "thinckness" as ums of enlargment around plaques can be defined in the macro)



### Analysis

Scripts 0 to 7 --> IMC data analysis

Scripts 8 to 12 --> IMC-RNAseq clusters matching and pathway analysis

///

SCRIPT: **0a. flipping XY coord of 90 180 270 degrees.R**

Type of analysis: Flip images to the correct orientation by changing XY coordinates (by rotating to either 90, 180 or 270 degrees clockwise)
Input data: 
   - clustered_cells.csv (SIMPLI output)
   - sample_rotation.csv (metadata - degrees of rotation for each image)

///

SCRIPT: **0b. pixel coordinates extraction from masks.R**

Type of analysis: Extract pixel coordinates from PNG files generated using the ImageJ macro
Input data: 
   - PNG files of plaques and rings masks from ImageJ macro

///

SCRIPT: **1. channel pixel area.R**

Type of analysis: Calculate average channel intensity (as pixel area positive for marker) per sample

Input data: 
   - area_measurement.csv (SIMPLI output)
   - vars_per_sample.csv (metadata)

///

SCRIPT: **2a. heatmap of clusters markers**

Type of analysis: Generate heatmap of markers expression in each cluster

Input data: 
   - clustered_cells.csv (SIMPLI output)

///

SCRIPT: **2b. UMAPs of clusters.R**

Type of analysis: Generate UMAP plots highlighting samples, markers and clusters (adapted from SIMPLI plotting functions)

Input data: 
   - clustered_cells.csv (SIMPLI output)

///

SCRIPT: **3. cells number or density per cluster.R**

Type of analysis: Calculate number and density of cells per cluster, perform statistical analysis and plotting

Input data: 
   - clustered_cells.csv (SIMPLI output)
   - clusters_labels.csv (metadata)

///

SCRIPT: **4. clusters distribution density.R**

Type of analysis: Calculate distribution of cell density from single clusters in cortical layers, perform statistical analysis and plotting

Input data: 
   - clustered_cells.csv (SIMPLI output)
   - vars_per_sample.csv (metadata)
   - clusters_labels.csv (metadata)

///

SCRIPT: **5. clustered cells positive for markers of interest.R**

Type of analysis: Calculate number or % of cells positive for markers of interest in each cluster

Input data: 
   - clustered_cells.csv (SIMPLI output)
   - markers_threshold.csv (metadata)
   - vars_per_sample.csv (metadata)
   - clusters_labels.csv (metadata)

///

SCRIPT: **6. clusters spatial interactions.R**

Type of analysis: Calculate spatial distirbution between clusters

Input data: 
   - all_cells-clusters.csv (SIMPLI output)
   - UMAPs.rds (generated in 2b. UMAPs of clusters)
   - clusters_labels.csv (metadata)

///

SCRIPT: **7. analysis of cells in plaques and rings.R**

Type of analysis: Calculate number, density and proportion of cells in plaques and rings masked area

Input data: 
   - clustered_cells.csv (SIMPLI output)
   - csv. files of plaques and rings masks (generated in 0b. pixel coordinates extraction from masks.R)
   - clusters_labels.csv (metadata)

///

SCRIPT: **8. calculate_clusters_size_quartiles.R**

Type of analysis: Calculate IMC to RNAseq clusters matching score based on size (cell number)

Input data: 
   - RNAseq_sample_metadata.tsv (metadata)
   - RNAseq_cluster_size.tsv (snRNAseq data)
   - cells_per_cluster_wide.csv (generated in 3. cells number or density per cluster)

///

SCRIPT: **9. Matching of IMC clusters to snRNASeq clusters.R**

Type of analysis:Calculate final match of IMC clusters to snRNASeq clusters based on markers expression and clusters size

Input data: 
   - IMC_markers_expression_in_RNAseq.csv (snRNAseq data)
   - clustered_cells.csv (SIMPLI output files)
   - clusters_labels.csv (metadata)
   - clusters_quartile_matching_score.csv (generated in 8.calculate_clusters_size_quartiles)

