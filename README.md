### Purpose of repository

This repository contains all R scripts used for performing data analysis for our project on neuronal vulnerability in Alzheimer's Disease (AD). A preprint of our study is available:

_Intra-cellular accumulation of amyloid is a marker of selective neuronal vulnerability in Alzheimer’s disease_
Alessia Caramello, Nurun Fancy, Clotilde Tournerie, Maxine Eklund, Vicky Chau, Emily Adair, Marianna Papageorgopoulou, Johanna Jackson, John Hardy, Paul M. Matthews
_medRxiv_ 2023.11.23.23298911; doi: https://doi.org/10.1101/2023.11.23.23298911

### Overview of the study

We used imaging mass cytometry (IMC) to identify the first neurons to be lost in human AD. IMC images were processed using the software [SIMPLI](https://github.com/ciccalab/SIMPLI) 
(the version used is available in my GitHub page). Raw and processed images, as well as SIMPLI output files are available on my
[figshare](https://figshare.com/account/home#/projects/228897/edit) repository (still private while paper is under review). 


### How to set up

1. Prepare "R code" folder.
   This should include:
   - "input" -> all necessary input files: SIMPLI output files (area_measurements.csv, clustered_cells.csv), metadata (clusters_labels.csv, sample_rotation.csv)
   - "R scripts" -> scripts in this repository
   - "plots" -> [empty] output for plots
   - "tables" -> [empty] output for tables
