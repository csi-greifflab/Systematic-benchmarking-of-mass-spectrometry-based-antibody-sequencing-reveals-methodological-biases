library(Rcpp)
library(ggplot2)
library(readr)
library(stringr)
library(reshape2) 
library(readxl)
library(rio)
library(seqinr)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(forcats)
library(circlize)
#library(googledrive)
library(combinat)

project_path <- "/storage/mariiac/msms_figures"

source(file.path(project_path, "data/preprocessing_functions.R"))
source(file.path(project_path, "plots.R"))
source(file.path(project_path, "utils.R"))

enzyme_cols <- c("All" = "#575A6C",
                 "AspN" = "#B4C540",
                 "Tryp" = "#E84D8A",
                 "Ct" = "#64C5EB",
                 "Ct+Tryp" = "#7F58AF")

HC_names <- c("h9C12-Q97A_HC", "h9C12-WT_HC", "Brimab_HC", "Umab_HC", "PGT121_HC", "PGDM1400_HC")
  
LC_names <- c("h9C12_LC", "Brimab_LC", "Umab_LC", "PGT121_LC", "PGDM1400_LC")
  
HC_and_LC_names <- c("h9C12-Q97A_HC", "h9C12_LC", "h9C12-WT_HC", "Brimab_HC", "Brimab_LC", "Umab_HC", "Umab_LC", 
                     "PGT121_HC", "PGT121_LC", "PGDM1400_HC", "PGDM1400_LC")

metadata_path <- file.path(project_path, "metadata")
data_path <- file.path(project_path, "data")

fig2_path <- file.path(project_path, "fig2")
fig3_path <- file.path(project_path, "fig3")
fig5_path <- file.path(project_path, "fig5")
fig6_path <- file.path(project_path, "fig6")
