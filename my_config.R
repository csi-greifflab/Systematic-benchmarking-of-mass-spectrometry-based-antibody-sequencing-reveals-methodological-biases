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
#library(DescTools)
library(circlize)
#library(googledrive)

project_path <- "/storage/mariiac/msms_figures"

source(file.path(project_path, "data/preprocessing_functions.R"))
source(file.path(project_path, "plots.R"))

enzyme_cols <- c("All" = "#575A6C",
                 "AspN" = "#B4C540",
                 "Tryp" = "#E84D8A",
                 "Ct" = "#64C5EB",
                 "Ct+Tryp" = "#7F58AF")



metadata_path <- file.path(project_path, "metadata")
data_path <- file.path(project_path, "data")

fig2_path <- file.path(project_path, "fig2")
fig5_path <- file.path(project_path, "fig5")
