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


source("/storage/mariiac/MSMS_paper_figures/data/preprocessing_functions.R")
source("/storage/mariiac/MSMS_paper_figures/plots.R")

enzyme_cols <- c("all" = "#575A6C",
                 "aspn" = "#B4C540",
                 "tryp" = "#E84D8A",
                 "ct" = "#64C5EB",
                 "ct+tryp" = "#7F58AF")

metadata_path <- "/storage/mariiac/MSMS_paper_figures/metadata"
data_path <- "/storage/mariiac/MSMS_paper_figures/data"

fig2_path <- "/storage/mariiac/MSMS_paper_figures/fig2"
fig5_path <- "/storage/mariiac/MSMS_paper_figures/fig5"
