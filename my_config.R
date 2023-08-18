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


source("/storage/mariiac/msms_paper/data/preprocessing_functions.R")
source("/storage/mariiac/msms_paper/plots.R")

enzyme_cols <- c("All" = "#575A6C",
                 "AspN" = "#B4C540",
                 "Tryp" = "#E84D8A",
                 "Ct" = "#64C5EB",
                 "Ct+Tryp" = "#7F58AF")

metadata_path <- "/storage/mariiac/msms_paper/metadata"
data_path <- "/storage/mariiac/msms_paper/data"

fig2_path <- "/storage/mariiac/msms_paper/fig2"
fig5_path <- "/storage/mariiac/msms_paper/fig5"
