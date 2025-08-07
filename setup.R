# Setup packages

required_packages <- c(
  "ggplot2",
  "readr",
  "stringr",
  "seqinr",
  "ComplexHeatmap",
  "circlize",
  "dplyr",
  "tidyr",
  "reshape2",
  "RColorBrewer",
  "here",
  "DescTools", 
  "rstatix", 
  "ggpubr", 
  "patchwork", 
  "ggbeeswarm", 
  "RVAideMemoire"
)
invisible(lapply(required_packages, library, character.only = TRUE))


# Suppress all warnings globally
options(warn = -1)


# Set knitr chunk options
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE
)


# mAb chain names
HC_names <- c("h9C12-Q97A_HC", "h9C12-WT_HC", "Brimab_HC", "Umab_HC", "PGT121_HC", "PGDM1400_HC")
LC_names <- c("h9C12_LC", "Brimab_LC", "Umab_LC", "PGT121_LC", "PGDM1400_LC")


# for plotting
match_ig_type_ordered_levels <- c(
  "h9C12-Q97A_HC", "h9C12-WT_HC", "Brimab_HC", "Umab_HC", "PGT121_HC", "PGDM1400_HC",
  "h9C12_LC", "Brimab_LC", "Umab_LC", "PGT121_LC", "PGDM1400_LC"
)

