plot_heatmap_per_mab_per_sample <- function(df, name = "", sample_info, concentrations) {
  
  df_HC <- df[, HC_names]
  df_LC <- df[, c("h9C12_LC", LC_names)]
  
  sample_info$has_blood[sample_info$has_blood] <- 1000
  sample_info$has_blood[c(69, 70)] <- 50000
  
  hmap <- Heatmap(concentratoins, 
                 cell_fun = function(j, i, x, y, width, height, fill) { 
                   if((df_HC[i, j] > 0)  & (df_LC[i, j] > 0)) {
                     grid.text(sprintf("H:%.f L:%.f", df_HC[i, j], df_LC[i, j]), x, y, gp = gpar(fontsize = 8))
                   } else if (df_HC[i, j] > 0) {
                     grid.text(sprintf("H:%.f", df_HC[i, j]), x, y, gp = gpar(fontsize = 8))
                   } else if (df_LC[i, j] > 0) {
                     grid.text(sprintf("L:%.f", df_LC[i, j]), x, y, gp = gpar(fontsize = 8))
                   }
                 } ,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 col=structure(c("white", brewer.pal(4,"Blues")), names = c("0", "1", "10", "100", "1000")), 
                 row_names_gp = grid::gpar(fontsize = 8),
                 column_title = name,   
                 name = "Concentration", 
                 left_annotation = rowAnnotation(
                   GingisKHAN = sample_info$has_GingisKHAN, 
                   blood = sample_info$has_blood, 
                   col = list(blood = c("0" = "white", "1000" = "#ff0000", "50000" = "#9b0000"), 
                              GingisKHAN = c("TRUE" = "black", "FALSE" = "white"))))
  
  return(hmap)
  
}








