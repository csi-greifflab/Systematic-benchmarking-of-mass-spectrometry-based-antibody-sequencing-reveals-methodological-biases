source(here("src", "plotting", "color_schemes.R"))

plot_total_coverage_per_mab <- function(df) {
  
  get_heatmap <- function(coverage, title) {
    
    m <- as.matrix(coverage[, -(1:2)])
    rownames(m) <- unlist(coverage[, 1])
    
    hmap <- Heatmap(m, 
                    cell_fun = function(j, i, x, y, width, height, fill) 
                    {grid.text(sprintf("%.2f", m[i, j]), x, y, gp = gpar(fontsize = 9))}, 
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    rect_gp = gpar(col = "white", lwd = 4), 
                    heatmap_legend_param = list(direction = "horizontal"),
                    col = colorRamp2(c(0, 70, 100), c("#6EB5FF", "white", "#FFABAB")), 
                    right_annotation = rowAnnotation(enzyme = rownames(m), col = list(enzyme = enzyme_cols)), 
                    row_title = title)
    
    return(hmap)
  }
  
  h_vdj <- get_heatmap(df[df$peptide_type == "vdj", ], "vdj_coverage")
  h_cdr3 <- get_heatmap(df[df$peptide_type == "cdr3", ], "cdr3_coverage")
  
  return(h_vdj %v% h_cdr3)
}