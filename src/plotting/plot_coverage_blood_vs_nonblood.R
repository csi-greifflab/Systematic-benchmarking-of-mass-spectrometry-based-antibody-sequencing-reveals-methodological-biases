
plot_coverage_blood_vs_nonblood <- function(df) {
  
  
  plot_coverage_blood_vs_nonblood_helper <- function(df, plot_name) {
    
    df <- df %>%
      mutate(match_ig_type = factor(match_ig_type, levels = match_ig_type_ordered_levels))
    
    
    med <- df %>%
      group_by(match_ig_type, concentration, has_blood) %>%
      summarise(median = median(coverage))
    
    
    mood_pvals <- df %>%
      group_by(concentration, match_ig_type) %>%
      group_modify(~ {
        
        if (length(unique(.x$has_blood)) < 2) 
          return(tibble(p = NA_real_))
        
        if (n_distinct(.x$coverage) < 2) 
          return(tibble(p = NA_real_))
        
        median_val <- median(.x$coverage)
        above_median <- .x$coverage > median_val
        cont_table <- table(above_median, .x$has_blood)
        if (nrow(cont_table) < 2 || ncol(cont_table) < 2) 
          return(tibble(p = NA_real_))
        
        test_result <- mood.medtest(coverage ~ has_blood, data = .x)
        tibble(p = test_result$p.value)
        
      }) %>%
      ungroup() %>%
      mutate(
        p.adj = p.adjust(p, method = "BH"),
        p.adj.signif = symnum(
          p.adj,
          corr = FALSE,
          na = FALSE,
          cutpoints = c(0, 0.001, 0.01, 0.05, 1),
          symbols = c("***", "**", "*", "ns")
        ),
        group1 = "TRUE",
        group2 = "FALSE",
        xmin = 1,
        xmax = 2,
        y.position = 105
      ) %>%
      filter(!is.na(p.adj))
    
    
    p <- ggplot(df, aes(x = has_blood, y = coverage, color = has_blood)) +
      geom_boxplot(aes(color = has_blood), width = 0.5, outlier.shape = NA, alpha = 0.5) +
      geom_beeswarm(aes(color = has_blood), pch = 1, size = 0.8, cex = 3, corral = "wrap", corral.width = 0.9) +
      stat_pvalue_manual(data = mood_pvals, hide.ns = TRUE, label = "p.adj.signif", label.size = 4) +
      geom_text(data = med, aes(y = -10, label = round(median, 2), color = has_blood), size = 2) +
      scale_color_manual(values = c("red3", "black")) +
      scale_y_continuous(breaks = seq (0, 100, by = 25), limits = c(-15, 115)) +
      facet_grid(match_ig_type ~ concentration) +
      ggtitle(plot_name) +
      theme_bw() +
      theme(legend.position = "bottom")
    
    return(p)
  }
  
  plot_vdj <- plot_coverage_blood_vs_nonblood_helper(df[df$peptide_type == "vdj", ], "VDJ coverage")
  plot_cdr3 <- plot_coverage_blood_vs_nonblood_helper(df[df$peptide_type == "cdr3", ], "CDR3 coverage")
  
  return(plot_vdj + plot_cdr3)
  
}