
source(here::here("src", "plotting", "color_schemes.R"))


plot_coverage_for_merged_runs <- function(df) {
  
  plot_coverage_for_merged_runs_helper <- function(df, plot_name) {
    
    df_CI <- df %>%
      group_by(n_runs, match_ig_type, concentration) %>%
      summarise(
        {
          if (length(unique(coverage)) == 1) {
            median  <- unique(coverage)
            lwr.ci  <- median
            upr.ci  <- median
          } else {
            method <- if (n() < 10) "boot" else "exact"
            ci     <- MedianCI(coverage, sides = "two.sided", method = method)
            median <- ci[1]
            lwr.ci <- ci[2]
            upr.ci <- ci[3]
          }
          tibble(median = median, lwr.ci = lwr.ci, upr.ci = upr.ci)
        },
        .groups = "drop"
      ) %>%
      mutate(
        concentration = as.factor(concentration), 
        match_ig_type = factor(match_ig_type, levels = match_ig_type_ordered_levels)
      ) %>%
      arrange(match_ig_type)
    
    
    kw_test_signif_keys <- df %>% 
      group_by(concentration, match_ig_type) %>% 
      kruskal_test(coverage ~ n_runs) %>%
      mutate(p_adj = p.adjust(p, method = "BH")) %>%
      filter(p_adj < 0.05) %>%
      group_by(concentration, match_ig_type) %>% 
      group_keys()
    
    facet_max <- df_CI %>%
      mutate(concentration = as.character(concentration)) %>%
      group_by(concentration, match_ig_type) %>%
      summarise(max_y = max(upr.ci, na.rm = TRUE), .groups = "drop")
    
    wilcox_test_pvals <- df %>% 
      semi_join(kw_test_signif_keys, join_by(concentration, match_ig_type)) %>%
      group_by(concentration, match_ig_type, n_runs) %>%
      filter(length(unique(coverage)) > 1 & n() > 1) %>%  
      ungroup() %>%
      group_by(concentration, match_ig_type) %>%
      group_modify(~ {
        runs_present <- unique(.x$n_runs)
        comps <- list(c(1,4), c(2,4), c(3,4))
        comps <- comps[sapply(comps, function(c) all(c %in% runs_present))]
        if (length(comps) > 0) {
          wilcox_test(.x, coverage ~ n_runs, comparisons = comps, p.adjust.method = "BH")
        } else {
          tibble()
        }
       }) %>%
      ungroup %>%
      mutate(concentration = as.character(concentration)) %>%
      left_join(facet_max, by = c("concentration", "match_ig_type")) %>%
      group_by(concentration, match_ig_type) %>%
      mutate(
        xmin = as.numeric(group1),
        xmax = as.numeric(group2),
        y.position = max_y + 5 + (row_number() - 1) * 8
      ) %>%
      ungroup() %>%
      mutate(
        match_ig_type = factor(match_ig_type, levels = match_ig_type_ordered_levels)
      ) %>%
      arrange(match_ig_type)
 
    
    p <- ggplot(df_CI, aes(x = n_runs, y = median)) + 
      geom_ribbon(aes(ymax = upr.ci, ymin = lwr.ci, fill = concentration), alpha = 0.25) +
      geom_line(aes(color = concentration), linewidth = 1) +
      geom_point(aes(color = concentration), size = 2) +
      stat_pvalue_manual(data = wilcox_test_pvals, hide.ns = TRUE, label = "p.adj.signif", label.size = 4) +
      geom_text(aes(x = n_runs, y = 2, label = round(median, 1)), size = 2) +
      facet_grid(match_ig_type ~ concentration) + 
      scale_fill_manual(values = concentration_cols) + 
      scale_color_manual(values = concentration_cols) +
      scale_x_continuous(limits = c(0.75, 4.25)) +
      scale_y_continuous(limits = c(0, 112), breaks = c(0, 25, 50, 75, 100)) +
      labs(x = "# of merged experimental replicates", y = "Median coverage, %") +
      ggtitle(plot_name) +
      theme_bw() +
      theme(legend.position = "bottom")
    
    return(p)
  }
  
  plot_vdj <- plot_coverage_for_merged_runs_helper(df[df$peptide_type == "vdj", ], "VDJ coverage")
  plot_cdr3 <- plot_coverage_for_merged_runs_helper(df[df$peptide_type == "cdr3", ], "CDR3 coverage")
  
  return(plot_vdj + plot_cdr3)
  
}