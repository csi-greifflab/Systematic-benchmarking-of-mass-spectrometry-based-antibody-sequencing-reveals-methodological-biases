
plot_coverage_single_vs_multiple_spikeins <- function(df) {
  
  plot_coverage_single_vs_multiple_spikeins_helper <- function(df, plot_name) {
    
    df <- df %>%
      mutate(match_ig_type = factor(match_ig_type, levels = match_ig_type_ordered_levels), 
             number_of_spike_ins = factor(number_of_spike_ins, levels = c("single", "multiple")))
    
    med <- df %>%
      group_by(match_ig_type, concentration, number_of_spike_ins) %>%
      summarise(median = median(coverage))
    
    p <- ggplot(df, aes(x = number_of_spike_ins, y = coverage, color = number_of_spike_ins)) +
      geom_beeswarm(pch = 1, size = 1, cex = 3) +
      geom_line(aes(group = run), color = "grey50") +
      scale_color_manual(values = c("#008080","#252525")) +
      geom_text(data = med, aes(y = -10, label = round(median, 2), color = number_of_spike_ins), size = 2) +
      facet_grid(match_ig_type ~ concentration) + 
      labs(x = "# of spike-in mAbs", y = "coverage, %") +
      ggtitle(plot_name) +
      theme_bw() +
      theme(legend.position = "bottom")

    return(p)
  }
  
  plot_vdj <- plot_coverage_single_vs_multiple_spikeins_helper(df[df$peptide_type == "vdj", ], "VDJ coverage")
  plot_cdr3 <- plot_coverage_single_vs_multiple_spikeins_helper(df[df$peptide_type == "cdr3", ], "CDR3 coverage")
  
  return(plot_vdj + plot_cdr3)
  
}