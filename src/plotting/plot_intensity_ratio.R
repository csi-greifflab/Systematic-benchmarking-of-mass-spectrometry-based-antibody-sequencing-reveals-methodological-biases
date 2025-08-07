
plot_intensity_ratio <- function(df) {
  
  plot_intensity_ratio_helper <- function(df, plot_name) {
    
    df <- df %>%
      mutate(Protease = factor(Protease, levels = c("Tryp", "Ct", "Ct+Tryp")), 
             match_ig_type = factor(match_ig_type, levels = match_ig_type_ordered_levels))
    
    med <- df %>% 
      group_by(match_ig_type, Protease, run, log_conc_ratio) %>%
      summarise(median_log = median(log_intensity_ratio), 
                median = median(intensity_ratio),
                n_points = n())
    
    coefs <- med %>%
      group_by(Protease, run, match_ig_type) %>%
      filter(n_distinct(log_conc_ratio) > 1) %>%
      do({
        model <- lm(median_log ~ log_conc_ratio, data = ., weights = n_points)
        data.frame(
          intercept = coef(model)[1],
          slope = coef(model)[2],
          r2 = summary(model)$r.squared
        )
      }) %>%
      ungroup() %>%
      mutate(intercept_eq = paste0("Intercept: ", format(round(intercept,2), digits = 2)), 
             slope_eq = paste0("Slope: ", format(round(slope,2), digits = 2)), 
             r2_eq = paste0("R^2: ", format(round(r2, 2), digits = 2)), 
             intercept_ground_truth = 0, 
             slope_ground_truth = 1)
    
    p <- ggplot(df, aes(x = log_conc_ratio, y = log_intensity_ratio)) + 
      geom_point(aes(fill = run), position = position_jitterdodge(), color = "black", shape = 21, size = 0.1, stroke = 0.03) +
      geom_boxplot(aes(fill = run, group = interaction(run, log_conc_ratio)), position = position_dodge(width = 0.6), width = 0.4, lwd = 0.3, outlier.shape = NA) +
      geom_text(data = med, aes(x = log_conc_ratio, y = 4.2, label = round(median, 1)), color = "black", size = 2) +
      geom_text(data = coefs, aes(-0.75, -3.25, label = intercept_eq), size = 2.25, color = "black", hjust = 0) +
      geom_text(data = coefs, aes(-0.75, -4, label = slope_eq), size = 2.25, color = "black", hjust = 0) +
      geom_text(data = coefs, aes(-0.75, -4.75, label = r2_eq), parse = TRUE, size = 2.25, color = "black", hjust = 0) +
      geom_abline(data = coefs, aes(intercept = intercept, slope = slope), color = "#DE583E", linewidth = 0.75, linetype = "dashed") +
      geom_abline(data = coefs, aes(intercept = intercept_ground_truth, slope = slope_ground_truth), color = "black", linewidth = 0.25) + 
      facet_grid(match_ig_type ~ Protease + run) + 
      scale_fill_manual(values = c("#969696", "#bdbdbd", "#d9d9d9", "#ffffff")) +
      scale_x_continuous(breaks = 0:3, labels = c(1, 10, 100, 1000)) +
      scale_y_continuous(breaks = -3:4, labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) + 
      theme_bw() +
      theme(legend.position = "bottom", 
            axis.text.x = element_text(angle = 90, vjust = 0.45), 
            axis.text.y = element_text(size = 7),
            panel.spacing = unit(0.1, "lines")) +
      ggtitle(plot_name) +
      labs(x = "mAb input ratio", y = "MS signal intensity ratio")
    
    return(p)
  }
  
  plot_blood <- plot_intensity_ratio_helper(df[df$has_blood, ], "Intensity ratio in blood samples")
  plot_nonblood <- plot_intensity_ratio_helper(df[!df$has_blood, ], "Intensity ratio in non-blood samples")
  
  return(plot_nonblood / plot_blood )
  
}