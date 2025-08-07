source(here::here("src", "plotting", "color_schemes.R"))

plot_total_number_of_peptides_per_mab <- function(n_peptides) {
  
  ggplot(n_peptides, aes(x = match_ig_type,  y = n_peptides, color = Protease, label = n_peptides)) +
    geom_point(size=2) +
    geom_text(hjust=1.8, size=2.5) +
    facet_grid(peptide_type ~ ., scales = "free_y") +
    scale_colour_manual(values = enzyme_cols) +
    theme_minimal() +
    theme(
      text = element_text(size = 8),
      legend.position = "bottom",
      panel.spacing = unit(1, "lines"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

}