plot_casanovo_threshold <- function(df) {
  p <- ggplot(df, aes(x = threshold, y = count)) +
    geom_vline(xintercept = 0.8, color = "red") +
    geom_line(aes(group = metric, linetype = metric), color = "#6a3e37", linewidth = 0.8) +
    geom_point(size = 1) +
    geom_text(aes(label = count), vjust = -0.45, hjust = -0.25, size = 2) +
    scale_y_log10(limits = c(14000, 3000000)) +
    scale_x_continuous(limits = c(0.5, 0.94)) +
    scale_linetype_manual(
      values = c("total_peptides" = "solid", "mab_peptides" = "dotted"),
      labels = c("Total Casanovo peptides", "mAb-related Casanovo peptides")
    ) +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom") +
    labs(x = "Casanovo prediction score threshold", y = "# of peptides", linetype = "Peptide type")
  
  return(p)
}