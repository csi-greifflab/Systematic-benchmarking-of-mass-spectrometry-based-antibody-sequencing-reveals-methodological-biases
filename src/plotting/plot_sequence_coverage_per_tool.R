source(here("src", "plotting", "color_schemes.R"))

plot_sequence_coverage_per_tool <- function(df, mabs_seq){
  
  df <- df %>%
    mutate(match_ig_type = factor(match_ig_type, levels = match_ig_type_ordered_levels), 
           tool = factor(tool, levels = c("MaxQuant", "MSFragger", "Casanovo")))
  
  cdr3_pos <- str_locate(string = mabs_seq$amino_acid_sequence_full, 
                         pattern = mabs_seq$amino_acid_sequence_cdr3) %>%
    as.data.frame() %>%
    bind_cols(match_ig_type = mabs_seq$sequence_name) %>%
    mutate(match_ig_type = factor(match_ig_type)) %>%
    filter(match_ig_type %in% df$match_ig_type)
  
  p <- ggplot(df, aes(x = position, y=counts, fill = tool)) +
    geom_rect(data = cdr3_pos, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
              fill = "yellow", alpha = 0.5, inherit.aes = FALSE) +
    geom_col() +
    facet_grid(match_ig_type ~ .) +
    scale_fill_manual(values=tools_cols) +
    scale_x_continuous(breaks = seq(0, max(df$position), by = 20)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", 
          strip.text.y = element_text(angle = 0)) +
    labs(x = "V(D)J region position", y = "Count")
  
  return(p)
}

