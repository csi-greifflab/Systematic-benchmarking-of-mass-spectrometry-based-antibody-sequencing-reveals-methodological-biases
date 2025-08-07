source(here::here("src", "data_preprocessing",  "additional_preprocessing_casanovo.R"))
source(here::here("src", "data_preprocessing",  "combine_peptides_casanovo.R"))

get_n_unique_peptides_per_mab <- function(df) {
  
  df_protease <- df %>%
    group_by(Protease, match_ig_type) %>%
    summarise(
      vdj  = length(unique(Sequence)),
      cdr3 = length(unique(Sequence[is_cdr3_related])),
      .groups = "drop"
    )
  
  df_all <- df %>%
    group_by(match_ig_type) %>%
    summarise(
      vdj  = length(unique(Sequence)),
      cdr3 = length(unique(Sequence[is_cdr3_related])),
      .groups = "drop"
    ) %>%
    mutate(Protease = "All")
  
  df_n_ab <- bind_rows(df_protease, df_all) %>%
    pivot_longer(
      cols = c(vdj, cdr3),
      names_to = "peptide_type",
      values_to = "n_peptides"
    )
  
  df_n_ab$match_ig_type <- factor(
    df_n_ab$match_ig_type,
    levels = c(
      "h9C12-Q97A_HC", "h9C12_LC", "h9C12-WT_HC", 
      "Brimab_HC", "Brimab_LC", "Umab_HC", "Umab_LC",
      "PGT121_HC", "PGT121_LC", "PGDM1400_HC", "PGDM1400_LC"
    )
  )
  
  df_n_ab$peptide_type <- factor(
    df_n_ab$peptide_type,
    levels = c("vdj", "cdr3")
  )
  
  return(df_n_ab)
    
}

get_n_peptides_per_ab_per_sample <- function(df, mode) {
  
  distinct_match_ig <- unique(df$match_ig_type)
  
  df_cov <- df %>%
    group_by(Sample, match_ig_type) %>%
    summarise(vdj = length(unique(Sequence)), 
              cdr3 = length(unique(Sequence[is_cdr3_related]))) %>%
    ungroup() %>%
    pivot_longer(
      cols = c(vdj, cdr3),
      names_to = "peptide_type",
      values_to = "n_peptides"
    ) %>%
    complete(
      Sample = 1:70,
      peptide_type = c("vdj", "cdr3"),
      match_ig_type = c("h9C12-Q97A_HC", "h9C12-WT_HC", "Brimab_HC", "Umab_HC", "PGT121_HC", "PGDM1400_HC", 
                        "h9C12_LC", "Brimab_LC", "Umab_LC", "PGT121_LC", "PGDM1400_LC"),
      fill = list(n_peptides = 0)  
    ) %>%
    pivot_wider(
      names_from = match_ig_type,
      values_from = n_peptides, 
      values_fill = 0
    ) %>%
    arrange(Sample)
  
  return(df_cov)
}


get_n_peptides_per_casanovo_threshold <- function(casanovo_path, thresholds, mabs_seq) {
  
  casanovo_peptides_raw <- combine_peptides_casanovo(casanovo_path) %>%
    additional_preprocessing_casanovo(min_search_engine_score = 0, min_peptide_length = 7) 
  
  do.call(rbind, lapply(thresholds, function(t) {
    df_filtered <- casanovo_peptides_raw %>%
      filter(search_engine_score > t)
    
    n_total <- nrow(df_filtered)
    n_mapped <- nrow(filter(keep_peptides_mapped_to_mabs(df_filtered, mabs_seq)))
    
    data.frame(
      threshold = t,
      total_peptides = n_total,
      mab_peptides = n_mapped
    )
  })) %>%
    pivot_longer(
      cols = c(total_peptides, mab_peptides),
      names_to = "metric",
      values_to = "count"
    )
}



