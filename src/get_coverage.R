get_coverage_per_chain <- function(peptides_vector, sequence_name, mabs_seq, 
                                   region = c("cdr3", "vdj"), 
                                   return_type = c("coverage", "counts")) {
  
  region <- match.arg(region)
  return_type <- match.arg(return_type)
  
  V_ref <- unlist(mabs_seq[mabs_seq$sequence_name == sequence_name, "amino_acid_sequence_Vregion"])
  VC_ref <- unlist(mabs_seq[mabs_seq$sequence_name == sequence_name, "amino_acid_sequence_full"])
  CDR3_ref <- unlist(mabs_seq[mabs_seq$sequence_name == sequence_name, "amino_acid_sequence_cdr3"])
  
  ref <- if (region == "cdr3") CDR3_ref else V_ref
  
  alignment_pos <- do.call(rbind, str_locate_all(string = VC_ref, pattern = peptides_vector)) |> as.data.frame()
  if (nrow(alignment_pos) == 0) return(0)
  
  counts <- integer(nchar(VC_ref))
  
  for (i in seq_len(nrow(alignment_pos))) {
    counts[alignment_pos[i, "start"]:alignment_pos[i, "end"]] <- 
      counts[alignment_pos[i, "start"]:alignment_pos[i, "end"]] + 1
  }
  
  if (return_type == "counts") {
    return(counts[1:nchar(V_ref)])
  }
  
  ref_location <- str_locate(string = VC_ref, pattern = ref)
  ref_length <- ref_location[2] - ref_location[1] + 1
  ref_positions <- ref_location[1]:ref_location[2]
  
  coverage <- sum(counts[ref_positions] > 0) / ref_length * 100
  
  return(coverage)
}


get_mabs_coverage <- function(df, group_vars, mabs_seq, complete_args = NULL) {

  chains <- c(
    "h9C12-Q97A_HC", "h9C12-WT_HC", "Brimab_HC", "Umab_HC", "PGT121_HC", "PGDM1400_HC",
    "h9C12_LC", "Brimab_LC", "Umab_LC", "PGT121_LC", "PGDM1400_LC"
  )
  
  df_cov <- df %>%
    group_by(across(all_of(group_vars)), match_ig_type) %>%
    summarise(peptides = list(unique(Sequence)), .groups = "drop") %>%
    rowwise() %>%
    mutate(
      vdj  = get_coverage_per_chain(peptides, match_ig_type, mabs_seq, region = "vdj"),
      cdr3 = get_coverage_per_chain(peptides, match_ig_type, mabs_seq, region = "cdr3")
    ) %>%
    ungroup() %>%
    select(-peptides) %>%
    pivot_longer(
      cols = c(vdj, cdr3),
      names_to = "peptide_type",
      values_to = "coverage"
    ) %>%
    complete(
      !!!syms(group_vars),
      match_ig_type = chains,
      peptide_type = c("vdj", "cdr3"),
      fill = list(coverage = 0)
    )
  
  # Fill missing combinations with 0
  if (!is.null(complete_args)) {
    df_cov <- df_cov %>%
      complete(!!!complete_args, 
               match_ig_type = chains,
               fill = list(coverage = 0))
  }

  
  return(df_cov)
}



get_coverage_per_ab_per_protease <- function(df, mabs_seq) {

  df_cov <- get_mabs_coverage(df, group_vars = "Protease", mabs_seq = mabs_seq)
  
  df_cov_all <- get_mabs_coverage(df, group_vars = NULL, mabs_seq = mabs_seq) %>%
    mutate(Protease = "All")
  
  bind_rows(df_cov, df_cov_all) %>%
    pivot_wider(
      names_from = match_ig_type,
      values_from = coverage
    ) %>%
    replace(is.na(.), 0) %>%
    mutate(Protease = factor(Protease, levels = c("All", "Ct+Tryp", "Tryp", "Ct", "AspN"))) %>%
    arrange(Protease)
}



get_coverage_per_ab_per_sample <- function(df, mabs_seq) {
  
  get_mabs_coverage(df, group_vars = "Sample", mabs_seq = mabs_seq,
    complete_args = list(Sample = 1:70, peptide_type = c("vdj", "cdr3"))) %>%
    pivot_wider(
      names_from = match_ig_type,
      values_from = coverage
    ) %>%
    replace(is.na(.), 0) %>%
    arrange(Sample)
}



get_coverage_per_concentration_per_merged_runs <- function(df, mabs_seq, sample_info) {
  
  replicate_subsets <- unlist(lapply(1:4, combn, x = paste0("run", 1:4), simplify = FALSE), recursive = FALSE)
  names(replicate_subsets) <- sapply(replicate_subsets, paste, collapse = "_")
  
  sample_info <- sample_info %>%
    select(c("Sample", str_subset(colnames(sample_info), "concentration"))) %>%
    mutate(`concentration_h9C12` = pmax(`concentration_h9C12-Q97A`, `concentration_h9C12-WT`)) %>%
    pivot_longer(
      cols = starts_with("concentration_"),
      names_to = "mab",
      values_to = "concentration"
    ) %>%
    mutate(
      mab = str_remove(mab, "^concentration_")
    ) %>%
    filter(concentration > 0)
   
  coverage <- lapply(names(replicate_subsets), function(name) {
    get_mabs_coverage(
      df[df$run %in% replicate_subsets[[name]], ], 
      group_vars = c("concentration", "Sample"),
      mabs_seq = mabs_seq, 
      complete_args = list(
        concentration = c(1, 10, 100, 1000), 
        Sample = 1:70,
        peptide_type = c("vdj", "cdr3"))
    ) %>% 
      mutate(subset = name)  
    }) %>%
    bind_rows() %>%
    mutate(n_runs = str_count(subset, "run"), 
           mab = str_remove(match_ig_type, "_[HL]C$")) %>%
    # Keep only those mabs present in the sample
    merge(sample_info, by = c("Sample", "mab", "concentration"))

  return(coverage)
}



get_coverage_per_concentration_per_run <- function(df, mabs_seq, sample_info) {
  
  sample_info <- sample_info %>%
    select(c("Sample", "has_blood", "n_abs", str_subset(colnames(sample_info), "concentration"))) %>%
    mutate(`concentration_h9C12` = pmax(`concentration_h9C12-Q97A`, `concentration_h9C12-WT`)) %>%
    pivot_longer(
      cols = starts_with("concentration_"),
      names_to = "mab",
      values_to = "concentration"
    ) %>%
    mutate(
      mab = str_remove(mab, "^concentration_")
    ) %>%
    filter(concentration > 0)
  
  coverage <- get_mabs_coverage(
    df, 
    group_vars = c("concentration", "Sample", "run"),
    mabs_seq = mabs_seq, 
    complete_args = list(
      concentration = c(1, 10, 100, 1000), 
      Sample = 1:70,
      run = c("run1", "run2", "run3", "run4"),
      peptide_type = c("vdj", "cdr3"))
    ) %>%
    mutate(mab = str_remove(match_ig_type, "_[HL]C$")) %>%
    merge(sample_info, by = c("Sample", "mab", "concentration"))
  
  return(coverage)
  
}


get_sequence_coverage_per_tool <- function(df, mabs_seq) {
  
  df %>%
    group_by(tool, match_ig_type) %>%
    summarise(peptides = list(unique(Sequence)), .groups = "drop") %>%
    rowwise() %>%
    mutate(counts = list(get_coverage_per_chain(
      peptides_vector = peptides,
      sequence_name = match_ig_type,
      mabs_seq = mabs_seq,
      region = "vdj",
      return_type = "counts"
    ))) %>%
    ungroup() %>%
    unnest_longer(col = counts) %>%     
    group_by(tool, match_ig_type) %>%
    mutate(position = row_number()) %>%
    ungroup()
  
}






