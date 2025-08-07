source(here::here("src", "align_peptides_to_references.R"))


keep_peptides_mapped_to_mabs <- function(df, mabs_seq, min_variable_region_overlap=5) {
  
  map_to_one_mab_chain <- function(df, i) {
    v_seq <- mabs_seq$amino_acid_sequence_Vregion[i]
    full_seq <- mabs_seq$amino_acid_sequence_full[i]
    mab_name <- mabs_seq$sequence_name[i]
    
    # Locate variable region in full-length mab sequence
    v_pos <- stringr::str_locate(full_seq, v_seq)
    # Locate peptides in full-length mab sequence
    peptide_pos <- stringr::str_locate(full_seq, df$Sequence)
    
    # Overlap with variable region
    df$overlap <- pmax(
      0,
      pmin(v_pos[1, 2], peptide_pos[, 2]) - pmax(v_pos[1, 1], peptide_pos[, 1]) + 1
    )

    matched <- df[which(df$overlap > min_variable_region_overlap), ]
    if (nrow(matched) == 0) return(NULL)
    
    # Trim peptides exceeding variable region
    idx <- nchar(matched$Sequence) > matched$overlap
    matched$Sequence[idx] <- substring(matched$Sequence[idx], 1, matched$overlap[idx])
    matched$match_ig_type <- mab_name
    return(matched)
  }
  
  results <- lapply(seq_len(nrow(mabs_seq)), function(i) map_to_one_mab_chain(df, i)) %>%
    bind_rows() %>%
    select(!overlap)
    
  return(results)
}



annotate_contaminations <- function(df, sample_info, rawfiles_info) {
  original_colnames <- colnames(df)
  rawfiles_info$sample_id <- paste(rawfiles_info$Rawfilenumber, rawfiles_info$Protease, rawfiles_info$run, sep="_")
  df <- merge(df, rawfiles_info[, c("Rawfilenumber", "Sample", "replicate", "sample_id")], by="Rawfilenumber") 
  df <- merge(df, sample_info, by="Sample", all.x = T)
  
  df <- df %>%
    mutate(concentration = case_when(
      match_ig_type %in% c("h9C12_LC")                   ~ pmax(`concentration_h9C12-Q97A`, `concentration_h9C12-WT`, na.rm = TRUE),
      match_ig_type %in% c("h9C12-Q97A_HC")              ~ `concentration_h9C12-Q97A`,
      match_ig_type %in% c("h9C12-WT_HC")                ~ `concentration_h9C12-WT`,
      match_ig_type %in% c("Brimab_HC", "Brimab_LC")     ~ concentration_Brimab,
      match_ig_type %in% c("Umab_HC", "Umab_LC")         ~ concentration_Umab,
      match_ig_type %in% c("PGT121_HC", "PGT121_LC")     ~ concentration_PGT121,
      match_ig_type %in% c("PGDM1400_HC", "PGDM1400_LC") ~ concentration_PGDM1400,
      TRUE ~ NA_real_)) %>%
    rowwise() %>%
    mutate(is_contamination = !(isTRUE(get(match_ig_type)))
    ) %>%
    ungroup()
  
  return(df)
}



annotate_cdr3_peptides <- function(df, mabs_seq, cdr3_min_overlap=3) {
  df_cdr3 <- df
  df_cdr3 <- merge(x = df_cdr3, y = mabs_seq, by.x = "match_ig_type", by.y = "sequence_name",  all.x = TRUE, sort = F)
  res_seq <- str_locate(df_cdr3$amino_acid_sequence_full, as.character(df_cdr3$Sequence))
  res_cdr <- str_locate(df_cdr3$amino_acid_sequence_full, as.character(df_cdr3$amino_acid_sequence_cdr3))
  df_cdr3$cdr3_overlap <- mapply(max, mapply(min, res_seq[,2], res_cdr[,2]) - mapply(max, res_seq[,1], res_cdr[,1]), 0)
  df_cdr3$is_cdr3_related <- df_cdr3$cdr3_overlap > cdr3_min_overlap
  return(df_cdr3[, c(colnames(df), "is_cdr3_related", "cdr3_overlap")])
}


annotate_imgt_peptides <- function(df, imgt_genes) {
  imgt_vector <- setNames(imgt_genes[[2]], imgt_genes[[1]])
  blood_contam_peptides <- unique(df[df$has_blood & df$is_contamination, "Sequence"])
  blood_alignment <- align_to_ref(reference = imgt_vector, 
                                  peptides = unlist(blood_contam_peptides))
  df$is_imgt_peptide <- df$has_blood & df$is_contamination & df$Sequence %in% blood_alignment$sequence
  return(df)
}


annotate_shared_peptides <- function(df, shared_mab_subseq) {
  
  rows_to_remove <- c()
  for (i in 1:nrow(shared_mab_subseq)) {
    ab1 <- shared_mab_subseq$Ab1[i]
    ab2 <- shared_mab_subseq$Ab2[i]
    peptide_candidates <- unique(df[(df$match_ig_type == ab1 & df[, ab2] & !df[, ab1]) | 
                                    (df$match_ig_type == ab2 & df[, ab1] & !df[, ab2]) , "Sequence"])
    ref <- unlist(str_split(shared_mab_subseq$sequence[i], ";"))
    names(ref) <- ref
    aligned_peptides <- align_to_ref(reference = ref, peptides = unlist(peptide_candidates))
    if (nrow(aligned_peptides) > 0) {
      to_remove <- which(((df$match_ig_type == ab1 & df[, ab2] & !df[, ab1]) | 
                           (df$match_ig_type == ab2 & df[, ab1] & !df[, ab2])) & (df$Sequence %in% aligned_peptides$sequence))
      rows_to_remove <- c(rows_to_remove, to_remove)
    }
  }
  
  df$is_shared_peptide <- FALSE
  df$is_shared_peptide[rows_to_remove] <- TRUE
  
  return(df)
}


annotate_carryover_peptides <- function(df, rawfiles_info) {
  
  assign_groups <- function(fnames) {
    group_id = 0
    fnames$group <- NA
    j <- 1
    while (j <= nrow(fnames)) {
      while ((fnames$Sample[j] != "blank") & (j <= nrow(fnames))) {
        fnames$group[j] <- group_id
        j <- j+1
      }
      group_id <- group_id+1
      while ((fnames$Sample[j] == "blank") & (j <= nrow(fnames))) {
        fnames$group[j] <- group_id
        j <- j+1
      }
    }
    return(fnames)
  }
  
  file_names_splitted <- split(rawfiles_info, list(rawfiles_info$run, rawfiles_info$Protease))
  rows_to_remove <- c()
  
  for(i in 1:length(file_names_splitted)) {
    fnames_by_group <- file_names_splitted[[i]] %>%
      assign_groups() %>%
      split(.$group)
    
    for (k in 1:length(fnames_by_group)) {
      fnames_group <- fnames_by_group[[k]]
      blank_ids <- fnames_group[fnames_group$Sample == "blank", ]$Rawfilenumber
      sample_ids <- fnames_group[fnames_group$Sample != "blank", ]$Rawfilenumber
      blank_peptides <- unique(df[df$Rawfilenumber %in% blank_ids, ]$Sequence)
      sample_contam_peptides <- unique(df[df$Rawfilenumber %in% sample_ids & df$is_contamination, ]$Sequence)

      if (length(blank_peptides) > 0 & length(sample_contam_peptides > 0)) {
        names(blank_peptides) <- blank_peptides
        aligned_peptides <- align_to_ref(reference = blank_peptides,
                                         peptides = sample_contam_peptides)
        to_remove <- which((df$run == fnames_group$run[1]) &
                             (df$Rawfilenumber %in% sample_ids) &
                             (df$Sequence %in% aligned_peptides$sequence) &
                             df$is_contamination)
        rows_to_remove <- c(rows_to_remove, to_remove)
      }
    }
  }
  
  df$is_carryover_peptide <- FALSE
  df$is_carryover_peptide[rows_to_remove] <- TRUE
  
  return(df)
  
}











