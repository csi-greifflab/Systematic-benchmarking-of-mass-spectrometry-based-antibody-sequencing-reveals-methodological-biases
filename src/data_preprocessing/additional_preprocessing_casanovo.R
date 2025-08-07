additional_preprocessing_casanovo <- function(df, min_search_engine_score=0.8, min_peptide_length=7) {
  # Filter peptides which casanovo_engine_score smaller than 0.8
  df <- df[!is.na(df$search_engine_score) & df$search_engine_score > min_search_engine_score, ]
  
  # Remove PTMs, e.g. YLTSM+15.995ASR will be replaced with YLTSMASR
  df$Sequence <- stringr::str_replace_all(df$Sequence, "[0-9.+-]", "")
  
  # Filter peptides shorter than 6 AA
  df <- df[sapply(df$Sequence, nchar) >= min_peptide_length, ]
  
  # Merge identical peptides and keep the highest search_engine_score for merged peptides
  df <- df%>%
    dplyr::group_by(Sequence, Rawfilenumber, run, Protease) %>% # for each unique peptide in a filename, run, and enzyme
    dplyr::arrange(dplyr::desc(search_engine_score)) %>% # sort by highest search_engine score
    dplyr::slice(1) %>% # take the peptide with the highest search_engine score
    dplyr::ungroup()

  df <- df[, c("Sequence", "Protease", "Rawfilenumber", "search_engine_score", "run", "tool")]
  
  return(df)
}


