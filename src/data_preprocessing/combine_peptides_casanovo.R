combine_peptides_casanovo <- function(input_folder) {
  peptide_files <- list.files(input_folder, pattern = ".mztab", full.names = TRUE, recursive = TRUE)
  
  combined <- lapply(peptide_files, function(pfile) {
    df <- readr::read_tsv(pfile, skip = 58) # skip mztab header
    df <- df[, c("sequence", "search_engine_score[1]")]
    colnames(df) <- c("Sequence", "search_engine_score")
    df$Rawfilenumber <- stringr::str_remove(basename(pfile), ".mztab")
    df$run <- stringr::str_extract(pfile, "run[0-9]")
    df$Protease <- stringr::str_extract(basename(dirname(pfile)), "(ct\\+tryp|ct|tryp|aspn)")
    df <- df[, c("Sequence", "Protease", "Rawfilenumber", "search_engine_score", "run")]
    
  }) |> dplyr::bind_rows()
  
  combined$tool <- "Casanovo"
  combined$Protease <- dplyr::recode(combined$Protease, "ct" = "Ct", "ct+tryp" = "Ct+Tryp", "tryp" = "Tryp", "aspn" = "AspN")
  return(combined)
}


