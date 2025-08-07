combine_peptides_mq <- function(input_folder) {
  peptide_files <- list.files(input_folder, pattern = ".txt", full.names = TRUE)
  
  combined <- lapply(peptide_files, function(pfile) {
    df <- readr::read_tsv(pfile, guess_max = 1e5)
    df$Protease <- strsplit(basename(pfile), "_")[[1]][2]
    df <- df[, c("Sequence", "Protease", colnames(df)[grepl("Intensity ", colnames(df), fixed = TRUE)])]
    df <- reshape2::melt(df, id.vars = c("Sequence", "Protease"))
    colnames(df) <- c("Sequence", "Protease", "Rawfilenumber", "intensity")
    df$run <- strsplit(basename(pfile), "_")[[1]][1]
    df$Rawfilenumber <- stringr::str_remove(df$Rawfilenumber, "Intensity ")
    df[df$intensity > 0, ]
  }) |> dplyr::bind_rows()
  
  combined$tool <- "MaxQuant"
  combined$Protease <- dplyr::recode(combined$Protease, "ct" = "Ct", "ct+tryp" = "Ct+Tryp", "tryp" = "Tryp", "aspn" = "AspN")
  
  return(combined)
}
