# data preprocessing functoins

preprocess <- function(df, cdr3, cdr3_min_overlap) {
  df <- align_to_ab_sequences(df, cdr3)
  print("Alignment done")
  df <- add_metainfo(df)
  print("Meta info added")
  df <- is_cdr3_related(df, cdr3, cdr3_min_overlap)
  print("Cdr3 peptides info added")
  return(df)
} 


align_to_ab_sequences <- function(df, cdr3) {
  
  # peptides that perfectly match to the variable region
  
  df_res <- data.frame()
  
  for (i in 1:nrow(cdr3)) {
    variable_positions <- str_locate(cdr3$amino_acid_sequence_full[i], cdr3$amino_acid_sequence_Vregion[i])
    constant_match <- str_locate(cdr3$amino_acid_sequence_full[i], df$Sequence)
    df$overlap <- pmax(0, pmin(variable_positions[1, 2], constant_match[, 2]) - 
                         pmax(variable_positions[1, 1], constant_match[, 1]) + 1)
    df_matched <- df[which(df$overlap > 5), ]
    
    # trim peptides by the end of the variable region
    idx <-( nchar(df_matched$Sequence) > df_matched$overlap) 
    df_matched$Sequence[idx] <- substring(df_matched$Sequence[idx], 1, df_matched$overlap[idx])
    df_matched$match_ig_type <- cdr3[i, "sequence_name"]
    df_res <- rbind(df_res, df_matched)
  }
  return(df_res)
}
    
  

add_metainfo <- function(df) {
  merge_rawfiles_description()
  generate_sample_info()
  sample_info <- read_tsv(file.path(metadata_path, "detailed_sample_description.tsv"))
  rawfiles_info <- read_tsv(file.path(metadata_path, "rawfiles_description_with_blanks.tsv"), guess_max = 5e3)
  
  rawfiles_info$sample_id <- paste(rawfiles_info$Rawfilenumber, rawfiles_info$Protease, rawfiles_info$run, sep="_")
  df <- merge(df, rawfiles_info[, c("Rawfilenumber", "Sample", "sample_id")], by="Rawfilenumber", all.x = T) 
  df <- merge(df, sample_info, by="Sample", all.x = T)
  return(df)
}


is_cdr3_related <- function(df, cdr3, cdr3_min_overlap) {
  df_cdr3 <- df
  df_cdr3 <- merge(x = df_cdr3, y = cdr3, by.x = "match_ig_type", by.y = "sequence_name",  all.x = TRUE, sort = F)
  res_seq <- str_locate(df_cdr3$amino_acid_sequence_full, as.character(df_cdr3$Sequence))
  res_cdr <- str_locate(df_cdr3$amino_acid_sequence_full, as.character(df_cdr3$amino_acid_sequence_cdr3))
  df_cdr3$overlap <- mapply(max, mapply(min, res_seq[,2], res_cdr[,2]) - mapply(max, res_seq[,1], res_cdr[,1]), 0)
  df_cdr3$is_cdr3_related <- df_cdr3$overlap > cdr3_min_overlap
  return(df_cdr3)
}



generate_sample_info <- function() {
  d <- data.frame(read_tsv(file.path(metadata_path, "sample_description.tsv")))
  d$Sample <- 1:nrow(d)
  d$has_blood <- str_detect(d$Sample_info, "blood-isolated hIgG1")
  d$has_GingisKHAN <- str_detect(d$Sample_info, "GingisKHAN")
  
  d$Sample_info <- str_replace_all(d$Sample_info, "1 μg", "1000 ng")
  
  long_to_short_name <- c("h9C12-Q97A", "h9C12-WT", "Brimab", "PGDM1400", "PGT121", "Umab")
  names(long_to_short_name) <- c("h9C12 Q97A", "h9C12 WT", "Briakinumab", "PGDM1400", "PGT121", "Ustekinumab")
  
  get_conc <- function(ab, sample_info) {
    res <- rep(0, length(sample_info))
    rgxpr <- regexpr(paste0("\\d+(?= ng ", ab, ")"), sample_info, perl=T)
    res[rgxpr!=-1] <- regmatches(sample_info, rgxpr)
    as.numeric(res)
  }
  
  for (i in 1:length(long_to_short_name)) {
    long_name <- names(long_to_short_name)[i]
    short_name <- long_to_short_name[i]
    d[, short_name] <- str_detect(d$Sample_info, long_name)
    d[, paste0(c(short_name, "HC"), collapse = "_")] <- d[, short_name]
    if (short_name == "h9C12-Q97A" | short_name == "h9C12-WT") {
      d[, "h9C12_LC"] <- str_detect(d$Sample_info, c("h9C12 Q97A", "h9C12 WT"))
    } else {
      d[, paste0(c(short_name, "LC"), collapse = "_")] <- d[, short_name]
    }
    d[, paste0(c("concentration", short_name), collapse = "_")] <- get_conc(long_name, d$Sample_info)
  }
  
  d$n_abs <-  d$'h9C12-Q97A' + d$'h9C12-WT' + d$Brimab + d$Umab + d$PGDM1400 + d$PGT121
  write_tsv(d, file.path(metadata_path, "detailed_sample_description.tsv"))
}

merge_rawfiles_description <- function() {
  df <- list()
  for (run in paste0("run", 1:4)) {
    df[[run]] <- data.frame(import_list(file.path(metadata_path, paste0("rawfiles_description_", run, ".xlsx")), 
                                           setclass = "tbl", rbind = TRUE))
    df[[run]]$run <- run
  }
  df <- do.call("rbind", df)
  df <- df[, colnames(df) != "X_file"]
  write_tsv(df, file.path(metadata_path, "rawfiles_description_with_blanks.tsv"))
  write_tsv(df[df$Sample != "blank", ], file.path(metadata_path, "rawfiles_description_wo_blanks.tsv"))
}


rename_enzymes <- function(df) {
  df$Protease[df$Protease == "ct"] <- "Ct"
  df$Protease[df$Protease == "ct+tryp"] <- "Ct+Tryp"
  df$Protease[df$Protease == "tryp"] <- "Tryp"
  df$Protease[df$Protease == "aspn"] <- "AspN"
  df
}

annotate_contaminations <- function(df) {
  df$is_not_contamination <- F
  ab_names <- c("h9C12", "Brimab", "PGDM1400", "PGT121", "Umab")
  
  for (ab in ab_names) {
    if (ab == "h9C12") {
      row_idx <- df$match_ig_type %in% c("h9C12_LC", "h9C12-Q97A_HC", "h9C12-WT_HC")
      df[row_idx, "is_not_contamination"] <- df[row_idx, "h9C12-WT"] | df[row_idx, "h9C12-Q97A"]
    } else {
      row_idx <- df$match_ig_type %in% paste0(ab, c("_HC", "_LC"))
      df[row_idx, "is_not_contamination"] <- df[row_idx, ab]  
    }
  }
  print(table(df$is_not_contamination))
  df$is_contamination <- !(df$is_not_contamination)  
  print(table(df$is_contamination))
  return(df)
}


get_imgt_genes <- function() {
  imgt <- list.files(file.path(metadata_path, "search_dbs"), full.names = T)
  imgt_genes <- vector()
  for (i in 1:length(imgt)) {
    tmp <- unlist(read.fasta(imgt[i], as.string = T))
    imgt_genes <- c(imgt_genes, tmp)
  }
  return(toupper(imgt_genes))
}


align_to_ref <- function(reference = reference, peptides = peptides) {
  res <- data.frame(sequence = peptides)
  for (i in 1:length(reference)) {
    res[, i+1] <- str_detect(string = reference[i], pattern = peptides)
    colnames(res)[i+1] <- names(reference)[i]
  }
  if (ncol(res) == 2) {
    res$is_aligned <- res[, 2]
  } else {
    res$is_aligned <- rowSums(res[, -1]) > 0
  }
  res <-  melt(res, id.vars = c("sequence"), variable.name = "ref_name", value.name = "is_aligned")
  res <- res[res$ref_name != "is_aligned" & res$is_aligned, ]
  if (nrow(res) > 0) {
    res <- res[order(res$sequence), ]
  }
  return(res)
}
