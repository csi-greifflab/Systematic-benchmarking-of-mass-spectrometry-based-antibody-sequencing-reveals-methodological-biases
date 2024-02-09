get_coverage_percent <- function(peptides_vector, sequence_name, annotation, mode) {
  V_ref <- unlist(annotation[annotation$sequence_name == sequence_name, "amino_acid_sequence_Vregion"])
  VC_ref <- unlist(annotation[annotation$sequence_name == sequence_name, "amino_acid_sequence_full"])
  CDR3_ref <- unlist(annotation[annotation$sequence_name == sequence_name, "amino_acid_sequence_cdr3"])
  
  alignment_pos <- as.data.frame(do.call(rbind,(str_locate_all(string = VC_ref, 
                                                             pattern = peptides_vector))))
  
  if (mode == "cdr3") {
    ref <- CDR3_ref
  } else if (mode == "vdj") {
    ref <- V_ref
  } else if (mode == "both") {
    return (c(vdj = get_coverage_percent(peptides_vector, sequence_name, annotation, "vdj"), 
              cdr3 = get_coverage_percent(peptides_vector, sequence_name, annotation, "cdr3")))
  } else {
    stop("Wrong mode parameter value. Should be cdr3 or vdj.")
  }
  
  if (nrow(alignment_pos) == 0)
    return(0)
    
  counts <- replicate(nchar(VC_ref), 0)
  for (i in 1:nrow(alignment_pos)) {
    counts[1:length(counts) %in% c(alignment_pos[i, ]$start : alignment_pos[i, ]$end)] <- 
      counts[alignment_pos[i, ]$start : alignment_pos[i, ]$end] + 1
  }
  
  ref_location <- str_locate(string = VC_ref, pattern = ref)
  ref_length <- ref_location[2] - ref_location[1] + 1
  ref_positions <- c(ref_location[1] : ref_location[2])
  coverage <- sum(counts[1:length(counts) %in% ref_positions] > 0) / (ref_length) * 100
  
  return(coverage)
}

get_coverage_per_ab_sample <- function(df, mode) { 
  coverage <- data.frame(matrix(ncol = length(HC_and_LC_names), nrow = 70))
  rownames(coverage) <- 1:70
  colnames(coverage) <- HC_and_LC_names
  
  for (sample_id in 1:70){
    for (seq_name_id in 1:(ncol(coverage))) {
      seq_name <- colnames(coverage)[seq_name_id]
      if (mode == "cdr3") {
        peptides <- unique(df[df$Sample == sample_id & df$match_ig_type == seq_name & df$is_cdr3_related, ]$Sequence)
        coverage[sample_id, seq_name_id] <- round(get_coverage_percent(peptides, seq_name, cdr3, "cdr3"), 0)
      } else if (mode == "vdj") {
        peptides <- unique(df[df$Sample == sample_id & df$match_ig_type == seq_name, ]$Sequence)
        coverage[sample_id, seq_name_id] <- round(get_coverage_percent(peptides, seq_name, cdr3, "vdj"), 0)
      } else {
        stop("Wrong mode parameter value. Should be cdr3 or vdj.")
      }
    }
  }
  coverage
}



