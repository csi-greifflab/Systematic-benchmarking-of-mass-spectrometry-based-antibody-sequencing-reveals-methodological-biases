get_coverage_percent <- function(peptides_vector, sequence_name, annotation, mode) {
  V_ref <- unlist(annotation[annotation$sequence_name == sequence_name, "amino_acid_sequence_Vregion"])
  VC_ref <- unlist(annotation[annotation$sequence_name == sequence_name, "amino_acid_sequence_full"])
  CDR3_ref <- unlist(annotation[annotation$sequence_name == sequence_name, "amino_acid_sequence_cdr3"])
  
  alignment_pos <- as.data.frame(do.call(rbind,(str_locate_all(string = VC_ref, 
                                                               pattern = peptides_vector))))
  if (nrow(alignment_pos) == 0)
    if(mode == "both") {
      return(c(vdj = 0, cdr3 = 0))
      }
    else {return(0)
      }
  
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

calculate_miscleavages_trypsin <- function(peptide_string) {
  # Initial condition: cleavage after "K" or "R" unless followed by "P"
  basic_cleavage_sites <- gregexpr("[KR](?![P])", peptide_string, perl=TRUE)[[1]]
  
  # Exception patterns where trypsin cannot cleave, including the end of the peptide where it's already cut
  block_patterns <- c("CKD", "DKD", "CKH", "CKY", "CRK", "RRR", "RRH", "K$", "R$")
  
  # Find all the exceptions and remove them from the basic cleavage sites
  block_sites <- unlist(lapply(block_patterns, function(pattern) {
    gregexpr(pattern, peptide_string, perl=TRUE)[[1]]}))
  
  # Merge and sort cleavage sites while excluding exceptions (-1 positions indicate no match)
  cleavage_sites <- setdiff(
    basic_cleavage_sites[basic_cleavage_sites > 0], 
    block_sites[block_sites > 0])
  
  # The number of miscleavages is the number of cleavage sites.
  miscleavages <- max(0, length(cleavage_sites))
  
  return(miscleavages)
}

calculate_miscleavages_chymotrypsin <- function(peptide_string) {
  # Initial condition: cleavage after "F" or "Y" or "W" unless followed by "P"
  basic_cleavage_sites <- gregexpr("[FYW](?![P])", peptide_string, perl=TRUE)[[1]]
  
  # Exception patterns where Chymotrypsin cannot cleave, including the end of the peptide where it's already cut
  block_patterns <- c("WM", "F$","W$", "Y$")
  
  # Find all the exceptions and remove them from the basic cleavage sites
  block_sites <- unlist(lapply(block_patterns, function(pattern) {
    gregexpr(pattern, peptide_string, perl=TRUE)[[1]]}))
  
  # Merge and sort cleavage sites while excluding exceptions (-1 positions indicate no match)
  cleavage_sites <- setdiff(
    basic_cleavage_sites[basic_cleavage_sites > 0], 
    block_sites[block_sites > 0])
  
  # The number of miscleavages is the number of cleavage sites.
  miscleavages <- max(0, length(cleavage_sites))
  
  return(miscleavages)
}

calculate_miscleavages_aspn <- function(peptide_string) {
  # Initial condition: cleavage before "D"
  basic_cleavage_sites <- gregexpr("(?<=D)", peptide_string, perl=TRUE)[[1]] - 1
  
  # Exception patterns where AspN cannot cleave, including the beginning of the peptide where it's already cut
  block_patterns <- c("^D")
  
  # Find all the exceptions and remove them from the basic cleavage sites
  block_sites <- unlist(lapply(block_patterns, function(pattern) {
    gregexpr(pattern, peptide_string, perl=TRUE)[[1]]}))
  
  # Merge and sort cleavage sites while excluding exceptions (-1 positions indicate no match)
  cleavage_sites <- setdiff(
    basic_cleavage_sites[basic_cleavage_sites > 0], 
    block_sites[block_sites > 0])
  
  # The number of miscleavages is the number of cleavage sites.
  miscleavages <- max(0, length(cleavage_sites))
  
  return(miscleavages)
}  

