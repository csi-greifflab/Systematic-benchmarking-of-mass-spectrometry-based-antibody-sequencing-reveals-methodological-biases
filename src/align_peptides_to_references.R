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
  res <-  reshape2::melt(res, id.vars = c("sequence"), variable.name = "ref_name", value.name = "is_aligned")
  res <- res[res$ref_name != "is_aligned" & res$is_aligned, ]
  if (nrow(res) > 0) {
    res <- res[order(res$sequence), ]
  }
  return(res)
}
