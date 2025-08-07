get_intensity_ratio <- function(df) {
  
  df_filtered <- df %>%
    filter(
      tool != "Casanovo",
      !match_ig_type %in% c("h9C12-WT_HC", "h9C12-Q97A_HC", "h9C12_LC"), 
      Protease != "AspN"
    ) %>%
    mutate(
      mab = str_remove(match_ig_type, "_[HL]C$")
    )
  
  res_ratio <- df_filtered %>%
    left_join(
      df_filtered,
      by = c("Sequence", "Protease", "tool", "match_ig_type", "mab", "is_cdr3_related", "run", "has_blood"),
      relationship = "many-to-many"
    ) %>%
    filter(Rawfilenumber.x != Rawfilenumber.y) %>%
    transmute(
      Sequence,
      Protease,
      tool,
      match_ig_type,
      mab,
      is_cdr3_related, 
      run, 
      has_blood, 
      intensity1 = intensity.x,
      intensity2 = intensity.y,
      concentration1 = concentration.x,
      concentration2 = concentration.y,
      conc_ratio = concentration.x / concentration.y,
      intensity_ratio = intensity.x / intensity.y
    ) 
  
  idx <- res_ratio$conc_ratio < 1
  res_ratio$conc_ratio[idx] <- 1/res_ratio$conc_ratio[idx]
  res_ratio$intensity_ratio[idx] <- 1/res_ratio$intensity_ratio[idx]
  
  res_ratio <- res_ratio %>%
    mutate(
      log_conc_ratio = log10(conc_ratio), 
      log_intensity_ratio = log10(intensity_ratio)
    )
  
  return(res_ratio)
}