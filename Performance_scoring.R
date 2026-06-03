## scoring of profiling depth (peptide and proteins)

total_summary$CV <- total_summary$sd_peptides / total_summary$mean_peptides
scoring <- total_summary

scoring$value <- scoring$mean_peptides / (1+scoring$CV)

library(dplyr)
scoring <- scoring %>%
  group_by(Gradient) %>%
  mutate(
    score = 1 + 4 * (value - min(value, na.rm = TRUE)) /
      (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
  ) %>%
  ungroup()


scoring_protein <- total_summary_protein
scoring_protein$CV <- scoring_protein$sd_protein / scoring_protein$mean_protein
scoring_protein$value <- scoring_protein$mean_protein / (1+scoring_protein$CV)

scoring_protein <- scoring_protein %>%
  group_by(Gradient, Origin) %>%
  mutate(
    score = 1 + 4 * (value - min(value, na.rm = TRUE)) /
      (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
  ) %>%
  ungroup()


spike_scoring <- spike_df %>%
  group_by(Method, Gradient, Group, Dilution_Replicate) %>%
  summarise(
    mean_peptides = mean(Peptide_Count, na.rm = TRUE),
    sd_peptides   = sd(Peptide_Count, na.rm = TRUE),
    CV            = sd_peptides / mean_peptides,
    value         = mean_peptides / (1 + CV),
    .groups = "drop"
  ) %>%
  group_by(Gradient, Group, Dilution_Replicate) %>%
  mutate(
    score = 1 + 4 * (value - min(value, na.rm = TRUE)) /
      (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)),
    score = round(score, 2)
  ) %>%
  ungroup()


spike_scoring

# For spike-in peptides 
# Replace missing scores (e.g., no identifications) with 1 before averaging
# This treats undetected conditions as lowest possible performance

spike_mean_scores <- spike_scoring %>%
  mutate(score = ifelse(is.na(score), 1, score)) %>%
  group_by(Method, Gradient, Group) %>%
  summarise(
    mean_score_raw = mean(score, na.rm = TRUE),
    sd_score       = sd(score, na.rm = TRUE),
    n_conditions   = n(),
    .groups = "drop"
  ) %>%
  group_by(Gradient, Group) %>%
  mutate(
    mean_score = 1 + 4 * (mean_score_raw - min(mean_score_raw, na.rm = TRUE)) /
      (max(mean_score_raw, na.rm = TRUE) - min(mean_score_raw, na.rm = TRUE)),
    mean_score = round(mean_score, 2)
  ) %>%
  ungroup() %>%
  arrange(Gradient, Group, desc(mean_score))


# scoring for quantitative accuracy for spiek-in species
## with coverage penalty 
library(dplyr)

# First calculate total unique peptides per gradient across all methods
peptide_universe <- data1[!grepl("background", data1$species),] %>%
  group_by(Gradient) %>%
  summarise(
    total_unique_peptides = n_distinct(ModifiedPeptide),
    .groups = "drop"
  )

accuracy_scoring <- data1[!grepl("background", data1$species),] %>%
  mutate(
    abs_error = abs(log2ratio - ExpectedRatio)
  ) %>%
  group_by(Method, Gradient) %>%
  summarise(
    mean_error = mean(abs_error, na.rm = TRUE),      # average deviation from expected
    sd_error   = sd(abs_error, na.rm = TRUE),        # variability of deviation
    CV_error   = sd_error / mean_error,              # stability of accuracy
    value      = mean_error * (1 + CV_error),        # lower is better
    n_peptides = n(),
    .groups = "drop"
  ) %>%
  left_join(peptide_universe, by = "Gradient") %>%
  group_by(Gradient) %>%
  mutate(
    raw_score = 1 + 4 * (max(value, na.rm = TRUE) - value) /
      (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)),
    
    coverage_factor = n_peptides / total_unique_peptides,
    
    coverage_adjusted_score = raw_score * coverage_factor,
    
    final_score = 1 + 4 * (coverage_adjusted_score - min(coverage_adjusted_score, na.rm = TRUE)) /
      (max(coverage_adjusted_score, na.rm = TRUE) - min(coverage_adjusted_score, na.rm = TRUE)),
    
    final_score = round(final_score, 2)
  ) %>%
  ungroup()


