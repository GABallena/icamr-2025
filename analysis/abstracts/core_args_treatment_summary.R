#!/usr/bin/env Rscript
# Create summary table matching core ARGs with treatment effects
# -------------------------------------------------------------------------
# Portfolio-ready note:
# - No hard-coded local paths.
# - Run from the project root (or set PROJECT_ROOT env var).
# -------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = getwd())
try(setwd(PROJECT_ROOT), silent = TRUE)


suppressPackageStartupMessages({
  library(tidyverse)
})

# Load core ARGs
core <- read_tsv("results/core_classes.tsv", show_col_types = FALSE)

# Load ΔCLR data
dclr <- read_tsv("abstract3/stage2/delta_clr_per_pair.tsv", show_col_types = FALSE) %>%
  pivot_longer(-site_core, names_to = "arg_class_short", values_to = "delta_clr") %>%
  filter(!is.na(delta_clr), is.finite(delta_clr))

# Calculate mean ΔCLR
dclr_summary <- dclr %>%
  group_by(arg_class_short) %>%
  summarise(
    mean_dclr = mean(delta_clr),
    median_dclr = median(delta_clr),
    sd_dclr = sd(delta_clr),
    n_sites = n(),
    .groups = "drop"
  ) %>%
  mutate(se_dclr = sd_dclr / sqrt(n_sites))

# Create mapping between naming conventions
# abstract3 uses shortened names, abstract2 uses full CARD names
name_mapping <- tribble(
  ~short, ~pattern,
  "aminoglycoside", "aminoglycoside",
  "beta_lactam", "lactam",
  "chloramphenicol", "chloramphenicol|phenicol",
  "florfenicol", "florfenicol|phenicol",
  "macrolide", "macrolide",
  "multidrug", "disinfecting|multidrug",
  "polymyxin", "polymyxin|peptide",
  "quinolone", "quinolone|fluoroquinolone",
  "sulfonamide", "sulfonamide|sulfone",
  "tetracycline", "tetracycline",
  "trimethoprim", "trimethoprim|diaminopyrimidine",
  "vancomycin", "vancomycin|glycopeptide"
)

# Match core ARGs with ΔCLR data
matched_results <- map_dfr(1:nrow(core), function(i) {
  core_class <- core$arg_class[i]
  
  # Try to find matching ΔCLR entries
  matches <- dclr_summary %>%
    filter(str_detect(tolower(core_class), 
                     pattern = paste(arg_class_short, collapse = "|")))
  
  if (nrow(matches) == 0) {
    # Try fuzzy matching via mapping table
    for (j in 1:nrow(name_mapping)) {
      if (str_detect(tolower(core_class), name_mapping$pattern[j])) {
        matches <- dclr_summary %>%
          filter(arg_class_short == name_mapping$short[j])
        break
      }
    }
  }
  
  if (nrow(matches) > 0) {
    tibble(
      arg_class = core_class,
      prevalence = core$prevalence_over_X[i],
      mean_cpm = core$mean_CPM[i],
      mean_dclr = matches$mean_dclr[1],
      se_dclr = matches$se_dclr[1],
      effect_direction = case_when(
        matches$mean_dclr[1] > 0.5 ~ "↑ Treated",
        matches$mean_dclr[1] < -0.5 ~ "↑ Untreated",
        TRUE ~ "~"
      )
    )
  } else {
    tibble(
      arg_class = core_class,
      prevalence = core$prevalence_over_X[i],
      mean_cpm = core$mean_CPM[i],
      mean_dclr = NA_real_,
      se_dclr = NA_real_,
      effect_direction = "No data"
    )
  }
})

# Create formatted summary table
summary_table <- matched_results %>%
  arrange(desc(abs(mean_dclr))) %>%
  mutate(
    prevalence_pct = sprintf("%.0f%%", prevalence * 100),
    cpm_formatted = sprintf("%.1f", mean_cpm),
    dclr_formatted = ifelse(is.na(mean_dclr), 
                           "-", 
                           sprintf("%.2f ± %.2f", mean_dclr, se_dclr))
  ) %>%
  select(
    `ARG Class` = arg_class,
    `Prevalence` = prevalence_pct,
    `Mean CPM` = cpm_formatted,
    `ΔCLR (Treated - Untreated)` = dclr_formatted,
    `Effect` = effect_direction
  )

# Save as TSV
write_tsv(matched_results, "combined_reads/core_args_with_treatment_effects.tsv")

# Create human-readable version
write_lines(
  c(
    "# Core ARGs with Treatment Effects",
    "",
    sprintf("Total core ARG classes: %d", nrow(core)),
    sprintf("Matched with treatment data: %d", sum(!is.na(matched_results$mean_dclr))),
    "",
    "Top 15 by absolute treatment effect:",
    ""
  ),
  "combined_reads/core_args_summary.txt"
)

summary_table %>%
  head(15) %>%
  format_tsv() %>%
  write_lines("combined_reads/core_args_summary.txt", append = TRUE)

# Also display to console
message("\n=== Core ARGs with Treatment Effects ===\n")
print(summary_table %>% head(20), n = 20)

message(sprintf("\n✓ Summary table saved to: combined_reads/core_args_with_treatment_effects.tsv"))
message(sprintf("✓ Formatted summary saved to: combined_reads/core_args_summary.txt"))
