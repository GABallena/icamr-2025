# -------------------------------------------------------------------------
# Portfolio-ready note:
# - No hard-coded local paths.
# - Run from the project root (or set PROJECT_ROOT env var).
# -------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = getwd())
try(setwd(PROJECT_ROOT), silent = TRUE)

#!/usr/bin/env Rscript
# Create a comprehensive combined figure for pathogenic MAGs
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(scales)

# Color palette
citystat_core12 <- c(
  Navy = "#233B5D", Blue = "#2F6DA8", DeepTeal = "#1C7C7D",
  Cyan = "#17BEBB", Forest = "#2E7D32", Olive = "#7A8F30",
  Mustard = "#C7A41A", Orange = "#E07A2D", Vermilion = "#D0432B",
  Crimson = "#B01E2F", Plum = "#7A3E9D", Indigo = "#3B4BA3"
)

classification_colors <- c(
  "High-quality MAG" = citystat_core12[["Blue"]],
  "Medium-quality MAG" = citystat_core12[["Orange"]],
  "Low-quality MAG" = citystat_core12[["Crimson"]]
)

# Load data
pathogen_summary <- read_tsv("mag_analysis/pathogen_summary.tsv", show_col_types = FALSE)
unified <- read_tsv("mag_analysis/MAG_quality_summary.tsv", show_col_types = FALSE)

# Load pathogen list
pat_lines <- readLines("clean_pathogens.txt", warn = FALSE)
pat_genera <- tibble(line = pat_lines) %>%
  filter(nchar(trimws(line)) > 0) %>%
  mutate(genus = stringr::str_extract(line, "^[A-Za-z]+")) %>%
  filter(!is.na(genus)) %>%
  distinct(genus)

# Identify pathogenic MAGs
pathogen_mags <- unified %>%
  filter(!is.na(genus), genus %in% pat_genera$genus) %>%
  mutate(
    completeness = coalesce(Completeness_chk2, Completeness_chk1),
    contamination = coalesce(Contamination_chk2, Contamination_chk1)
  )

# Panel A: Top pathogenic genera bar chart
p1 <- pathogen_summary %>%
  arrange(desc(n_mags)) %>%
  slice_head(n = 10) %>%
  mutate(pathogen_genus = factor(pathogen_genus, levels = rev(pathogen_genus))) %>%
  ggplot(aes(x = pathogen_genus, y = n_mags)) +
  geom_col(fill = citystat_core12[["Vermilion"]], alpha = 0.85) +
  geom_text(aes(label = n_mags), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(
    title = "A) Top 10 Pathogenic Genera",
    x = NULL,
    y = "Number of MAGs"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(face = "bold", size = 11)
  )

# Panel B: Quality distribution stacked bar
pathogen_quality_long <- pathogen_summary %>%
  arrange(desc(n_mags)) %>%
  slice_head(n = 10) %>%
  tidyr::pivot_longer(
    cols = c(hq_count, mq_count, lq_count),
    names_to = "quality_class",
    values_to = "count"
  ) %>%
  mutate(
    quality_label = case_when(
      quality_class == "hq_count" ~ "High-quality",
      quality_class == "mq_count" ~ "Medium-quality",
      quality_class == "lq_count" ~ "Low-quality"
    ),
    quality_label = factor(quality_label, levels = c("High-quality", "Medium-quality", "Low-quality"))
  )

p2 <- pathogen_quality_long %>%
  mutate(pathogen_genus = factor(pathogen_genus, levels = rev(unique(pathogen_quality_long$pathogen_genus)))) %>%
  ggplot(aes(x = pathogen_genus, y = count, fill = quality_label)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = classification_colors, name = "Quality") +
  coord_flip() +
  labs(
    title = "B) Quality Distribution",
    x = NULL,
    y = "Number of MAGs"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 11)
  )

# Panel C: Completeness vs Contamination scatter
p3 <- pathogen_mags %>%
  ggplot(aes(x = completeness, y = contamination, color = classification_mimag)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = classification_colors, name = "Quality") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey50", alpha = 0.6) +
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey50", alpha = 0.6) +
  labs(
    title = "C) Quality Assessment",
    x = "Completeness (%)",
    y = "Contamination (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 11)
  )

# Panel D: Sample distribution
pathogen_per_sample <- pathogen_mags %>%
  filter(!is.na(sample)) %>%
  count(sample) %>%
  arrange(desc(n)) %>%
  slice_head(n = 15)

p4 <- pathogen_per_sample %>%
  mutate(sample = factor(sample, levels = rev(sample))) %>%
  ggplot(aes(x = sample, y = n)) +
  geom_col(fill = citystat_core12[["Crimson"]], alpha = 0.85) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.5) +
  coord_flip() +
  labs(
    title = "D) Top 15 Samples by Pathogen Count",
    x = NULL,
    y = "Pathogenic MAGs"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11)
  )

# Combine all panels
combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Pathogenic MAG Analysis",
    subtitle = sprintf("%d pathogenic MAGs across %d genera from %d total MAGs",
                      nrow(pathogen_mags), 
                      n_distinct(pathogen_mags$genus),
                      nrow(unified)),
    caption = "High-quality: ≥90% complete, ≤5% contamination | Medium: ≥50% complete, <10% contamination",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  )

# Save
ggsave("mag_analysis/pathogen_mags_combined_figure.png", 
       combined, width = 14, height = 10, dpi = 300, bg = "white")

cat("✓ Created combined pathogenic MAG figure: mag_analysis/pathogen_mags_combined_figure.png\n")
cat("  - Panel A: Top 10 pathogenic genera\n")
cat("  - Panel B: Quality distribution\n")
cat("  - Panel C: Quality assessment scatter\n")
cat("  - Panel D: Top 15 samples by pathogen count\n")
