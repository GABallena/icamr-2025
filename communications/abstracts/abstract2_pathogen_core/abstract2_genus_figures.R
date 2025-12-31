#!/usr/bin/env Rscript
# Portfolio-safe script (paths/IDs generalized; inputs not included)

# ------- Working directory handling -------
# When run via `Rscript`, set working directory to the script's directory.
# This keeps the script portable (no machine-specific paths).
try({
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if(length(file_arg) == 1){
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    setwd(dirname(script_path))
  }
}, silent = TRUE)
cat("Working directory:", getwd(), "\n")

# Genus figures: core heatmap, prevalence bar, burden–richness scatter
# Expects a table with either:
#  (A) long format:  sample, genus, rel_abundance
#  (B) wide format:  first column = sample, other columns = genera (relative abundance)
#
# Outputs:
#  results/figs/genus_core_heatmap.png
#  results/figs/genus_core_prevalence_bar.png
#  results/figs/pathogen_burden_vs_richness.png
#  results/figs/genus_core_metrics.tsv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(scales)
})

# ------------------ CONFIG (edit as needed) ------------------
setwd("<PROJECT_ROOT>")

infile <- "results/tables/pathogen_genus_abundance.tsv"  # change if needed
# Optional curated pathogen list (single column "genus"); set to NULL to skip
pathogen_list_file <- NULL  # e.g., "metadata/pathogen_genus_list.tsv"

min_prev  <- 0.70      # core prevalence threshold (≥70% sites)
min_abund <- 0.001     # 0.1% relative abundance threshold
fig_dir   <- "results/figs"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
# -------------------------------------------------------------

message("Reading: ", infile)
stopifnot(file.exists(infile))
dat0 <- suppressMessages(read_tsv(infile, show_col_types = FALSE))

# Detect long vs wide; normalise to long: sample, genus, rel_abundance
if (all(c("sample","genus","rel_abundance") %in% names(dat0))) {
  dat <- dat0 %>%
    select(sample, genus, rel_abundance)
} else {
  # assume first column is sample, others are genera
  stopifnot(ncol(dat0) >= 2)
  sample_col <- names(dat0)[1]
  dat <- dat0 %>%
    rename(sample = !!sample_col) %>%
    pivot_longer(-sample, names_to = "genus", values_to = "rel_abundance")
}

# Clean up types and NAs
dat <- dat %>%
  mutate(
    sample = as.character(sample),
    genus  = as.character(genus),
    rel_abundance = suppressWarnings(as.numeric(rel_abundance))
  ) %>%
  filter(!is.na(sample), !is.na(genus), !is.na(rel_abundance))

# Optional: restrict to curated pathogen genera
if (!is.null(pathogen_list_file)) {
  stopifnot(file.exists(pathogen_list_file))
  patholist <- read_tsv(pathogen_list_file, show_col_types = FALSE) %>%
    transmute(genus = as.character(genus))
  dat <- dat %>% semi_join(patholist, by = "genus")
}

# Basic per-sample summaries
site_list <- dat %>% distinct(sample) %>% arrange(sample) %>% pull(sample)
n_sites   <- length(site_list)

burden_richness <- dat %>%
  group_by(sample) %>%
  summarise(
    pathogen_burden = sum(rel_abundance, na.rm = TRUE),
    pathogen_richness = sum(rel_abundance > 0, na.rm = TRUE),
    .groups = "drop"
  )

# Core definition metrics per genus
core_metrics <- dat %>%
  group_by(genus) %>%
  summarise(
    prevalence = mean(rel_abundance >= min_abund, na.rm = TRUE), # fraction of sites
    median_abund = median(rel_abundance, na.rm = TRUE),
    mean_abund   = mean(rel_abundance, na.rm = TRUE),
    n_sites_nonzero = sum(rel_abundance >= min_abund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(prevalence), desc(median_abund))

core_genera <- core_metrics %>%
  filter(prevalence >= min_prev) %>%
  arrange(desc(prevalence), desc(median_abund)) %>%
  pull(genus)

message("Detected ", length(core_genera), " core genera at ≥",
        100*min_prev, "% prevalence and ≥", 100*min_abund, "% abundance.")

# Order sites by pathogen burden to get a nice heatmap ordering
site_order <- burden_richness %>%
  arrange(desc(pathogen_burden)) %>%
  pull(sample)

# Heatmap data: only core genera
hm <- dat %>%
  filter(genus %in% core_genera) %>%
  mutate(
    sample = factor(sample, levels = site_order),
    genus  = factor(genus,  levels = core_genera)
  )

# If you have a lot of dynamic range, you can optionally transform:
# hm <- hm %>% mutate(rel_abundance = asinh(rel_abundance * 1000))  # comment out by default

p_heat <- ggplot(hm, aes(x = sample, y = genus, fill = rel_abundance)) +
  geom_tile() +
  scale_fill_viridis_c(labels = percent_format(accuracy = 0.1),
                       name = "Rel. abundance") +
  labs(x = "Site", y = "Core pathogen genera",
       title = "Core pathogen genera across sites",
       subtitle = paste0("Core = prevalence ≥", 100*min_prev,
                         "% at abundance ≥", 100*min_abund, "% (n=",
                         length(core_genera), " genera; ", n_sites, " sites)")) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(file.path(fig_dir, "genus_core_heatmap.png"),
       p_heat, width = 10, height = 6, dpi = 600)

# Prevalence bar with median abundance dot
core_for_bar <- core_metrics %>%
  filter(genus %in% core_genera) %>%
  mutate(genus = factor(genus, levels = core_genera))

p_bar <- ggplot(core_for_bar, aes(x = genus, y = prevalence)) +
  geom_col() +
  geom_point(aes(y = pmin(1, median_abund / max(median_abund + 1e-12))),  # scaled dot (optional)
             inherit.aes = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Core pathogen genera",
       y = "Prevalence across sites",
       title = "Core genus prevalence (dot indicates scaled median abundance)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave(file.path(fig_dir, "genus_core_prevalence_bar.png"),
       p_bar, width = 10, height = 5, dpi = 600)

# Burden vs richness scatter
p_scatter <- ggplot(burden_richness, aes(x = pathogen_richness, y = pathogen_burden)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(x = "Pathogen richness (genus count)",
       y = "Pathogen burden (sum of pathogen genera, % of reads)",
       title = "Pathogen burden vs richness across sites") +
  theme_minimal(base_size = 11)

ggsave(file.path(fig_dir, "pathogen_burden_vs_richness.png"),
       p_scatter, width = 6.5, height = 5, dpi = 600)

# Save metrics table
out_metrics <- file.path(fig_dir, "genus_core_metrics.tsv")
write_tsv(core_metrics, out_metrics)

# Console summary (handy for poster text)
summary_lines <- c(
  "=== GENUS FIGURES SUMMARY ===",
  paste0("Sites: ", n_sites),
  paste0("Core prevalence threshold: ", 100*min_prev, "%; abundance threshold: ", 100*min_abund, "%"),
  paste0("Core genera (n=", length(core_genera), "): ", paste(core_genera, collapse = ", ")),
  paste0("Median burden: ", percent(median(burden_richness$pathogen_burden), accuracy = 0.1),
         " (IQR ", percent(quantile(burden_richness$pathogen_burden, 0.25), accuracy = 0.1),
         "–", percent(quantile(burden_richness$pathogen_burden, 0.75), accuracy = 0.1), ")"),
  paste0("Median richness: ", signif(median(burden_richness$pathogen_richness), 3),
         " (IQR ", signif(quantile(burden_richness$pathogen_richness, 0.25), 3),
         "–", signif(quantile(burden_richness$pathogen_richness, 0.75), 3), ")"),
  paste0("Metrics saved: ", out_metrics),
  paste0("Figures saved:\n  - ", file.path(fig_dir, "genus_core_heatmap.png"),
         "\n  - ", file.path(fig_dir, "genus_core_prevalence_bar.png"),
         "\n  - ", file.path(fig_dir, "pathogen_burden_vs_richness.png"))
)
cat(paste0(summary_lines, collapse = "\n"), "\n")
