#!/usr/bin/env Rscript
# Regenerate MAG figures at 600 DPI for HD poster
# This script sources the relevant sections from abstract6.R
cat("\n\n")
cat("   Regenerating MAG Figures at 600 DPI                   \n")
cat("\n\n")
# -------------------------------------------------------------------------
# Portfolio-ready note:
# - No hard-coded local paths.
# - Run from the project root (or set PROJECT_ROOT env var).
# -------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = getwd())
try(setwd(PROJECT_ROOT), silent = TRUE)


suppressPackageStartupMessages({
  library(tidyverse)
  library(ggraph)
  library(igraph)
})

# Color palette from abstract6.R
citystat_core12 <- c(
  Navy      = "#233B5D", Blue      = "#2F6DA8", DeepTeal  = "#1C7C7D", Cyan      = "#17BEBB",
  Forest    = "#2E7D32", Olive     = "#7A8F30", Mustard   = "#C7A41A", Orange    = "#E07A2D",
  Vermilion = "#D0432B", Crimson   = "#B01E2F", Plum      = "#7A3E9D", Indigo    = "#3B4BA3"
)

HD_DPI <- 600
OUT_DIR <- "HD_poster_figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Figure 1: MAG Taxonomy Sunburst (600 DPI)
# ============================================================================
cat("[1/2] Regenerating MAG Taxonomy Sunburst at 600 DPI...\n")

nodes_path <- file.path("mag_analysis", "mag_taxonomy_sunburst_nodes.tsv")
if (!file.exists(nodes_path)) {
  cat("   Warning: mag_taxonomy_sunburst_nodes.tsv not found.\n")
  cat("  -> Please run abstract6.R first to generate MAG analysis data.\n")
  cat("  -> Skipping sunburst regeneration.\n\n")
} else {
  nodes_df <- read_tsv(nodes_path, show_col_types = FALSE)
  
  # Build edges from id and parent columns
  edges_list <- nodes_df %>%
    filter(!is.na(parent), parent != "") %>%
    select(from = parent, to = id)
  
  # Add rank_label if not present
  if (!"rank_label" %in% names(nodes_df)) {
    nodes_df <- nodes_df %>%
      mutate(
        rank_label = case_when(
          depth == 1 ~ "Domain",
          depth == 2 ~ "Phylum",
          depth == 3 ~ "Class",
          depth == 4 ~ "Order",
          depth == 5 ~ "Family",
          depth >= 6 ~ "Genus",
          TRUE ~ "Unknown"
        )
      )
  }
  
  g_sun <- graph_from_data_frame(edges_list, directed = TRUE, vertices = nodes_df)
  
  p_sun <- ggraph(g_sun, layout = 'partition', circular = TRUE) +
    geom_node_arc_bar(aes(fill = rank_label), linewidth = 0.2, color = "white") +
    scale_fill_manual(
      values = c(
        "Domain" = "#2E3A59", "Phylum" = "#2B7A78", "Class" = "#17BEBB",
        "Order" = "#A8DADC", "Family" = "#F4A261", "Genus" = "#E76F51"
      ),
      name = "Taxonomic Rank",
      na.value = "grey80"
    ) +
    coord_fixed() +
    theme_void(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = "MAG Taxonomic Composition (Sunburst)")
  
  out_png <- file.path(OUT_DIR, "mag_taxonomy_sunburst_HD.png")
  ggsave(out_png, p_sun, width = 10, height = 10, dpi = HD_DPI, bg = "white")
  cat(sprintf("  ✓ Saved: %s [%.1f MB]\n", basename(out_png), file.size(out_png) / 1024^2))
}

# ============================================================================
# Figure 2: MAG Drug Class Stack (600 DPI)
# ============================================================================
cat("\n[2/2] Regenerating MAG Drug Class Stack at 600 DPI...\n")

arg_strict_path <- "mag_analysis/arg_to_mag_strict.tsv"
if (!file.exists(arg_strict_path)) {
  cat("   Warning: arg_to_mag_strict.tsv not found.\n")
  cat("  -> Please run abstract6.R first to generate MAG ARG analysis.\n")
  cat("  -> Skipping drug class stack regeneration.\n\n")
} else {
  # Read ARG data
  arg_raw <- read_tsv(arg_strict_path, show_col_types = FALSE)
  
  # Clean bin IDs (from abstract6.R logic)
  clean_bin_id <- function(x) {
    x <- gsub("\\.orig\\.fa$", "", x)
    x <- gsub("\\.fa$", "", x)
    x <- gsub("\\.fasta$", "", x)
    x <- gsub("\\.fna$", "", x)
    x
  }
  
  # Rename bin_id to mag_id and clean
  arg_raw <- arg_raw %>%
    rename(mag_id = bin_id) %>%
    mutate(mag_id = clean_bin_id(mag_id)) %>%
    filter(!is.na(mag_id), mag_id != "", mag_id != "unbinned")
  
  # Get top MAGs by ARG richness
  top_n_mags <- 30
  arg_top <- arg_raw %>%
    group_by(mag_id) %>%
    summarise(
      n_distinct_ARO = n_distinct(ARO_accession),
      drug_classes = paste(unique(na.omit(drug_classes)), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_distinct_ARO)) %>%
    slice_head(n = top_n_mags)
  
  # Create drug class long format
  drug_long <- arg_top %>%
    separate_rows(drug_classes, sep = ";") %>%
    filter(!is.na(drug_classes), drug_classes != "NA", drug_classes != "") %>%
    mutate(drug_classes_plot = drug_classes)
  
  # Assign colors
  palette_classes <- setNames(
    rep_len(citystat_core12, length.out = length(unique(drug_long$drug_classes_plot))),
    sort(unique(drug_long$drug_classes_plot))
  )
  
  # Create plot
  p_drug <- drug_long %>%
    count(mag_id, drug_classes_plot) %>%
    left_join(arg_top %>% select(mag_id, n_distinct_ARO), by = "mag_id") %>%
    arrange(desc(n_distinct_ARO)) %>%
    mutate(mag_label = factor(mag_id, levels = unique(mag_id))) %>%
    ggplot(aes(mag_label, n, fill = drug_classes_plot)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = palette_classes, name = "Drug Class") +
    labs(
      title = "ARG Drug Class Composition per MAG (Top 30)",
      y = "Hit count (drug class)",
      x = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(face = "italic", size = 9),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0)
    )
  
  out_png <- file.path(OUT_DIR, "mag_drug_class_stack_HD.png")
  ggsave(out_png, p_drug, width = 12, height = 6, dpi = HD_DPI, bg = "white")
  cat(sprintf("  ✓ Saved: %s [%.1f MB]\n", basename(out_png), file.size(out_png) / 1024^2))
}

cat("\n\n")
cat("   MAG HD Figure Regeneration Complete                   \n")
cat("\n")
cat("  Location: HD_poster_figures/                              \n")
cat("  Resolution: 600 DPI                                       \n")
cat("\n\n")

# List all HD figures
cat("All HD Poster Figures (600 DPI):\n")
hd_files <- list.files(OUT_DIR, pattern = "_HD\\.png$", full.names = TRUE)
for (f in hd_files) {
  fsize <- file.size(f) / 1024^2
  dims <- system(sprintf("file %s | grep -oP '\\d+ x \\d+'", shQuote(f)), intern = TRUE)
  cat(sprintf("  ✓ %-50s [%.1f MB, %s]\n", basename(f), fsize, dims))
}
cat("\n")
