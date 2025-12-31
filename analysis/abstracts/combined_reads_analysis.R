#!/usr/bin/env Rscript
# Combined Read-Level ARG Analysis
# Integrates core ARG detection (abstract2) with treatment effects (abstract3)
#
# Input:
#   - results/core_classes.tsv (from abstract2)
#   - results/class_abundance_cpm.tsv (from abstract2)
#   - abstract3/stage2/* (ΔCLR data from abstract3)
#
# Output:
#   - combined_reads/integrated_arg_analysis.png (4-panel figure)
# -------------------------------------------------------------------------
# Portfolio-ready note:
# - No hard-coded local paths.
# - Run from the project root (or set PROJECT_ROOT env var).
# -------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = getwd())
try(setwd(PROJECT_ROOT), silent = TRUE)


suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

# Create output directory
dir.create("combined_reads", showWarnings = FALSE, recursive = TRUE)

# ========== 1. Load Core ARGs from abstract2 ==========
message("Loading core ARG classes from abstract2...")
core_args <- read_tsv("results/core_classes.tsv", show_col_types = FALSE)

# Load full CPM data
arg_cpm <- read_tsv("results/class_abundance_cpm.tsv", show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  filter(!is.na(cpm), cpm > 0)

# Summarize core ARG stats
core_summary <- arg_cpm %>%
  filter(arg_class %in% core_args$arg_class) %>%
  group_by(arg_class) %>%
  summarise(
    n_samples = n_distinct(sample),
    mean_cpm = mean(cpm, na.rm = TRUE),
    median_cpm = median(cpm, na.rm = TRUE),
    prevalence = n_samples / n_distinct(arg_cpm$sample),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cpm))

message(sprintf("Found %d core ARG classes detected in %d samples",
                nrow(core_summary), n_distinct(arg_cpm$sample)))

# ========== 2. Load Treatment Effects from abstract3 ==========
message("Loading treatment effect data from abstract3...")

# Load ΔCLR data from abstract3/stage2
dclr_file <- "abstract3/stage2/delta_clr_per_pair.tsv"

if (file.exists(dclr_file)) {
  dclr_wide <- read_tsv(dclr_file, show_col_types = FALSE)
  message(sprintf("Loaded ΔCLR data: %d sites, %d ARG classes", 
                  nrow(dclr_wide), ncol(dclr_wide) - 1))
  
  # Convert from wide to long format
  dclr_data <- dclr_wide %>%
    pivot_longer(-site_core, names_to = "arg_class", values_to = "delta_clr") %>%
    filter(!is.na(delta_clr), is.finite(delta_clr))
  
  message(sprintf("Reshaped to long format: %d observations", nrow(dclr_data)))
} else {
  dclr_data <- NULL
  message("ΔCLR data file not found at: ", dclr_file)
}

# ========== 3. Load Sample Metadata ==========
meta_file <- "metadata/sample_metadata.tsv"
if (file.exists(meta_file)) {
  metadata <- read_tsv(meta_file, show_col_types = FALSE) %>%
    janitor::clean_names() %>%
    mutate(
      treatment_status = case_when(
        str_detect(tolower(sample_type), "untreated") ~ "Untreated",
        str_detect(tolower(sample_type), "treated") ~ "Treated",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(treatment_status)) %>%
    select(sample_code, treatment_status, sample_type, sample_description)
  
  # Join with CPM data
  arg_cpm_meta <- arg_cpm %>%
    left_join(metadata, by = c("sample" = "sample_code"))
  
  message(sprintf("Joined metadata: %d samples with treatment info",
                  sum(!is.na(arg_cpm_meta$treatment_status))))
} else {
  arg_cpm_meta <- arg_cpm %>% mutate(treatment_status = NA_character_)
  message("Metadata file not found, proceeding without treatment labels")
}

# ========== 4. Create Visualizations ==========

# Define color palette
p4_pal <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

# Panel A: Core ARG Prevalence vs Abundance
message("Creating Panel A: Core ARG prevalence vs abundance...")

top_core <- core_summary %>% head(15)

panel_a <- ggplot(top_core, aes(x = prevalence, y = mean_cpm)) +
  geom_point(aes(size = median_cpm), color = p4_pal[2], alpha = 0.7) +
  geom_text_repel(
    aes(label = str_trunc(arg_class, 30)),
    size = 3, max.overlaps = 15, segment.size = 0.2
  ) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_log10(labels = scales::comma) +
  scale_size_continuous(name = "Median CPM", range = c(2, 10)) +
  labs(
    title = "Core ARG Classes: Prevalence vs Abundance",
    subtitle = sprintf("Top %d core ARGs (≥0.1 CPM in ≥70%% sites)", nrow(top_core)),
    x = "Prevalence Across Samples",
    y = "Mean Abundance (CPM, log scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

# Panel B: Detection Frequency by Treatment Status
message("Creating Panel B: Detection by treatment status...")

if (any(!is.na(arg_cpm_meta$treatment_status))) {
  detection_summary <- arg_cpm_meta %>%
    filter(arg_class %in% top_core$arg_class, !is.na(treatment_status)) %>%
    mutate(detected = cpm >= 0.1) %>%
    group_by(arg_class, treatment_status) %>%
    summarise(
      detection_rate = mean(detected, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    filter(n_samples >= 3)
  
  panel_b <- ggplot(detection_summary, 
                    aes(x = reorder(str_trunc(arg_class, 25), detection_rate),
                        y = detection_rate, fill = treatment_status)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("Treated" = p4_pal[3], "Untreated" = p4_pal[6])) +
    coord_flip() +
    labs(
      title = "ARG Detection by Treatment Status",
      subtitle = "Proportion of samples with ≥0.1 CPM",
      x = NULL,
      y = "Detection Rate",
      fill = "Sample Type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
} else {
  # Fallback: show overall detection rates
  detection_overall <- arg_cpm %>%
    filter(arg_class %in% top_core$arg_class) %>%
    group_by(arg_class) %>%
    summarise(
      detection_rate = sum(cpm >= 0.1) / n_distinct(sample),
      n_detections = sum(cpm >= 0.1),
      .groups = "drop"
    ) %>%
    arrange(desc(detection_rate)) %>%
    head(15)
  
  panel_b <- ggplot(detection_overall, 
                    aes(x = reorder(str_trunc(arg_class, 25), detection_rate),
                        y = detection_rate)) +
    geom_col(fill = p4_pal[2], alpha = 0.7) +
    geom_text(aes(label = n_detections), hjust = -0.3, size = 3) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
    coord_flip() +
    labs(
      title = "ARG Detection Frequency",
      subtitle = "Proportion of samples with detection (≥0.1 CPM)",
      x = NULL,
      y = "Detection Rate"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

# Panel C: Read Coverage Distribution
message("Creating Panel C: Read coverage distribution...")

core_cpm_dist <- arg_cpm %>%
  filter(arg_class %in% top_core$arg_class) %>%
  mutate(arg_label = str_trunc(arg_class, 30))

# Get top 8 classes by mean CPM for clarity
top_8_classes <- top_core %>% head(8) %>% pull(arg_class)

panel_c <- core_cpm_dist %>%
  filter(arg_class %in% top_8_classes) %>%
  ggplot(aes(x = reorder(str_trunc(arg_class, 30), cpm, FUN = median),
             y = cpm, fill = arg_class)) +
  geom_violin(alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, show.legend = FALSE) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = rep(p4_pal, length.out = 8)) +
  coord_flip() +
  labs(
    title = "Read Coverage Distribution",
    subtitle = sprintf("Top %d core ARG classes by mean abundance", length(top_8_classes)),
    x = NULL,
    y = "Abundance (CPM, log scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

# Panel D: Treatment Effect (if ΔCLR data available)
message("Creating Panel D: Treatment effects...")

if (!is.null(dclr_data) && nrow(dclr_data) > 0) {
  # Standardize ARG class names for matching
  # abstract3 uses underscores, abstract2 uses spaces/full names
  # Map common patterns
  arg_name_map <- tribble(
    ~abstract3_name, ~abstract2_pattern,
    "aminoglycoside", "aminoglycoside",
    "beta_lactam", "lactam",
    "chloramphenicol", "phenicol",
    "florfenicol", "phenicol",
    "macrolide", "macrolide",
    "multidrug", "disinfecting",
    "polymyxin", "peptide",
    "quinolone", "fluoroquinolone",
    "sulfonamide", "sulfonamide",
    "tetracycline", "tetracycline",
    "trimethoprim", "diaminopyrimidine",
    "vancomycin", "glycopeptide"
  )
  
  # Calculate mean ΔCLR across sites for each ARG class
  dclr_summary <- dclr_data %>%
    group_by(arg_class) %>%
    summarise(
      mean_dclr = mean(delta_clr, na.rm = TRUE),
      median_dclr = median(delta_clr, na.rm = TRUE),
      sd_dclr = sd(delta_clr, na.rm = TRUE),
      n_sites = n(),
      .groups = "drop"
    ) %>%
    mutate(
      se_dclr = sd_dclr / sqrt(n_sites),
      # Try to match with core ARGs
      matches_core = map_lgl(arg_class, ~{
        any(str_detect(tolower(core_args$arg_class), 
                      fixed(tolower(.x))))
      })
    ) %>%
    filter(matches_core | abs(mean_dclr) > 0.5) %>%  # Keep if matches core or has large effect
    arrange(desc(abs(mean_dclr)))
  
  if (nrow(dclr_summary) > 0) {
    top_effects <- dclr_summary %>% head(15)
    
    panel_d <- top_effects %>%
      mutate(
        direction = case_when(
          mean_dclr > 0.5 ~ "Higher in Treated",
          mean_dclr < -0.5 ~ "Higher in Untreated",
          TRUE ~ "Minimal Change"
        ),
        arg_label = str_replace_all(arg_class, "_", " ") %>% str_to_title()
      ) %>%
      ggplot(aes(x = reorder(arg_label, mean_dclr),
                 y = mean_dclr, fill = direction)) +
      geom_col(alpha = 0.8) +
      geom_errorbar(aes(ymin = mean_dclr - se_dclr, 
                       ymax = mean_dclr + se_dclr),
                   width = 0.3, alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", 
                 color = "gray60", alpha = 0.5) +
      scale_fill_manual(
        values = c("Higher in Treated" = p4_pal[3], 
                  "Higher in Untreated" = p4_pal[6],
                  "Minimal Change" = "gray60")
      ) +
      coord_flip() +
      labs(
        title = "Treatment Effects on ARG Classes",
        subtitle = sprintf("ΔCLR (Treated - Untreated) | %d site pairs", 
                          max(dclr_summary$n_sites)),
        x = NULL,
        y = "Mean ΔCLR ± SE",
        fill = "Direction"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom"
      )
  } else {
    panel_d <- NULL
  }
} else {
  panel_d <- NULL
}

# If Panel D couldn't be created, make a summary table instead
if (is.null(panel_d)) {
  message("Creating Panel D fallback: Core ARG summary table...")
  
  summary_text <- top_core %>%
    head(10) %>%
    mutate(
      label = sprintf("%s\n%.1f CPM | %.0f%% prev",
                     str_trunc(arg_class, 35),
                     mean_cpm,
                     prevalence * 100)
    ) %>%
    pull(label) %>%
    paste(collapse = "\n\n")
  
  panel_d <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = summary_text,
             hjust = 0.5, vjust = 0.5, size = 3.5, family = "mono") +
    labs(title = "Top 10 Core ARG Classes",
         subtitle = "Abundance | Prevalence summary") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10))
    )
}

# ========== 5. Combine and Save ==========
message("Combining panels...")

combined_figure <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(
    title = "Integrated Read-Level ARG Analysis",
    subtitle = sprintf("Core ARG detection and treatment effects | %d samples | %d core classes",
                      n_distinct(arg_cpm$sample), nrow(core_args)),
    caption = "Data sources: abstract2.R (ShortBRED quantification) + abstract3.R (compositional analysis)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
      plot.caption = element_text(size = 9, hjust = 0.5, color = "gray40")
    )
  )

# Save figure
output_file <- "combined_reads/integrated_arg_analysis.png"
ggsave(output_file, combined_figure, 
       width = 14, height = 10, dpi = 600, bg = "white")

message(sprintf("\n✓ Combined figure saved to: %s", output_file))

# Save summary statistics
summary_stats <- list(
  n_samples = n_distinct(arg_cpm$sample),
  n_core_args = nrow(core_args),
  n_treated = sum(arg_cpm_meta$treatment_status == "Treated", na.rm = TRUE),
  n_untreated = sum(arg_cpm_meta$treatment_status == "Untreated", na.rm = TRUE),
  top_core_class = core_summary$arg_class[1],
  top_core_mean_cpm = core_summary$mean_cpm[1]
)

write_lines(
  c(
    "# Integrated Read-Level ARG Analysis Summary",
    "",
    sprintf("Total samples: %d", summary_stats$n_samples),
    sprintf("Core ARG classes: %d", summary_stats$n_core_args),
    sprintf("Treated samples: %d", summary_stats$n_treated),
    sprintf("Untreated samples: %d", summary_stats$n_untreated),
    "",
    sprintf("Top core ARG: %s (%.1f CPM)", 
            summary_stats$top_core_class, 
            summary_stats$top_core_mean_cpm),
    "",
    "Output files:",
    sprintf("  - %s", output_file)
  ),
  "combined_reads/ANALYSIS_SUMMARY.txt"
)

message("✓ Analysis complete!")
