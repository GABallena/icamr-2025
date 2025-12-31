#!/usr/bin/env Rscript
# Generate HD (600 DPI) versions of key figures for poster
# Target figures:
# 1. figs/Fig1_ARG_prevalence_abundance.png
# 2. figs/Fig7_prevalence_methods.png
# 3. mag_analysis/mag_taxonomy_sunburst.png
# 4. mag_analysis/arg_viz_strict/mag_drug_class_stack.png
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
  library(viridis)
  library(scales)
  library(ggrepel)
})

HD_DPI <- 600
OUT_DIR <- "HD_poster_figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n\n")
cat("   Generating HD Poster Figures (600 DPI)                \n")
cat("\n\n")

# Custom theme
theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.8),
      strip.text = element_text(face = "bold", size = rel(0.9)),
      axis.title = element_text(face = "bold", size = rel(1.1)),
      axis.text = element_text(color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0),
      plot.subtitle = element_text(size = rel(0.95), color = "grey30", hjust = 0),
      plot.caption = element_text(size = rel(0.8), color = "grey50", hjust = 1)
    )
}

# Helper function
simplify_arg_for_viz <- function(df, col_name = "arg_class") {
  df %>%
    mutate(
      arg_class_display = case_when(
        str_detect(!!sym(col_name), ";") ~ "Multi-drug resistance",
        TRUE ~ !!sym(col_name)
      )
    )
}

# Load data
cat("[1/4] Loading ARG data...\n")
prev_over_X <- read_tsv("results/class_prevalence_over_X.tsv", show_col_types = FALSE)
prev_LOD <- read_tsv("results/class_prevalence_LOD.tsv", show_col_types = FALSE)

# Apply simplification
prev_over_X <- simplify_arg_for_viz(prev_over_X) %>%
  group_by(arg_class_display) %>%
  summarise(
    prevalence_over_X = max(prevalence_over_X),
    mean_CPM = sum(mean_CPM),
    .groups = "drop"
  ) %>%
  rename(arg_class = arg_class_display)

prev_LOD <- simplify_arg_for_viz(prev_LOD) %>%
  group_by(arg_class_display) %>%
  summarise(prevalence_LOD = max(prevalence_LOD), .groups = "drop") %>%
  rename(arg_class = arg_class_display)

core_prev_threshold <- 0.70
cpm_core_threshold <- 0.1
lod_min_reads <- 3

cat("[2/4] Generating Figure 1: ARG Prevalence & Abundance (HD)...\n")
fig1_data <- prev_over_X %>%
  arrange(desc(prevalence_over_X)) %>%
  slice(1:20) %>%
  mutate(arg_class = forcats::fct_reorder(arg_class, prevalence_over_X),
         is_core = prevalence_over_X >= core_prev_threshold)

p1a <- ggplot(fig1_data, aes(x = prevalence_over_X, y = arg_class, fill = is_core)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = core_prev_threshold, linetype = "dashed", 
             color = "red", linewidth = 0.8, alpha = 0.7) +
  scale_x_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = c("TRUE" = "#2E7D32", "FALSE" = "#757575"),
                    labels = c("TRUE" = "Core", "FALSE" = "Non-core"),
                    name = "Classification") +
  labs(x = "Prevalence (%)", y = NULL,
       title = "Top 20 ARG Classes by Prevalence",
       subtitle = glue::glue("Core threshold: ≥{round(100*core_prev_threshold)}% prevalence at ≥{cpm_core_threshold} CPM")) +
  theme_publication(base_size = 14) +
  theme(legend.position.inside = c(0.85, 0.15))

p1b <- ggplot(fig1_data, aes(x = mean_CPM, y = arg_class, fill = log10(mean_CPM + 1))) +
  geom_col(width = 0.7) +
  scale_x_log10(labels = comma_format()) +
  scale_fill_viridis_c(option = "plasma", name = "log10(CPM)", guide = "none") +
  labs(x = "Mean CPM (log scale)", y = NULL,
       title = "Mean Abundance") +
  theme_publication(base_size = 14)

fig1 <- p1a + p1b + plot_layout(widths = c(2, 1))
ggsave(file.path(OUT_DIR, "Fig1_ARG_prevalence_abundance_HD.png"), 
       fig1, width = 12, height = 7, dpi = HD_DPI, bg = "white")
cat("  ✓ Saved: Fig1_ARG_prevalence_abundance_HD.png\n\n")

cat("[3/4] Generating Figure 7: Prevalence Methods Comparison (HD)...\n")
prev_grouped <- prev_over_X %>%
  select(arg_class, prevalence_over_X)

prev_lod_grouped <- prev_LOD %>%
  select(arg_class, prevalence_LOD)

comp_prev <- prev_grouped %>%
  left_join(prev_lod_grouped, by = "arg_class") %>%
  filter(!is.na(prevalence_over_X), !is.na(prevalence_LOD)) %>%
  rename(prev_CPM = prevalence_over_X, prev_LOD = prevalence_LOD) %>%
  mutate(is_core_CPM = prev_CPM >= core_prev_threshold,
         is_core_LOD = prev_LOD >= core_prev_threshold,
         category = case_when(
           is_core_CPM & is_core_LOD ~ "Core (both methods)",
           is_core_CPM & !is_core_LOD ~ "Core (CPM only)",
           !is_core_CPM & is_core_LOD ~ "Core (LOD only)",
           TRUE ~ "Non-core"
         ))

labels_to_show <- comp_prev %>%
  filter(is_core_CPM | is_core_LOD | prev_CPM >= 0.5 | prev_LOD >= 0.5) %>%
  distinct(arg_class, .keep_all = TRUE)

p7 <- ggplot(comp_prev, aes(x = prev_CPM, y = prev_LOD, color = category)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = core_prev_threshold, linetype = "dotted", color = "red", alpha = 0.5) +
  geom_vline(xintercept = core_prev_threshold, linetype = "dotted", color = "red", alpha = 0.5) +
  geom_point(size = 3.5, alpha = 0.7) +
  ggrepel::geom_label_repel(
    data = labels_to_show,
    aes(label = arg_class),
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    max.overlaps = 20,
    min.segment.length = 0.1,
    force = 2,
    seed = 42
  ) +
  scale_x_continuous(labels = percent_format()) +
  scale_y_continuous(labels = percent_format()) +
  scale_color_manual(values = c("Core (both methods)" = "#1B5E20",
                                "Core (CPM only)" = "#FFA726",
                                "Core (LOD only)" = "#42A5F5",
                                "Non-core" = "#9E9E9E"),
                    name = "Classification") +
  labs(x = "Prevalence by CPM Method", y = "Prevalence by LOD Method",
       title = "ARG Prevalence: CPM vs. LOD Detection Methods",
       subtitle = "Multi-drug resistances grouped; core classes and high-prevalence ARGs labeled",
       caption = glue::glue("Core threshold: {round(100*core_prev_threshold)}% | LOD: ≥{lod_min_reads} reads | n={nrow(comp_prev)} classes")) +
  theme_publication(base_size = 14) +
  theme(legend.position.inside = c(0.2, 0.75),
        legend.background = element_rect(fill = "white", color = "black"))

ggsave(file.path(OUT_DIR, "Fig7_prevalence_methods_HD.png"), 
       p7, width = 11, height = 9, dpi = HD_DPI, bg = "white")
cat("  ✓ Saved: Fig7_prevalence_methods_HD.png\n\n")

cat("[4/4] Checking for MAG analysis figures...\n")

# For MAG figures, we need to check if they exist and if the source scripts are available
mag_sunburst_exists <- file.exists("mag_analysis/mag_taxonomy_sunburst.png")
mag_drug_stack_exists <- file.exists("mag_analysis/arg_viz_strict/mag_drug_class_stack.png")

if(mag_sunburst_exists || mag_drug_stack_exists) {
  cat("  -> MAG figure sources found. Checking for generation scripts...\n")
  
  # Check for abstract6.R or similar MAG analysis scripts
  if(file.exists("abstract6.R")) {
    cat("  -> abstract6.R found. Note: MAG figures require running full MAG analysis.\n")
    cat("  -> To regenerate MAG HD figures, please run the generate_remaining_hd_figures.R script\n")
    cat("     or source the relevant sections of abstract6.R manually.\n")
  } else {
    cat("  -> MAG analysis script not found in current directory.\n")
    cat("  -> MAG HD figures will need to be generated separately.\n")
  }
} else {
  cat("  -> MAG figures not found. Skipping MAG HD generation.\n")
  cat("  -> If you need MAG figures, please run the MAG analysis pipeline first.\n")
}

cat("\n")
cat("\n")
cat("   HD Poster Figure Generation Complete                  \n")
cat("\n")
cat("  Location: HD_poster_figures/                              \n")
cat("  Resolution: 600 DPI                                       \n")
cat("\n\n")

# List generated files
hd_files <- list.files(OUT_DIR, pattern = "_HD\\.png$", full.names = FALSE)
cat("Generated HD figures:\n")
for(f in hd_files) {
  fpath <- file.path(OUT_DIR, f)
  fsize <- file.info(fpath)$size / 1024^2
  cat(sprintf("  ✓ %-50s [%.1f MB]\n", f, fsize))
}

cat("\nNote: For MAG analysis HD figures (sunburst and drug class stack),\n")
cat("      please ensure abstract6.R has been run or use the dedicated\n")
cat("      MAG HD figure generation script if available.\n\n")
