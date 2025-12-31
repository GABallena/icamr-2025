## Figure 3 Sample QC and Purification — combined into one data frame

# Input data (portfolio-safe) -------------------------------------------------
# Expected file format (tab-delimited):
#   sample_code    sample_type    pre_ng_ul    post_ng_ul
#
# If the file is missing, the script generates a small demo dataset so that the
# plotting code can be run without any project-specific values.
#
# Override path:
#   QUBIT_TSV=/path/to/qubit_pre_post.tsv Rscript qc_purification_and_depth_arg_figures.R

qubit_path <- Sys.getenv("QUBIT_TSV", "data/qubit_pre_post.tsv")

if (file.exists(qubit_path)) {
  qubit_tbl <- read.delim(qubit_path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("sample_code","sample_type","pre_ng_ul","post_ng_ul") %in% names(qubit_tbl)))
  sample_code    <- as.character(qubit_tbl$sample_code)
  sample_type    <- as.character(qubit_tbl$sample_type)
  pre_qubit_raw  <- as.character(qubit_tbl$pre_ng_ul)
  post_qubit_raw <- as.character(qubit_tbl$post_ng_ul)
} else {
  set.seed(42)
  n <- 12
  sample_code <- sprintf("SAMPLE-%02d", seq_len(n))
  sample_type <- rep(c("GroupA","GroupB","GroupC","Control"), length.out = n)
  # keep as *strings* to exercise the parser below (including '<' cases)
  pre_qubit_raw  <- sprintf("%.3f", runif(n, 0.05, 5))
  post_qubit_raw <- sprintf("%.2f",  runif(n, 1, 100))
  pre_qubit_raw[sample(1:n, 2)] <- "< 0.10"
  message("NOTE: QUBIT_TSV not found; using demo values (no project data).")
}

# Helpers --------------------------------------------------------------------
parse_qubit <- function(x) {
	# Returns a numeric value with '<' removed and a logical flag for below detection
	below <- grepl("^\\s*<", x)
	# strip any leading '<' and spaces
	num <- as.numeric(gsub("^\\s*<\\s*", "", x))
	list(value = num, below_dl = below)
}

pre_parsed <- parse_qubit(pre_qubit_raw)
post_parsed <- parse_qubit(post_qubit_raw)

# Final merged data frame -----------------------------------------------------
qc_df <- data.frame(
	sample_code = sample_code,
	sample_type = sample_type,
	qubit_pre_raw = pre_qubit_raw,
	qubit_pre = pre_parsed$value,
	qubit_pre_below_dl = pre_parsed$below_dl,
	qubit_post = post_parsed$value,
	stringsAsFactors = FALSE
)

# Optional: order columns nicely and show a quick preview
qc_df <- qc_df[, c("sample_code","sample_type","qubit_pre_raw","qubit_pre_below_dl","qubit_pre","qubit_post")]
print(qc_df)

#### Qubit plot in the style of the template (two lines + difference shading) ####

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Add a numeric index per sample to allow a continuous x scale (required by stat_difference)
qc_df <- qc_df %>%
  mutate(
    sample_num = as.integer(gsub("^SAMPLE-", "", sample_code)),
    .sample_id = row_number()
  )

# Choose colors similar to template
col_pre  <- "#C32E5A"  # red-ish
col_post <- "#3D85F7"  # blue-ish

has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
if (!has_ggh4x) {
  # Try to install from CRAN to enable filled intersections
  try({
    install.packages("ggh4x", repos = "https://cloud.r-project.org", quiet = TRUE)
  }, silent = TRUE)
  has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
}
if (has_ggh4x) {
  library(ggh4x)
}
has_colorspace <- requireNamespace("colorspace", quietly = TRUE)
lighten_or_same <- function(col) {
  if (has_colorspace) return(colorspace::lighten(col))
  col
}

# Build plot
plt <- ggplot(qc_df, aes(x = .sample_id)) +
  geom_line(aes(y = qubit_post, color = "post"), linewidth = 0.7) +
  geom_line(aes(y = qubit_pre,  color = "pre"),  linewidth = 0.7)

if (has_ggh4x) {
  plt <- plt +
    ggh4x::stat_difference(
      aes(ymin = qubit_pre, ymax = qubit_post, fill = after_stat(sign)),
      alpha = 0.3
    )
}


plt <- plt +
  scale_color_manual(
    values = c(pre = col_pre, post = col_post),
    breaks = c("pre", "post"),
    labels = c("Pre-extraction", "Post-purification")
  ) +
  # Fill colours for stat_difference: order corresponds to (top, bottom, equal)
  (if (has_ggh4x) scale_fill_manual(
    values = c(
      lighten_or_same(col_post),
      lighten_or_same(col_pre),
      "grey70"
    ),
    labels = c("post > pre", "pre > post", "same"),
    name = NULL
  ) else NULL) +
  scale_x_continuous(
    breaks = qc_df$.sample_id,
    labels = qc_df$sample_code,
    expand = expansion(add = c(0.2, 0.2))
  ) +
  labs(
    title = "Qubit concentrations: pre vs post purification",
    y = "ng/µL",
    x = NULL,
    caption = NULL
  ) +
  guides(color = guide_legend(order = 1), fill = "none") +
  theme_minimal() +
  theme(
  legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.title = element_blank(),
  plot.background = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 16, 12, 16),
    plot.title = element_text(size = 16, face = "bold", color = "grey25"),
    plot.caption = element_text(size = 9),
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "grey40", angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.text.y = element_text(color = "grey40"),
    strip.text = element_text(face = "bold", color = "grey20")
  )

print(plt)

# Save a copy
try({
  ggsave(filename = "qubit_pre_post_template_style.png", plot = plt, width = 12, height = 7, dpi = 300, bg="white")
}, silent = TRUE)



### Figure 6 Genus and ARG Identification
# Library insert sizes (bp) ---------------------------------------------------
# Optional file format (tab-delimited):
#   sample_code    insert_bp
#
# Override path:
#   LIB_INSERT_TSV=/path/to/library_insert_bp.tsv
lib_path <- Sys.getenv("LIB_INSERT_TSV", "data/library_insert_bp.tsv")

if (file.exists(lib_path)) {
  lib_tbl <- read.delim(lib_path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c("sample_code","insert_bp") %in% names(lib_tbl)))
  library_bp <- as.numeric(lib_tbl$insert_bp[match(qc_df$sample_code, lib_tbl$sample_code)])
} else {
  # demo values (no project data); adjust to your dataset for real runs
  set.seed(7)
  library_bp <- round(runif(nrow(qc_df), min = 250, max = 650))
  message("NOTE: LIB_INSERT_TSV not found; using demo insert sizes (no project data).")
}


# Attach to qc_df for convenience and compute per-kb factor
stopifnot(length(library_bp) == nrow(qc_df))
qc_df$library_bp <- library_bp
qc_df$library_kb <- qc_df$library_bp / 1000

# Helper to normalize a numeric vector by library size (counts per kb of library)
normalize_by_library <- function(x, lib_kb) {
  # Returns counts per kb of library length
  x / lib_kb
}

# Example (placeholder): when loading ARG counts for Figure 6, prefer unnormalized counts and
# convert to library-size-normalized counts instead of using RPKM.
# Pseudocode sketch:
#   counts_vec <- c(...)  # total counts per feature for a given sample
#   lib_norm   <- normalize_by_library(counts_vec, qc_df$library_kb[idx_of_sample])
# Then proceed to plotting/aggregation using lib_norm.


# === Figure 6 actual panel: depth vs ARG richness with bubble size = ARG load (CPM) ===

arg_samples_dir <- file.path(getwd(), "arg_oap_work", "samples")
sample_counts_path <- file.path(getwd(), "diversity_results", "sample_counts.tsv")

read_arg_type_counts <- function(sdir) {
  # Prefer merged/ tables; fallback to base sample dir; fallback to gene-level if type-level missing
  candidates <- c(
    file.path(sdir, "merged", "unnormalized_count.type.txt"),
    file.path(sdir, "unnormalized_count.type.txt")
  )
  fp <- candidates[file.exists(candidates)][1]
  fallback_gene <- NULL
  if (is.na(fp) || is.null(fp)) {
    # Try gene-level as last resort
    cand_g <- c(
      file.path(sdir, "merged", "unnormalized_count.gene.txt"),
      file.path(sdir, "unnormalized_count.gene.txt")
    )
    fallback_gene <- cand_g[file.exists(cand_g)][1]
  }
  if (!is.na(fp) && length(fp) == 1 && !is.null(fp)) {
    df <- tryCatch(suppressMessages(read.delim(fp, check.names = FALSE)), error = function(e) NULL)
    if (!is.null(df) && ncol(df) >= 3) {
      num_cols <- grep("_R[12]_trimmed", names(df))
      if (length(num_cols) < 1) num_cols <- 2:ncol(df)  # fallback if headers are simpler
      df$sum_counts <- rowSums(df[, num_cols, drop = FALSE], na.rm = TRUE)
      total <- sum(df$sum_counts, na.rm = TRUE)
      richness <- sum(df$sum_counts > 0, na.rm = TRUE)
      return(data.frame(total_ARG_counts = total, arg_richness_types = richness))
    }
  }
  if (!is.na(fallback_gene) && length(fallback_gene) == 1 && !is.null(fallback_gene)) {
    df <- tryCatch(suppressMessages(read.delim(fallback_gene, check.names = FALSE)), error = function(e) NULL)
    if (!is.null(df) && ncol(df) >= 3) {
      num_cols <- grep("_R[12]_trimmed", names(df))
      if (length(num_cols) < 1) num_cols <- 2:ncol(df)
      df$sum_counts <- rowSums(df[, num_cols, drop = FALSE], na.rm = TRUE)
      total <- sum(df$sum_counts, na.rm = TRUE)
      richness <- sum(df$sum_counts > 0, na.rm = TRUE)
      return(data.frame(total_ARG_counts = total, arg_richness_types = richness))
    }
  }
  return(NULL)
}

list_dirs <- function(path) {
  xs <- list.dirs(path, full.names = TRUE, recursive = FALSE)
  # Keep only base sample dirs (exclude .part_ and non-sample folders)
  xs[grepl("/SAMPLE-.*(_S|_RUN)[0-9]+$", xs)]
}

metrics <- NULL
if (dir.exists(arg_samples_dir)) {
  sdirs <- list_dirs(arg_samples_dir)
  rows <- lapply(sdirs, function(sd) {
    mm <- read_arg_type_counts(sd)
    if (is.null(mm)) return(NULL)
    sample_id <- basename(sd)
    cbind(data.frame(sample = sample_id, stringsAsFactors = FALSE), mm)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) > 0) metrics <- do.call(rbind, rows)
}

if (!is.null(metrics)) {
  # Attach sequencing depth N (reads)
  if (file.exists(sample_counts_path)) {
    sc <- suppressMessages(read.delim(sample_counts_path, check.names = FALSE))
    # Normalize columns: lowercase then rename 'n' -> 'N'
    names(sc) <- tolower(names(sc))
    # try to identify sample column
    if (!"sample" %in% names(sc)) names(sc)[1] <- "sample"
    # identify N column and rename to 'N' for downstream
    n_col <- if ("n" %in% names(sc)) "n" else if ("N" %in% names(sc)) "N" else names(sc)[2]
    names(sc)[names(sc) == n_col] <- "N"
    sc$sample <- as.character(sc$sample)
    metrics <- merge(metrics, sc[, c("sample","N")], by = "sample", all.x = TRUE)
  } else {
    metrics$N <- NA_real_
  }

  # Map sample_type from qc_df via prefix before underscore (SAMPLE-XX)
  metrics$sample_prefix <- sub("_.*$", "", metrics$sample)
  st_map <- qc_df[, c("sample_code", "sample_type", "library_kb")]
  names(st_map)[1] <- "sample_prefix"
  metrics <- merge(metrics, st_map, by = "sample_prefix", all.x = TRUE)

  # Compute log10 depth and ARG load CPM (per million reads)
  metrics$log10N <- ifelse(is.na(metrics$N) | metrics$N <= 0, NA_real_, log10(metrics$N))
  metrics$arg_cpm <- ifelse(is.na(metrics$N) | metrics$N <= 0, NA_real_, (metrics$total_ARG_counts / metrics$N) * 1e6)

  # Bubble plot
  suppressPackageStartupMessages({
    has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
  })
  p6 <- ggplot(metrics, aes(x = log10N, y = arg_richness_types)) +
    geom_hline(yintercept = 0, color = "grey92") +
    geom_vline(xintercept = 6, linetype = "dashed", color = "grey70") +
    geom_point(aes(size = arg_cpm, color = sample_type), alpha = 0.85) +
    scale_size_continuous(name = "ARG load (CPM)", range = c(2, 10)) +
    scale_color_manual(values = c(
      "H-untreated" = "#E64B35",
      "H-treated" = "#4DBBD5",
      "CW-untreated" = "#00A087",
      "CW-treated" = "#3C5488",
      "SW-freshwater" = "#F39B7F",
      "SW-marine water" = "#7E6148",
      "Control" = "#91D1C2",
      "Optimization" = "#8491B4"
    ), na.value = "#8A8A8A", guide = guide_legend(override.aes = list(size = 5))) +
    labs(
      title = "Figure 6: Depth vs ARG richness (bubble size = ARG load per million reads)",
      x = expression(log[10](Sequencing~depth~(reads))),
      y = "ARG types detected"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  # Label a few interesting points: lowest depth and highest ARG CPM
  if (exists("has_ggrepel") && isTRUE(has_ggrepel)) {
    library(ggrepel)
    lab_df <- metrics
    # pick 3 lowest depth and 3 highest CPM
    idx_lowN <- order(metrics$N, na.last = NA)[seq_len(min(3, sum(!is.na(metrics$N))))]
    idx_hiCPM <- order(-metrics$arg_cpm, na.last = NA)[seq_len(min(3, sum(!is.na(metrics$arg_cpm))))]
    sel <- sort(unique(c(idx_lowN, idx_hiCPM)))
    lab_df <- metrics[sel, , drop = FALSE]
    p6 <- p6 + ggrepel::geom_text_repel(data = lab_df, aes(label = sample), size = 3, seed = 42, max.overlaps = 50)
  }

  print(p6)
  try({
    ggsave(filename = "figure6_arg_depth_bubble.png", plot = p6, width = 12, height = 7, dpi = 300, bg = "white")
  }, silent = TRUE)
}
