#!/usr/bin/env Rscript
# Generate HD CRISPR Ribbon Diagram (600 DPI)
# Simplified version with ribbon width = number of taxa
# Removes redundancies from the ultra circos version
# -------------------------------------------------------------------------
# Portfolio-ready note:
# - No hard-coded local paths.
# - Run from the project root (or set PROJECT_ROOT env var).
# -------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = getwd())
try(setwd(PROJECT_ROOT), silent = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(circlize)
  library(scales)
})

HD_DPI <- 600
OUT_DIR <- "HD_poster_figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n\n")
cat("   CRISPR Host-Plasmid Ribbon Diagram (600 DPI)          \n")
cat("\n\n")

# CityCore palette
citystat_core12 <- c(
  Navy      = "#233B5D", Blue      = "#2F6DA8", DeepTeal  = "#1C7C7D", Cyan      = "#17BEBB",
  Forest    = "#2E7D32", Olive     = "#7A8F30", Mustard   = "#C7A41A", Orange    = "#E07A2D",
  Vermilion = "#D0432B", Crimson   = "#B01E2F", Plum      = "#7A3E9D", Indigo    = "#3B4BA3"
)

# Read data
cat("[1/4] Loading CRISPR spacer data...\n")
spacer_edges <- read_tsv("results/tables/spacer_edges_contig.tsv", show_col_types = FALSE)

# Load taxonomy mappings
contig2bin <- if(file.exists("mag_analysis/full_contig_bin_map.tsv")) {
  read_tsv("mag_analysis/full_contig_bin_map.tsv", show_col_types = FALSE) %>%
    mutate(contig_full = paste0(sample, "_", contig)) %>%
    select(contig_full, mag_id)
} else {
  tibble(contig_full = character(), mag_id = character())
}

# Load PLSDB taxonomy for plasmid labels
plsdb_nuccore <- if(file.exists("PLSDB_db/nuccore.csv")) {
  read_csv("PLSDB_db/nuccore.csv", show_col_types = FALSE) %>%
    select(NUCCORE_ACC, TAXONOMY_UID)
} else {
  tibble(NUCCORE_ACC = character(), TAXONOMY_UID = character())
}

plsdb_taxonomy <- if(file.exists("PLSDB_db/taxonomy.csv")) {
  read_csv("PLSDB_db/taxonomy.csv", show_col_types = FALSE) %>%
    select(TAXONOMY_UID, TAXONOMY_genus, TAXONOMY_species) %>%
    mutate(
      plasmid_taxa = case_when(
        !is.na(TAXONOMY_species) ~ str_replace(TAXONOMY_species, "_", " "),
        !is.na(TAXONOMY_genus) ~ paste(TAXONOMY_genus, "sp."),
        TRUE ~ NA_character_
      )
    )
} else {
  tibble(TAXONOMY_UID = character(), plasmid_taxa = character())
}

# Create plasmid accession to taxa mapping
plasmid2taxa <- plsdb_nuccore %>%
  left_join(plsdb_taxonomy, by = "TAXONOMY_UID") %>%
  mutate(
    acc_base = sub("\\.\\d+$", "", NUCCORE_ACC)
  ) %>%
  select(acc_base, plasmid_taxa) %>%
  distinct()

gtdbtk_files <- list.files("gtdbtk_drep_out", pattern = "gtdbtk.bac120.summary.tsv$", 
                           recursive = TRUE, full.names = TRUE)

gtdbtk_tax <- if(length(gtdbtk_files) > 0) {
  map_dfr(gtdbtk_files, ~{
    df <- read_tsv(.x, show_col_types = FALSE)
    df %>% select(user_genome, classification)
  }) %>%
    distinct() %>%
    mutate(
      # Extract full taxonomic path
      genus = str_match(classification, ";g__([^;]+)")[,2],
      species = str_match(classification, ";s__([^;]+)")[,2],
      genus = str_replace(genus, "_[A-Z]$", ""),
      species = str_replace(species, "^[^_]+_", ""),  # Remove genus prefix
      species = str_replace(species, "_[A-Z]$", ""),
      # Filter out placeholder IDs
      genus = ifelse(str_detect(genus, "^[A-Z0-9]+$"), NA_character_, genus),
      species = ifelse(str_detect(species, "^[A-Z0-9]+$"), NA_character_, species),
      # Create proper binomial nomenclature
      taxa_name = case_when(
        !is.na(genus) & !is.na(species) & species != "" ~ paste(genus, species),
        !is.na(genus) ~ paste(genus, "sp."),
        TRUE ~ "Unknown"
      ),
      mag_id = user_genome
    ) %>%
    select(mag_id, genus, taxa_name)
} else {
  tibble(mag_id = character(), genus = character(), taxa_name = character())
}

cat("[2/4] Processing and aggregating edges...\n")

# Create contig to genus mapping
contig2genus <- contig2bin %>%
  left_join(gtdbtk_tax, by = "mag_id") %>%
  mutate(
    genus = coalesce(genus, "Unknown"),
    taxa_name = coalesce(taxa_name, "Unknown")
  ) %>%
  select(contig_full, mag_id, genus, taxa_name) %>%
  distinct()

# Debug: check mapping coverage
cat(sprintf("  -> Loaded %d contig-to-MAG mappings\n", nrow(contig2genus)))
cat(sprintf("  -> MAGs with taxonomy: %d\n", sum(contig2genus$taxa_name != "Unknown")))

# Aggregate edges and get genus counts
edges_agg <- spacer_edges %>%
  group_by(contig, plasmid_id) %>%
  summarise(
    n_spacers = sum(weight),
    .groups = "drop"
  ) %>%
  # Join with genus info - contig IDs should match directly
  left_join(contig2genus, by = c("contig" = "contig_full")) %>%
  mutate(
    genus = coalesce(genus, "Unknown"),
    taxa_name = coalesce(taxa_name, "Unknown"),
    mag_id = coalesce(mag_id, "Unbinned")
  )

# Debug: check join success
n_mapped <- sum(edges_agg$taxa_name != "Unknown")
cat(sprintf("  -> Contigs successfully mapped to taxonomy: %d/%d (%.1f%%)\n", 
            n_mapped, nrow(edges_agg), 100*n_mapped/nrow(edges_agg)))

# Select top hosts and plasmids by interaction strength
top_n <- 40  # Increased for better coverage
top_hosts <- edges_agg %>%
  count(contig, wt = n_spacers, sort = TRUE) %>%
  slice_head(n = top_n) %>%
  pull(contig)

top_plasmids <- edges_agg %>%
  count(plasmid_id, wt = n_spacers, sort = TRUE) %>%
  slice_head(n = top_n) %>%
  pull(plasmid_id)

# Filter to top interactions and count unique genera per plasmid
edges_filtered <- edges_agg %>%
  filter(contig %in% top_hosts, plasmid_id %in% top_plasmids)

# Calculate ribbon width = number of distinct genera targeting each plasmid
ribbon_widths <- edges_filtered %>%
  group_by(plasmid_id) %>%
  summarise(
    n_hosts = n_distinct(contig),  # Number of different host contigs
    n_genera = n_distinct(genus[genus != "Unknown"]),
    total_spacers = sum(n_spacers),
    .groups = "drop"
  ) %>%
  mutate(
    # Use n_hosts as ribbon width (more reliable than genera)
    ribbon_width = pmax(1, n_hosts)
  )

# Join back to edges
edges_with_taxa <- edges_filtered %>%
  left_join(ribbon_widths %>% select(plasmid_id, ribbon_width, n_hosts), 
            by = "plasmid_id") %>%
  arrange(desc(ribbon_width), desc(n_spacers))

cat(sprintf("  -> Filtered to %d hosts × %d plasmids = %d interactions\n", 
    length(top_hosts), length(top_plasmids), nrow(edges_filtered)))
cat(sprintf("  -> Ribbon widths (hosts/plasmid) range from %d to %d\n", 
    min(edges_with_taxa$ribbon_width), max(edges_with_taxa$ribbon_width)))

# Remove redundant edges: for each plasmid, keep only top N hosts by spacer count
max_hosts_per_plasmid <- 10
edges_dedup <- edges_with_taxa %>%
  group_by(plasmid_id) %>%
  slice_max(n_spacers, n = max_hosts_per_plasmid, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("  -> Removed redundancies: %d -> %d edges\n", 
    nrow(edges_with_taxa), nrow(edges_dedup)))

# Prepare chord diagram data
cat("[3/4] Preparing chord diagram matrix...\n")

# Create adjacency matrix
all_hosts <- unique(edges_dedup$contig)
all_plasmids <- unique(edges_dedup$plasmid_id)
all_nodes <- c(all_hosts, all_plasmids)

mat <- matrix(0, nrow = length(all_nodes), ncol = length(all_nodes))
rownames(mat) <- colnames(mat) <- all_nodes

# Fill matrix: use ribbon_width as the value (# of genera)
for (i in seq_len(nrow(edges_dedup))) {
  h <- edges_dedup$contig[i]
  p <- edges_dedup$plasmid_id[i]
  # Use ribbon_width (n_genera) as the value
  mat[h, p] <- edges_dedup$ribbon_width[i]
}

# Assign colors: hosts = forest green, plasmids = plum purple
node_colors <- ifelse(all_nodes %in% all_hosts, 
                     citystat_core12[["Forest"]], 
                     citystat_core12[["Plum"]])
names(node_colors) <- all_nodes

cat("[4/4] Generating ribbon diagram...\n")

# Create node labels: genus names for hosts, plasmid IDs for plasmids
node_labels <- setNames(rep("", length(all_nodes)), all_nodes)

# For hosts, use MAG taxonomy names (from bins) or extract sample info
host_labels <- edges_dedup %>%
  filter(contig %in% all_hosts) %>%
  select(contig, taxa_name, mag_id) %>%
  distinct() %>%
  mutate(
    # Extract sample and contig number for unbinned contigs
    sample_id = sub("^([A-Za-z0-9]+-\\d+)_.*", "\\1", contig),
    contig_num = sub(".*k141_(\\d+).*", "\\1", contig),
    # Use taxa name if available, otherwise create a readable label
    label = case_when(
      taxa_name != "Unknown" ~ taxa_name,
      mag_id != "Unbinned" ~ paste0("MAG: ", sub(".*bin\\.(\\d+).*", "Bin\\1", mag_id)),
      TRUE ~ paste0(sample_id, " Ctg", contig_num)  # e.g., "PROJECT-12 Ctg58020"
    )
  )

for (i in seq_len(nrow(host_labels))) {
  node_labels[host_labels$contig[i]] <- host_labels$label[i]
}

# For plasmids, use taxa from PLSDB
plasmid_labels <- tibble(plasmid_id = all_plasmids) %>%
  mutate(
    acc_base = sub("\\.\\d+$", "", plasmid_id)
  ) %>%
  left_join(plasmid2taxa, by = "acc_base") %>%
  mutate(
    label = coalesce(
      plasmid_taxa,
      sub("^([A-Z]+_?[0-9]+).*", "\\1", plasmid_id)
    )
  )

for (i in seq_len(nrow(plasmid_labels))) {
  node_labels[plasmid_labels$plasmid_id[i]] <- plasmid_labels$label[i]
}

# Generate high-res plot
png(file.path(OUT_DIR, "crispr_ribbon_diagram_HD.png"), 
    width = 10, height = 10, units = "in", res = HD_DPI, bg = "white")

# Set circlize parameters for cleaner layout
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 2,
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.01, 0.01)
)

# Create chord diagram with custom labels
chordDiagram(
  mat,
  grid.col = node_colors,
  transparency = 0.3,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow",
  diffHeight = 0.04,
  annotationTrack = "grid",
  annotationTrackHeight = 0.03,
  link.border = alpha("grey30", 0.5),
  link.lwd = 0.5,
  preAllocateTracks = list(
    track.height = 0.1
  )
)

# Add custom labels with proper formatting
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  label = node_labels[sector.name]
  if (!is.na(label) && label != "") {
    # Format text: both hosts and plasmids in italic (both are Linnean taxa)
    is_host <- sector.name %in% all_hosts
    circos.text(
      CELL_META$xcenter, 
      CELL_META$ylim[1], 
      label,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.45,
      font = 3  # italic for all taxonomic names
    )
  }
}, bg.border = NA)

# Add title and legend
title("CRISPR Host-Plasmid Interactions", cex.main = 1.5, font.main = 2)

# Add legend
legend("bottomleft", 
       legend = c("Host bacteria (CRISPR)", "Plasmid (targeted)", 
                  "Ribbon width ∝ host diversity"),
       fill = c(citystat_core12[["Forest"]], citystat_core12[["Plum"]], NA),
       border = c("black", "black", NA),
       bty = "n",
       cex = 0.85)

# Add statistics text
total_interactions <- sum(mat > 0)
unique_hosts <- length(all_hosts)
unique_plasmids <- length(all_plasmids)
text(-0.9, -0.9, 
     sprintf("Hosts: %d | Plasmids: %d\nInteractions: %d", 
             unique_hosts, unique_plasmids, total_interactions),
     adj = c(0, 0), cex = 0.8, col = "grey30")

circos.clear()
dev.off()

cat(sprintf("  ✓ Saved: crispr_ribbon_diagram_HD.png [%.1f MB]\n",
    file.size(file.path(OUT_DIR, "crispr_ribbon_diagram_HD.png")) / 1024^2))

cat("\n\n")
cat("   CRISPR Ribbon Diagram Complete                        \n")
cat("\n")
cat("  Ribbon width represents taxonomic diversity              \n")
cat("  Redundant interactions removed                           \n")
cat("  Resolution: 600 DPI                                      \n")
cat("\n\n")
