#!/usr/bin/env Rscript
# Generate CRISPR host-plasmid ribbon diagram at 600 DPI
# Using spacer_edges.tsv which already has host_lineage taxonomy
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
})

HD_DPI <- 600
OUT <- "HD_poster_figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("\n=== Generating CRISPR Ribbon Diagram (600 DPI) ===\n\n")

# Load CRISPR spacer edges with host_lineage already included
cat("Loading CRISPR spacer edges with host lineage...\n")
spacer <- read_tsv("results/tables/spacer_edges.tsv", show_col_types = FALSE)
cat("  Loaded", nrow(spacer), "spacer edges\n")

# Count how many have taxonomy
n_with_tax <- sum(!is.na(spacer$host_lineage) & spacer$host_lineage != "NA")
cat("  Edges with host taxonomy:", n_with_tax, "/", nrow(spacer), 
    sprintf("(%.1f%%)\n", 100*n_with_tax/nrow(spacer)))

# Load PLSDB taxonomy for plasmid labels
cat("\nLoading PLSDB taxonomy for plasmids...\n")
plsdb_nuccore <- read_csv("PLSDB_db/nuccore.csv", show_col_types = FALSE) %>%
  select(NUCCORE_ACC, TAXONOMY_UID)

plsdb_taxonomy <- read_csv("PLSDB_db/taxonomy.csv", show_col_types = FALSE) %>%
  select(TAXONOMY_UID, TAXONOMY_genus, TAXONOMY_species)

plsdb_tax <- plsdb_nuccore %>%
  left_join(plsdb_taxonomy, by = "TAXONOMY_UID") %>%
  mutate(
    genus = str_to_sentence(TAXONOMY_genus),
    species = str_remove(TAXONOMY_species, "^[^_]+_"),  # Remove genus prefix
    taxa_name = case_when(
      !is.na(species) & species != "" & species != genus ~ paste(genus, species),
      !is.na(genus) & genus != "" ~ paste(genus, "sp."),
      TRUE ~ NA_character_
    )
  ) %>%
  select(plasmid_id = NUCCORE_ACC, genus, taxa_name)

cat("  Loaded taxonomy for", sum(!is.na(plsdb_tax$taxa_name)), "plasmids\n")

# Extract genus and species from host_lineage
# Format: sk:Bacteria;p:Firmicutes;c:Negativicutes;o:Selenomonadales;f:Selenomonadaceae;g:Megamonas;s:Megamonas rupellensis
cat("\nExtracting Linnean nomenclature from host lineages...\n")

extract_taxon <- function(lineage, rank) {
  if (is.na(lineage) || lineage == "NA") return(NA_character_)
  parts <- str_split(lineage, ";")[[1]]
  rank_part <- parts[str_detect(parts, paste0("^", rank, ":"))]
  if (length(rank_part) == 0) return(NA_character_)
  str_remove(rank_part[1], paste0("^", rank, ":"))
}

spacer <- spacer %>%
  mutate(
    genus = map_chr(host_lineage, ~extract_taxon(.x, "g")),
    species_full = map_chr(host_lineage, ~extract_taxon(.x, "s")),
    # Extract just species epithet (remove genus prefix if present)
    species = case_when(
      is.na(species_full) ~ NA_character_,
      str_detect(species_full, " ") ~ str_extract(species_full, " (.+)$") %>% str_trim(),
      TRUE ~ species_full
    ),
    # Create binomial nomenclature
    host_taxa = case_when(
      !is.na(species) & species != "" & species != genus ~ paste(genus, species),
      !is.na(genus) & genus != "" ~ paste(genus, "sp."),
      TRUE ~ "Unbinned"
    ),
    # Flag unbinned contigs
    is_binned = !is.na(genus) & genus != ""
  )

cat("  Created", sum(!is.na(spacer$host_taxa)), "host taxonomy labels\n")
cat("  Binned hosts:", sum(spacer$is_binned), "/ Unbinned:", sum(!spacer$is_binned), "\n")

# Aggregate edges by host-plasmid pairs
cat("\nAggregating spacer edges...\n")
edges_agg <- spacer %>%
  group_by(host_taxa, plasmid_id, is_binned) %>%
  summarise(
    n_spacers = sum(n_unique_spacers),
    n_samples = n_distinct(sample),
    .groups = "drop"
  )

cat("  Aggregated to", nrow(edges_agg), "unique host-plasmid pairs\n")

# Calculate ribbon widths (number of hosts per plasmid)
cat("\nCalculating ribbon widths...\n")
plasmid_diversity <- edges_agg %>%
  group_by(plasmid_id) %>%
  summarise(n_hosts = n_distinct(host_taxa), .groups = "drop")

edges_agg <- edges_agg %>%
  left_join(plasmid_diversity, by = "plasmid_id")

cat("  Ribbon widths range:", min(edges_agg$n_hosts), "-", max(edges_agg$n_hosts), "hosts per plasmid\n")

# Filter to top hosts and plasmids to reduce redundancy
cat("\nFiltering to reduce redundancy...\n")

# Get top binned hosts only (exclude unbinned)
top_binned_hosts <- edges_agg %>%
  filter(is_binned) %>%
  group_by(host_taxa) %>%
  summarise(total_spacers = sum(n_spacers), .groups = "drop") %>%
  arrange(desc(total_spacers)) %>%
  slice_head(n = 35) %>%
  pull(host_taxa)

# Keep unbinned as a single category
top_hosts <- c(top_binned_hosts, "Unbinned")

top_plasmids <- edges_agg %>%
  group_by(plasmid_id) %>%
  summarise(total_spacers = sum(n_spacers), .groups = "drop") %>%
  arrange(desc(total_spacers)) %>%
  slice_head(n = 40) %>%
  pull(plasmid_id)

# Keep only top 10 hosts per plasmid to remove redundancy
edges_filtered <- edges_agg %>%
  filter(host_taxa %in% top_hosts, plasmid_id %in% top_plasmids) %>%
  group_by(plasmid_id) %>%
  arrange(desc(n_spacers)) %>%
  slice_head(n = 10) %>%
  ungroup()

cat("  Filtered to", length(unique(edges_filtered$host_taxa)), "hosts x", 
    length(unique(edges_filtered$plasmid_id)), "plasmids =", 
    nrow(edges_filtered), "interactions\n")
cat("  Binned hosts:", sum(edges_filtered$is_binned %>% unique()), 
    "/ Unbinned:", sum(!edges_filtered$is_binned %>% unique()), "\n")

# Add plasmid taxonomy labels
edges_filtered <- edges_filtered %>%
  left_join(plsdb_tax, by = "plasmid_id") %>%
  mutate(
    plasmid_label = coalesce(taxa_name, plasmid_id)
  )

# Prepare chord diagram data
cat("\nPreparing chord diagram...\n")

# Create adjacency matrix with ribbon widths
chord_data <- edges_filtered %>%
  select(host_taxa, plasmid_label, n_hosts) %>%
  distinct()

# Set up colors - binned hosts in color, unbinned in grey
unique_hosts_final <- unique(edges_filtered$host_taxa)
binned_hosts <- unique_hosts_final[unique_hosts_final != "Unbinned"]
n_binned <- length(binned_hosts)

host_colors <- colorRampPalette(c("#2F6DA8", "#17BEBB", "#2E7D32"))(n_binned)
names(host_colors) <- binned_hosts
if ("Unbinned" %in% unique_hosts_final) {
  host_colors["Unbinned"] <- "#808080"  # Grey for unbinned
}

plasmid_colors <- colorRampPalette(c("#E07A2D", "#D0432B", "#7A3E9D"))(length(top_plasmids))
plasmid_labels <- edges_filtered %>% 
  distinct(plasmid_id, plasmid_label) %>% 
  filter(plasmid_id %in% top_plasmids) %>%
  pull(plasmid_label)
names(plasmid_colors) <- plasmid_labels

all_colors <- c(host_colors, plasmid_colors)

# Generate chord diagram
cat("\nGenerating chord diagram at 600 DPI...\n")
png(file.path(OUT, "crispr_ribbon_diagram_HD.png"), 
    width = 10*HD_DPI, height = 10*HD_DPI, res = HD_DPI)

par(mar = c(2,2,2,2))

# Get unique sectors
unique_hosts <- unique(edges_filtered$host_taxa)
unique_plasmids <- unique(edges_filtered$plasmid_label)
n_hosts_final <- length(unique_hosts)
n_plasmids_final <- length(unique_plasmids)

circos.clear()
circos.par(start.degree = 90)

chordDiagram(
  chord_data,
  grid.col = all_colors,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow",
  link.arr.length = 0.1,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.15)
)

# Add labels with italics for scientific names (skip unbinned)
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  
  # Skip label for unbinned
  if (sector.name == "Unbinned") return()
  
  # Determine if it's a host (has genus in name) or plasmid (accession)
  is_taxa <- str_detect(sector.name, "^[A-Z][a-z]+( |$)")
  font_face <- if(is_taxa) 3 else 1  # 3 = italic for taxa names
  
  circos.text(
    mean(xlim), ylim[1],
    sector.name,
    facing = "clockwise",
    niceFacing = TRUE,
    adj = c(0, 0.5),
    cex = 0.6,
    font = font_face
  )
}, bg.border = NA)

# Add legend
legend(
  "topright",
  legend = c("Hosts (CRISPR-bearing contigs)", "Plasmids (targets)", "Unbinned hosts"),
  fill = c("#2E7D32", "#D0432B", "#808080"),
  border = "black",
  bty = "n",
  cex = 1.0,
  title = "Category"
)

circos.clear()
dev.off()

file_size <- file.info(file.path(OUT, "crispr_ribbon_diagram_HD.png"))$size / 1024^2
cat(sprintf("\n✓ crispr_ribbon_diagram_HD.png (%.1f MB)\n", file_size))
cat("\n=== CRISPR Ribbon Diagram Complete ===\n\n")
