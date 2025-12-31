#!/usr/bin/env Rscript
# ------- Working directory handling -------
# When run via Rscript, switch working directory to the script's directory
# so relative paths work consistently.
try({
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if(length(file_arg) == 1){
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    setwd(dirname(script_path))
  }
}, silent = TRUE)

pkgs <- c("readr","dplyr","stringr","tidyr","lme4","broom.mixed","purrr","ggplot2","blme","scales","rlang")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(lme4);  library(broom.mixed); library(purrr); library(ggplot2)
  library(scales); library(rlang)
})

# --------------------- config ---------------------
paths <- list(
  orfs_plc   = "results/tables/arg_orfs_with_plc.tsv",
  contigs    = "results/tables/contig_metrics.tsv",
  contig_cov = "results/tables/contig_coverage.tsv",              # optional
  plc_table  = "results/tables/plc_class_by_contig.tsv",
  pip_cov    = "results/tables/pipdb_read_coverage_annot.tsv",    # annotated
  pip_len    = "results/tables/pipdb_ref_lengths.tsv",
  pip_call   = "results/tables/pipdb_backbone_presence.tsv",      # optional
  spacer     = "results/tables/spacer_edges.tsv",
  meta       = "metadata/sample_meta.tsv"                         # optional
)

min_breadth <- 0.70
min_depth   <- 1.0
min_contig_len_for_glmm <- 1000
min_rows_per_class <- 30
fdr_method <- "BH"
alpha_fdr  <- 0.05
# --------------------------------------------------

dir.create("results/stats", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figs",  showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(paste0(..., collapse=""), "\n")

# Normalize sample IDs coming from filenames or BAM/RG-derived names to a canonical form
normalize_sample <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- stringr::str_replace(x, "#.*$", "")       # drop RG fragments after '#'
  x <- stringr::str_replace(x, "@.*$", "")       # drop RG fragments after '@'
  x <- stringr::str_replace(x, "\\.bam$", "")    # drop file extension
  x <- stringr::str_replace(x, "\\.cram$", "")
  x <- stringr::str_replace(x, "\\.sorted$", "")
  x <- stringr::str_replace(x, "_dedup(ed)?$", "")
  x <- stringr::str_replace(x, "_L\\d{3}", "")   # lane
  x <- stringr::str_replace(x, "-L\\d{3}", "")
  x <- stringr::str_replace(x, "_R[12](_001)?", "")
  x <- stringr::str_replace(x, "_\\d{3}$", "")
  x <- stringr::str_replace(x, "\\..*$", "")     # drop anything after first dot
  x <- stringr::str_trim(x)
  x
}

# Sanitize display names for plasmid backbones
clean_display_name <- function(x, fallback) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  bad <- is.na(x) | x == "" | x == "-" | tolower(x) %in% c("na","nan","none","null","0","0.0","0.000000")
  ifelse(bad, fallback, x)
}

need <- c(paths$orfs_plc, paths$contigs, paths$pip_cov, paths$pip_len, paths$spacer)
missing <- need[!file.exists(need)]
if(length(missing)) stop("Missing required files:\n", paste(missing, collapse="\n"))

# ---------- load core tables ----------
orfs <- read_tsv(paths$orfs_plc, show_col_types = FALSE) %>%
  mutate(
    orf_length_aa = suppressWarnings(readr::parse_number(as.character(orf_length_aa))),
    gc            = suppressWarnings(readr::parse_number(as.character(gc))),
    length        = suppressWarnings(readr::parse_number(as.character(length))),
    arg_class     = ifelse(is.na(arg_class) | arg_class=="", "unclassified", arg_class)
  )

# if plc_class is missing in orfs, try to recover from plc table
if (!has_name(orfs, "plc_class")) {
  if (!file.exists(paths$plc_table)) stop("orfs file lacks 'plc_class' and plc table not found.")
  plc_fix <- read_tsv(paths$plc_table, show_col_types = FALSE) %>%
    select(contig_id, plc_class)
  orfs <- orfs %>% left_join(plc_fix, by = "contig_id")
}
# derive binary flag
orfs <- orfs %>% mutate(plasmid_like = as.integer(plc_class != "chromosomal"))

contigs <- read_tsv(paths$contigs, show_col_types = FALSE) %>%
  mutate(length = as.integer(length), gc = as.numeric(gc))

cov_tbl <- if (file.exists(paths$contig_cov)) {
  read_tsv(paths$contig_cov, show_col_types = FALSE)
} else tibble(sample=character(), contig_id=character(), mean_depth=numeric(), breadth=numeric())

plc_features <- if (file.exists(paths$plc_table)) {
  read_tsv(paths$plc_table, show_col_types = FALSE) %>%
    mutate(has_relaxase=as.integer(has_relaxase),
           has_orit=as.integer(has_orit))
} else tibble(contig_id=character(), has_relaxase=integer(), has_orit=integer(),
              pipdb_hit=integer(), mobrecon_pos=integer(), plc_class=character())

# relaxase/oriT within plasmid-like subset
plc_plasmid <- plc_features %>% filter(plc_class != "chromosomal")
relax_sum <- plc_plasmid %>%
  summarise(n_contigs = n(),
            n_relaxase = sum(has_relaxase, na.rm=TRUE),
            n_orit     = sum(has_orit,     na.rm=TRUE)) %>%
  mutate(p_relaxase = n_relaxase / pmax(n_contigs,1),
         p_orit     = n_orit / pmax(n_contigs,1))

# ---------- join covariates ----------
dat <- orfs %>%
  left_join(contigs %>% select(contig_id, contig_len=length, contig_gc=gc), by="contig_id") %>%
  left_join(cov_tbl, by=c("sample","contig_id")) %>%
  mutate(
    gc = ifelse(!is.na(contig_gc), contig_gc, gc),
    mean_depth = ifelse(is.na(mean_depth), 0, mean_depth)
  )

# sample factors (fallback parse if no meta)
if (file.exists(paths$meta)) {
  meta <- read_tsv(paths$meta, show_col_types = FALSE) %>%
    select(sample, plant, campaign) %>%
    mutate(meta_norm = normalize_sample(sample))
} else {
  meta <- dat %>%
    distinct(sample) %>%
    mutate(
      plant    = stringr::str_match(sample, "P\\d+-(\\d+)_")[,2],
      plant    = ifelse(is.na(plant), "UNK", paste0("PL", plant)),
      campaign = stringr::str_match(sample, "_(S\\d+)")[,2],
      campaign = ifelse(is.na(campaign), "UNK", campaign),
      meta_norm = normalize_sample(sample)
    )
}
dat <- dat %>% left_join(meta, by="sample") %>% mutate(plant=factor(plant), campaign=factor(campaign))
if (min_contig_len_for_glmm > 0) dat <- dat %>% filter(!is.na(contig_len), contig_len >= min_contig_len_for_glmm)

# ---------- GLMM: class-specific enrichment ----------
dat <- dat %>% mutate(arg_class_norm = str_replace_all(arg_class, "\\s*;\\s*", ";"))

tokens <- dat$arg_class_norm %>%
  str_split(";", simplify = FALSE) %>% unlist() %>% str_trim() %>%
  { x <- .; x[x != "" & x != "unclassified"] } %>% unique()

has_token <- function(x, tok) {
  pat <- paste0("(^|;)", stringr::fixed(tok), "($|;)")
  stringr::str_detect(x, pat)
}

cls_tab <- tibble(arg_class = tokens) %>%
  rowwise() %>%
  mutate(n = sum(has_token(dat$arg_class_norm, arg_class)),
         pos = sum(dat$plasmid_like[has_token(dat$arg_class_norm, arg_class)]),
         neg = n - pos) %>%
  ungroup()

classes <- cls_tab %>%
  filter(n >= min_rows_per_class, pos > 0, neg > 0) %>%
  arrange(desc(n)) %>% pull(arg_class)

fit_one <- function(cls, use_depth = any(dat$mean_depth > 0, na.rm = TRUE)) {
  d <- dat %>% mutate(is_cls = as.integer(has_token(.data$arg_class_norm, cls)))
  form <- if (use_depth) {
    plasmid_like ~ is_cls + scale(orf_length_aa) + scale(gc) + scale(log1p(mean_depth)) + (1|plant) + (1|campaign)
  } else {
    plasmid_like ~ is_cls + scale(orf_length_aa) + scale(gc) + (1|plant) + (1|campaign)
  }
  tryCatch({
    m <- glmer(form, data = d, family = binomial(),
               control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE,
                                      optCtrl = list(maxfun = 2e5)))
  fx <- broom.mixed::tidy(m, effects = "fixed") %>% filter(.data$term == "is_cls")
    if(!nrow(fx)) return(NULL)
    est <- fx$estimate[1]; se <- fx$std.error[1]
    tibble(arg_class = cls,
           n_total = nrow(d), n_pos = sum(d$plasmid_like), n_neg = nrow(d)-sum(d$plasmid_like),
           logOR = est, SE = se, odds_ratio = exp(est),
           CI_low = exp(est - 1.96*se), CI_high = exp(est + 1.96*se),
           p = if("p.value" %in% names(fx)) fx$p.value[1] else 2*pnorm(-abs(est/se)))
  }, error = function(e) NULL)
}

glmm_res <- map_dfr(classes, fit_one) %>%
  mutate(p_adj = p.adjust(p, method=fdr_method)) %>%
  arrange(p_adj)
write_tsv(glmm_res, "results/stats/arg_class_plasmid_enrich_glmm.tsv")

# ---------- PIPdb presence ----------
pip_cov <- read_tsv(paths$pip_cov, show_col_types = FALSE) %>%
  filter(sample != "sample") %>%                                 # drop any repeated header row
  mutate(
    mean_depth       = suppressWarnings(readr::parse_number(as.character(mean_depth))),
    covered_fraction = suppressWarnings(readr::parse_number(as.character(covered_fraction)))
  )
if ("plasmid_name" %in% names(pip_cov)) {
  pip_cov <- pip_cov %>% mutate(plasmid_name = as.character(plasmid_name))
}

pip_presence <- if (file.exists(paths$pip_call)) {
  read_tsv(paths$pip_call, show_col_types = FALSE) %>%
    mutate(present = as.integer(present)) %>%
    mutate(sample_norm = normalize_sample(sample))
} else {
  pip_cov %>%
    mutate(present = as.integer(covered_fraction >= min_breadth & mean_depth >= min_depth)) %>%
    select(sample, plasmid_ref, present, mean_depth, covered_fraction, covered_bases, ref_len) %>%
    mutate(sample_norm = normalize_sample(sample))
}
write_tsv(pip_presence, "results/stats/pipdb_backbone_presence.tsv")

if (exists("meta")) {
  # Restrict to samples present in metadata
  meta_norm <- unique(meta$meta_norm)
  pip_presence <- pip_presence %>% filter(sample_norm %in% meta_norm)
}
total_samples <- n_distinct(pip_presence$sample_norm)

pip_names <- tryCatch({
  tmp <- pip_cov %>% select(plasmid_ref, plasmid_name) %>% distinct()
  if (!"plasmid_name" %in% names(tmp)) tmp$plasmid_name <- NA_character_
  # pick the most frequent non-empty name per plasmid_ref
  tmp %>%
    mutate(plasmid_name = as.character(plasmid_name)) %>%
    mutate(valid = !(is.na(plasmid_name) | plasmid_name == "" | plasmid_name == "-")) %>%
    group_by(plasmid_ref, plasmid_name) %>% summarise(n = n(), .groups = "drop") %>%
    arrange(plasmid_ref, desc(n)) %>%
    group_by(plasmid_ref) %>% slice_head(n = 1) %>% ungroup() %>%
    select(plasmid_ref, plasmid_name)
}, error = function(e) tibble(plasmid_ref=character(), plasmid_name=character()))

pip_summary <- pip_presence %>%
  filter(present == 1) %>%
  distinct(sample_norm, plasmid_ref) %>%
  left_join(pip_names, by = "plasmid_ref") %>%
  mutate(display_name = clean_display_name(plasmid_name, plasmid_ref)) %>%
  distinct(sample_norm, display_name) %>%
  count(display_name, name = "n_samples") %>%
  arrange(desc(n_samples))
any_present <- pip_presence %>%
  group_by(sample_norm) %>% summarise(any = any(present == 1, na.rm=TRUE), .groups="drop") %>%
  summarise(n = sum(any, na.rm=TRUE)) %>% pull(n)
top_refs_labels <- pip_summary %>% slice_head(n = 9) %>%
  mutate(lbl = sprintf("%s in %d/%d sites", display_name, n_samples, .env$total_samples)) %>% pull(lbl)

# ---------- CRISPR spacer edges ----------
sp_edges <- read_tsv(paths$spacer, show_col_types = FALSE)
sp_n_edges   <- nrow(sp_edges)
sp_n_spacers <- sum(sp_edges$n_unique_spacers)

top_hosts <- sp_edges %>%
  filter(!is.na(host_lineage), host_lineage != "NA") %>%
  group_by(host_lineage) %>% summarise(spacers=sum(n_unique_spacers), .groups="drop") %>%
  arrange(desc(spacers)) %>% slice_head(n=5)

top_plasmids <- sp_edges %>%
  group_by(plasmid_id) %>% summarise(spacers=sum(n_unique_spacers), .groups="drop") %>%
  arrange(desc(spacers)) %>% slice_head(n=5)

# ---------- Global counts ----------
total_orfs   <- nrow(orfs)
plasmid_orfs <- sum(orfs$plasmid_like, na.rm=TRUE)
pct_plasmid  <- 100 * plasmid_orfs / total_orfs

# protect against missing plc_class
if (!has_name(orfs, "plc_class")) orfs <- orfs %>% mutate(plc_class = "unknown")
plc_counts <- orfs %>%
  mutate(plc_class = if_else(is.na(plc_class) | plc_class=="", "unknown", plc_class)) %>%
  count(plc_class, name="n") %>% arrange(desc(n))

# If relax/orit summary failed (empty), fill with NA row
if(!nrow(relax_sum)) {
  relax_sum <- tibble(n_contigs=NA_integer_, n_relaxase=NA_integer_, n_orit=NA_integer_,
                      p_relaxase=NA_real_, p_orit=NA_real_)
}

n_sig <- sum(glmm_res$p_adj <= alpha_fdr, na.rm=TRUE)
top_sig <- glmm_res %>%
  filter(p_adj <= alpha_fdr) %>%
  arrange(desc(odds_ratio)) %>% slice_head(n=5)

# ---------- Figure (log scale needs >0) ----------
if(nrow(glmm_res)){
  plt_dat <- glmm_res %>% filter(is.finite(odds_ratio), is.finite(CI_low), is.finite(CI_high),
                                 odds_ratio > 0, CI_low > 0, CI_high > 0)
  if(nrow(plt_dat)){
    p_enrich <- ggplot(plt_dat, aes(x=odds_ratio, y=reorder(arg_class, odds_ratio))) +
      geom_point() +
      geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0) +
      geom_vline(xintercept=1, linetype="dashed") +
      scale_x_log10() +
      labs(x="Odds ratio (plasmid-like carriage)", y=NULL,
           title="ARG class enrichment on plasmid-like contigs",
           subtitle=paste0("FDR: ", fdr_method, " (α=", alpha_fdr, ")")) +
      theme_bw()
    ggsave("results/figs/arg_class_plasmid_enrich_forest.png", p_enrich, width=7, height=6, dpi=300)
    # HD poster version
    dir.create("HD_poster_figures", showWarnings = FALSE, recursive = TRUE)
    ggsave("HD_poster_figures/abstract5_ARG_class_enrichment_HD.png", p_enrich, width=7, height=6, dpi=900)
    ggsave("HD_poster_figures/abstract5_ARG_class_enrichment_HD.pdf", p_enrich, width=7, height=6)
  }
}

# ---------- Copy-paste abstract lines ----------
blk <- c(
  "\n=== ABSTRACT5: copy-paste lines ===\n",
  sprintf("We profiled %s ARG ORFs across Metro Manila sites; %s (%.1f%%) were on plasmid-like contigs (non-chromosomal class).",
          scales::comma(total_orfs), scales::comma(plasmid_orfs), pct_plasmid),
  paste0("Plasmid-context breakdown: ",
         paste0(plc_counts$plc_class, "=", plc_counts$n, collapse="; "), "."),
  if (!any(is.na(relax_sum$n_contigs)))
    sprintf("Across %s plasmid-like contigs with MOB-suite calls, we detected relaxases on %s (%.2f%%) and oriT on %s (%.2f%%).",
            scales::comma(relax_sum$n_contigs),
            scales::comma(relax_sum$n_relaxase), 100*relax_sum$p_relaxase,
            scales::comma(relax_sum$n_orit),     100*relax_sum$p_orit) else NULL,
  sprintf("%s ARG classes met inclusion for mixed models; %s were significantly enriched on plasmid-like contigs (FDR≤%.2f).",
          length(classes), n_sig, alpha_fdr),
  if (n_sig > 0) paste0(
    "Top enriched classes (OR [95% CI], FDR): ",
    paste0(
      top_sig %>%
        mutate(lbl = sprintf("%s %.2f [%.2f–%.2f], %.3g", arg_class, odds_ratio, CI_low, CI_high, p_adj)) %>%
        pull(lbl),
      collapse = " ; "),
    ".") else NULL,
  {
    n_samples <- total_samples
    any_present2 <- pip_presence %>%
      group_by(sample_norm) %>%
      summarise(any = any(present==1, na.rm=TRUE), .groups="drop") %>%
      summarise(n = sum(any, na.rm=TRUE)) %>% pull(n)
    c(
      sprintf("Read recruitment to PIPdb backbones detected presence in %s/%s samples (breadth≥%.0f%% & depth≥%.1fx).",
              any_present2, n_samples, 100*min_breadth, min_depth),
      if(length(top_refs_labels)) paste0("Most prevalent backbones: ",
                                         paste(top_refs_labels, collapse=" ; "), ".") else NULL
    )
  },
  sprintf("CRISPR evidence connected hosts and plasmids through %s edges supported by %s unique spacers.",
          scales::comma(sp_n_edges), scales::comma(sp_n_spacers)),
  if(nrow(top_hosts)) paste0("Top spacer-rich hosts: ",
                             paste(sprintf("%s (%s)", top_hosts$host_lineage, top_hosts$spacers), collapse=" ; "), ".") else NULL,
  if(nrow(top_plasmids)) paste0("Top spacer-hit plasmids: ",
                                paste(sprintf("%s (%s)", top_plasmids$plasmid_id, top_plasmids$spacers), collapse=" ; "),
                                ".") else NULL,
  "=== /ABSTRACT5 ===\n"
)

cat(paste(blk, collapse="\n"))
writeLines(blk, "results/stats/abstract5_text.txt")

if(nrow(glmm_res)) {
  glmm_res %>%
    mutate(across(c(odds_ratio, CI_low, CI_high, p, p_adj), ~signif(., 3))) %>%
    write_tsv("results/stats/abstract5_glmm_pretty.tsv")
}

say("Done. Copy-paste text saved at results/stats/abstract5_text.txt")