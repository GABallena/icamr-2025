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

# ------- Working directory handling -------
# Original had: setwd(~"<PROJECT_ROOT>") which is invalid R syntax (formula) and typo-prone.
# We'll attempt to set wd to the script's directory when run via Rscript.
try({
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if(length(file_arg) == 1){
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    setwd(dirname(script_path))
  }
}, silent = TRUE)
cat("Working directory:", getwd(),"\n")
# ====== 0) Setup ======
pkgs_core <- c(
  "tidyverse","readr","stringr","glue","janitor","fs",
  "zCompositions","compositions","vegan","brms","loo",
  "Matrix","patchwork","data.table","sf","readxl"
)
for(p in pkgs_core){ if(!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org") }
invisible(lapply(pkgs_core, library, character.only=TRUE))

# Resolve common dplyr verb masking (e.g., MASS::select)
if(requireNamespace("conflicted", quietly=TRUE)){
  conflicted::conflict_prefer("select","dplyr", quiet=TRUE)
  conflicted::conflict_prefer("filter","dplyr", quiet=TRUE)
  conflicted::conflict_prefer("lag","dplyr", quiet=TRUE)
  conflicted::conflict_prefer("var","stats", quiet=TRUE)
} else {
  # Fallback: explicitly bind dplyr verbs to local environment
  select   <- dplyr::select
  filter   <- dplyr::filter
  mutate   <- dplyr::mutate
  summarise<- dplyr::summarise
  transmute<- dplyr::transmute
  group_by <- dplyr::group_by
  ungroup  <- dplyr::ungroup
  arrange  <- dplyr::arrange
}

# brms backend handling: prefer cmdstanr, else rstan, else skip
if(requireNamespace("cmdstanr", quietly=TRUE)){
  brms_backend <- "cmdstanr"
} else if(requireNamespace("rstan", quietly=TRUE)) {
  brms_backend <- "rstan"
  message("cmdstanr missing; using rstan backend for brms.")
} else {
  brms_backend <- NA_character_
  message("WARNING: Neither cmdstanr nor rstan installed; Bayesian modeling will be skipped.")
}

# If you want to use cmdstan for brms (fast & reliable), uncomment this once:
# cmdstanr::install_cmdstan()  # one-time setup
options(mc.cores = parallel::detectCores())

# ====== 1) Paths & parameters ======
dir_shortbred <- Sys.getenv("SHORTBRED_DIR", unset = "shortbred")
dir_markers   <- "markers"
dir_ref       <- "ref"
dir_card      <- "ARG_db/card_db"

marker_faa    <- file.path(dir_markers, "CARD_ARG_markers.faa")
map_sb        <- file.path(dir_ref, "CARD_ARG_proteins.sb.map.tsv")   # sanitized->original ARO mapping
aro_idx       <- file.path(dir_card, "aro_index.tsv")
aro_cat       <- file.path(dir_card, "aro_categories.tsv")
neff_file     <- file.path(dir_shortbred, "N_eff.tsv")

# Your metadata file: one row per site/sample
# Columns example (rename in code below if your headers differ):
# sample, pop_density, poverty, literacy, healthcare_dist, sewerage, drainage, land_urban, land_ind, rain_wk, temp_wk, district, x, y
meta_file     <- "metadata/sites_covariates.tsv"  # optional; script will continue if absent
masterlist_file <- "P4 DATA MASTERLIST - P4 Y1 DATA.tsv"  # drives exclusion when present

# Core definition & LOD parameters
cpm_core_threshold <- 0.1   # X in abstract ("≥ X CPM")
core_prev_threshold <- 0.70 # ≥ 70% sites
lod_min_reads <- 3          # minimum reads to be "detected" at marker-level

# Modeling knobs
set.seed(42)
loo_kfold <- 5

# ====== 2) Helpers ======
read_shortbred <- function(path){
  # ShortBRED-Quantify common columns: Marker, Reads, Aligned_Bases, Marker_Length
  # We coerce column names gently; keep what exists.
  df <- suppressMessages(read_tsv(path, show_col_types = FALSE)) %>% janitor::clean_names()
  # Try to guess columns
  nm <- names(df)
  # New observed schema: family, count, hits, tot_marker_length
  marker_col <- nm[which(nm %in% c("marker","marker_id","id","family"))][1]
  reads_col  <- nm[which(nm %in% c("reads","count","n_reads","hits"))][1]
  len_col    <- nm[which(nm %in% c("marker_length","marker_length_bp","length","len","tot_marker_length"))][1]
  if(is.na(marker_col) || is.na(reads_col)){
    stop(glue("Could not find expected columns in {path}. Have: {paste(nm, collapse=', ')}"))
  }
  if(is.na(len_col)) len_col <- NA_character_
  df %>%
    transmute(marker_id = .data[[marker_col]],
              reads     = as.numeric(.data[[reads_col]]),
              marker_len_bp = if(!is.na(len_col)) as.numeric(.data[[len_col]]) else NA_real_)
}

read_fasta_lengths <- function(faa){
  # fallback if ShortBRED TSV lacks marker lengths
  lens <- c()
  cur  <- NA_character_
  con <- file(faa, "r")
  on.exit(close(con))
  while(TRUE){
    ln <- readLines(con, n=1)
    if(length(ln)==0) break
    if(startsWith(ln, ">")){
      cur <- sub("^>","",ln)
      lens[cur] <- 0
    } else {
      lens[cur] <- lens[cur] + nchar(ln)
    }
  }
  tibble(marker_id = names(lens), marker_len_bp = as.numeric(lens))
}

fmt_pct <- function(x, digits=1) paste0(sprintf(paste0("%.",digits,"f"), 100*x), "%")

# ====== 3) Load inputs ======
# Auto-detect alternative root (e.g., shortbred_transfer) if defaults missing ----
required_core <- c(dir_shortbred, dir_markers, dir_ref, dir_card,
                   marker_faa, map_sb, aro_idx, aro_cat)
if(any(!fs::file_exists(required_core) & !fs::dir_exists(required_core))){
  candidate_root <- "shortbred_transfer"
  if(fs::dir_exists(candidate_root)){
    alt_marker <- fs::path(candidate_root, "markers", "CARD_ARG_markers.faa")
    if(fs::file_exists(alt_marker)){
      message("Auto-detected alternative data root: ", candidate_root)
      dir_shortbred <- fs::path(candidate_root, "shortbred")
      dir_markers   <- fs::path(candidate_root, "markers")
      dir_ref       <- fs::path(candidate_root, "ref")
      dir_card      <- fs::path(candidate_root, "ARG_db/card_db")
      marker_faa    <- fs::path(dir_markers, "CARD_ARG_markers.faa")
      map_sb        <- fs::path(dir_ref, "CARD_ARG_proteins.sb.map.tsv")
      aro_idx       <- fs::path(dir_card, "aro_index.tsv")
      aro_cat       <- fs::path(dir_card, "aro_categories.tsv")
      neff_file     <- fs::path(dir_shortbred, "N_eff.tsv")
    }
  }
}

# Recompute required_core after auto-detect
required_core <- c(dir_shortbred, dir_markers, dir_ref, dir_card,
                   marker_faa, map_sb, aro_idx, aro_cat)
core_missing <- required_core[!fs::file_exists(required_core) & !fs::dir_exists(required_core)]
if(length(core_missing) > 0){
  cat("Missing required core inputs (still not found):\n")
  writeLines(paste0("  - ", core_missing))
  stop("Aborting: required core ARG data missing.")
}

meta_available <- fs::file_exists(meta_file)
if(!meta_available){
  message("Metadata file not found (", meta_file, "). Downstream db-RDA & Bayesian sections will be skipped.")
}
stopifnot(dir_exists(dir_shortbred))
shortbred_files <- dir_ls(dir_shortbred, glob = "*.shortbred.tsv")
stopifnot(length(shortbred_files) > 0)

message(glue("Found {length(shortbred_files)} ShortBRED result files."))

# Build marker length reference
marker_len_tbl <- read_fasta_lengths(marker_faa)

# Merge per-sample counts
counts_long <- map_df(shortbred_files, function(f){
  s <- basename(f) %>% str_remove("\\.shortbred\\.tsv$")
  read_shortbred(f) %>% mutate(sample = s)
})

# Attach marker lengths if missing
counts_long <- counts_long %>%
  select(marker_id, sample, reads, marker_len_bp) %>%
  group_by(marker_id) %>%
  mutate(marker_len_bp = coalesce(marker_len_bp, dplyr::first(na.omit(marker_len_bp)))) %>%
  ungroup() %>%
  left_join(marker_len_tbl, by="marker_id", suffix=c("",".fa")) %>%
  mutate(marker_len_bp = coalesce(marker_len_bp, marker_len_bp.fa)) %>%
  { if("marker_len_bp.fa" %in% names(.)) select(., -marker_len_bp.fa) else . }

# N_eff (for CPM denominator)
# Optional: if not found, create a placeholder (all NA) so downstream joins work.
neff <- if (fs::file_exists(neff_file)) {
  read_tsv(neff_file, col_names = c("sample","N_eff"), show_col_types = FALSE) %>% distinct()
} else {
  message("N_eff.tsv not found at ", neff_file, "; proceeding without it (CPM will use raw ShortBRED counts).")
  tibble(sample = unique(counts_long$sample), N_eff = NA_real_)
}

# Metadata (optional)
if(meta_available){
  meta <- read_tsv(meta_file, show_col_types = FALSE) %>% janitor::clean_names()
  stopifnot("sample" %in% names(meta))
} else {
  meta <- tibble(sample = unique(counts_long$sample))
}

# ====== Sample exclusion via masterlist ======
if(fs::file_exists(masterlist_file)){
  master_raw <- suppressMessages(read_tsv(masterlist_file, show_col_types = FALSE))
  # Clean column names but preserve originals for mapping sample IDs
  master <- master_raw
  orig_names <- names(master)
  names(master) <- janitor::make_clean_names(orig_names)
  # Identify likely columns
  sample_code_col <- which(str_detect(names(master), '^sample_code$|^sample_id$|^sample$'))[1]
  sample_type_col <- which(str_detect(names(master), 'sample_type'))[1]
  sample_desc_col <- which(str_detect(names(master), 'sample_description'))[1]
  if(!is.na(sample_code_col)){
    master <- master %>% rename(sample = all_of(names(master)[sample_code_col]))
  } else {
    # Fallback: try first column if it matches pattern like SAMPLE-XX
    if('x' %in% names(master) && all(str_detect(master$x[!is.na(master$x)], '^SAMPLE-'))) master <- master %>% rename(sample = x)
  }
  if(!'sample' %in% names(master)){
    warning("Masterlist present but no sample code column detected; reverting to regex exclusion.")
  } else {
    types <- if(!is.na(sample_type_col)) master[[sample_type_col]] else NA_character_
    descs <- if(!is.na(sample_desc_col)) master[[sample_desc_col]] else NA_character_
    df_excl <- tibble(sample = master$sample,
                      type = types,
                      desc = descs) %>%
      mutate(across(c(type, desc), as.character))
    excl <- df_excl %>%
      filter(
        str_detect(coalesce(type,''), regex('control', ignore_case=TRUE)) |
        str_detect(coalesce(type,''), regex('optimization|opt\\.', ignore_case=TRUE)) |
        str_detect(coalesce(desc,''), regex('control', ignore_case=TRUE)) |
        str_detect(coalesce(desc,''), regex('optimization|opt\\.', ignore_case=TRUE)) |
        str_detect(sample, regex('^SAMPLE-01$', ignore_case=TRUE)) |
        str_detect(sample, regex('^SAMPLE-23$', ignore_case=TRUE))
      ) %>% pull(sample) %>% unique()
    if(length(excl) > 0){
      # Normalize ShortBRED sample IDs to base form (strip trailing _S\d+)
      counts_long <- counts_long %>% mutate(sample_base = sub("_S[0-9]+$","", sample))
      meta        <- meta %>% mutate(sample_base = sub("_S[0-9]+$","", sample))
      # Determine which actual sample IDs match exclusion list via base
      excl_bases <- excl
      to_drop <- counts_long %>% filter(sample_base %in% excl_bases) %>% pull(sample) %>% unique()
      if(length(to_drop) > 0){
        message("Excluding via masterlist annotations (matched with _S suffix handling): ", paste(to_drop, collapse=", "))
        counts_long <- counts_long %>% filter(!sample %in% to_drop) %>% select(-sample_base)
        meta <- meta %>% filter(!sample %in% to_drop) %>% select(-sample_base)
      } else {
        message("Exclusion codes found but no sample IDs matched after suffix normalization.")
        counts_long <- counts_long %>% select(-sample_base)
        meta <- meta %>% select(-sample_base)
      }
      excluded_samples_summary <- if(exists("to_drop") && length(to_drop)>0) paste(to_drop, collapse=", ") else "None"
    } else {
      excluded_samples_summary <- "None"
    }
  }
}

# Derive covariates from masterlist (physicochemical) for remaining samples -----------------
covariate_summary <- NULL
if(fs::file_exists(masterlist_file)){
  suppressMessages({ master_cov_raw <- read_tsv(masterlist_file, show_col_types = FALSE) })
  m2 <- master_cov_raw
  names(m2) <- janitor::make_clean_names(names(m2))
  sc_col2 <- which(str_detect(names(m2), '^sample_code$|^sample$'))[1]
  if(!is.na(sc_col2)) names(m2)[sc_col2] <- 'sample_base'
  # Numeric detection (columns convertible to numeric with some variation)
  numeric_cols <- names(m2)[sapply(m2, function(x){ suppressWarnings({ all(!nzchar(setdiff(na.omit(as.character(x)), as.character(as.numeric(as.character(x)))))) }) })]
  m2[numeric_cols] <- lapply(m2[numeric_cols], function(x) suppressWarnings(as.numeric(x)))
  num_meta <- m2 %>% select(any_of(c('sample_base', numeric_cols)))
  if('sample_base' %in% names(num_meta)){
    var_cols <- setdiff(names(num_meta), 'sample_base')
    keep_cols <- var_cols[sapply(num_meta[var_cols], function(v) sum(!is.na(v)) >= 5 && dplyr::n_distinct(na.omit(v)) >= 4)]
    drop_patterns <- c('total_dna','qubit','agarose','tapestation','library','band','smear')
    keep_cols <- keep_cols[!str_detect(keep_cols, paste(drop_patterns, collapse='|'))]
    meta_cov <- num_meta %>% select(sample_base, all_of(keep_cols)) %>% distinct()
    # Attach
    if('sample_base' %in% names(meta)){
      meta <- meta %>% left_join(meta_cov, by='sample_base')
    }
    if(length(keep_cols)>0){ covariate_summary <- glue('{length(keep_cols)} covariates loaded: {paste(keep_cols, collapse=", ")}') }
  }
}

# ====== Spatial Administrative Boundary Mapping (1 km buffers + diagnostics) ======
admin_mapping <- NULL
admin_summary_line <- NULL
try({
  if(fs::file_exists(masterlist_file)){
    ml_coords <- suppressMessages(read_tsv(masterlist_file, show_col_types = FALSE))
    mlc <- ml_coords; names(mlc) <- janitor::make_clean_names(names(mlc))
    sc_col3 <- which(str_detect(names(mlc), '^sample_code$|^sample$'))[1]
    if(!is.na(sc_col3)) names(mlc)[sc_col3] <- 'sample_base'
    lat_candidates <- names(mlc)[str_detect(names(mlc), 'lat')] ; lon_candidates <- names(mlc)[str_detect(names(mlc), 'lon|long')]
    lat_col <- if(length(lat_candidates)>0) lat_candidates[which.max(nchar(lat_candidates))] else NA_character_
    lon_col <- if(length(lon_candidates)>0) lon_candidates[which.max(nchar(lon_candidates))] else NA_character_
    if(!is.na(lat_col) && !is.na(lon_col) && 'sample_base' %in% names(mlc)){
      coord_df <- mlc %>% select(sample_base, !!lon_col, !!lat_col) %>% rename(longitude=!!lon_col, latitude=!!lat_col) %>% mutate(across(c(longitude,latitude), as.numeric)) %>% filter(!is.na(longitude), !is.na(latitude)) %>% distinct()
      kept <- if('sample_base' %in% names(counts_long)) counts_long %>% distinct(sample, sample_base) else counts_long %>% distinct(sample) %>% mutate(sample_base=sub('_S[0-9]+$','',sample))
      pts <- kept %>% left_join(coord_df, by='sample_base') %>% filter(!is.na(longitude), !is.na(latitude))
      message('[SPATIAL] Coordinate detection: lat_col=', lat_col, ' lon_col=', lon_col, ' points_with_coords=', nrow(pts))
      if(nrow(pts) > 0){
        sf_pts <- sf::st_as_sf(pts, coords=c('longitude','latitude'), crs=4326, remove=FALSE)
        # Derive padded bbox for optional pre-filtering of polygons
        bb <- sf::st_bbox(sf_pts)
        pad <- 0.15
        crop_bbox <- c(xmin=bb$xmin-pad, ymin=bb$ymin-pad, xmax=bb$xmax+pad, ymax=bb$ymax+pad)
        message('[SPATIAL] Sample bbox: ', paste(signif(unlist(bb),6), collapse=', '))
        message('[SPATIAL] Crop bbox (pad 0.15°): ', paste(signif(crop_bbox,6), collapse=', '))
        utm51 <- 32651
        sf_pts_utm <- try(sf::st_transform(sf_pts, utm51), silent=TRUE)
        if(inherits(sf_pts_utm,'try-error')) sf_pts_utm <- sf::st_transform(sf_pts, 3857)
        sf_buf <- sf::st_buffer(sf_pts_utm, dist=1000) %>% sf::st_transform(4326)
        read_layer <- function(lv){
          f <- file.path('phl_adm_psa_namria_20231106_shp', sprintf('phl_admbnda_adm%d_psa_namria_20231106.shp', lv))
          if(!fs::file_exists(f)) return(NULL)
          suppressMessages(suppressWarnings(sf::st_read(f, quiet=TRUE)))
        }
        admin_levels <- list()
        join_stats <- list()
        for(lv in 1:4){
          message('[SPATIAL] --- Processing ADM', lv, ' ---')
          tryCatch({
            lyr_path <- file.path('phl_adm_psa_namria_20231106_shp', sprintf('phl_admbnda_adm%d_psa_namria_20231106.shp', lv))
            message('[SPATIAL] ADM', lv, ' path: ', lyr_path, ' exists=', fs::file_exists(lyr_path))
            lyr <- read_layer(lv)
            if(is.null(lyr)) { message('[SPATIAL] ADM', lv, ' layer missing (NULL).'); next }
            message('[SPATIAL] ADM', lv, ' n_polygons pre-crop=', nrow(lyr))
            suppressWarnings({ lyr <- sf::st_make_valid(lyr) })
            bbx <- try(sf::st_as_sfc(sf::st_bbox(c(xmin=crop_bbox['xmin'], ymin=crop_bbox['ymin'], xmax=crop_bbox['xmax'], ymax=crop_bbox['ymax']), crs = sf::st_crs(lyr))), silent=TRUE)
            lyr_crop <- try(suppressWarnings(sf::st_intersection(lyr, bbx)), silent=TRUE)
            if(!inherits(lyr_crop,'try-error') && !is.null(lyr_crop) && nrow(lyr_crop) > 0){
              message('[SPATIAL] ADM', lv, ' cropped to ', nrow(lyr_crop), ' polygons within bbox.')
              lyr <- lyr_crop
            } else {
              message('[SPATIAL] ADM', lv, ' crop skipped (try-error or empty).')
            }
            # Intersection-based and area-weighted selection per sample
            buf <- sf_buf %>% dplyr::select(sample)
            jo_int <- NULL
            s2_prev <- try(sf::sf_use_s2(), silent=TRUE)
            jo_try <- try(suppressMessages(suppressWarnings(sf::st_intersection(buf, lyr))), silent=TRUE)
            if(!inherits(jo_try, 'try-error')){
              jo_int <- jo_try
            } else {
              message('[SPATIAL] ADM', lv, ' st_intersection failed; disabling s2 and retrying...')
              try(sf::sf_use_s2(FALSE), silent=TRUE)
              jo_try2 <- try(suppressMessages(suppressWarnings(sf::st_intersection(buf, lyr))), silent=TRUE)
              if(!inherits(jo_try2, 'try-error')) jo_int <- jo_try2
              if(!inherits(s2_prev, 'try-error')) try(sf::sf_use_s2(s2_prev), silent=TRUE)
            }
            if(!is.null(jo_int) && nrow(jo_int) > 0){
              jo_int$._area <- as.numeric(sf::st_area(jo_int))
              name_pref <- sprintf('ADM%d_EN', lv)
              name_col <- if(name_pref %in% names(jo_int)) name_pref else grep(sprintf('ADM%d', lv), names(jo_int), value=TRUE)[1]
              jo_df <- jo_int %>%
                sf::st_drop_geometry() %>%
                dplyr::group_by(sample) %>%
                dplyr::slice_max(order_by = ._area, n = 1, with_ties = FALSE) %>%
                dplyr::summarise(!!glue('adm{lv}_unit') := dplyr::first(.data[[name_col]]), .groups='drop')
              admin_levels[[lv]] <- jo_df
              join_stats[[as.character(lv)]] <- list(level=lv, name_col=name_col, matched=sum(!is.na(jo_df[[glue('adm{lv}_unit')]])))
              message('[SPATIAL] ADM', lv, ' area-based assignment success: name_col=', name_col, ' matched_samples=', join_stats[[as.character(lv)]]$matched)
            } else {
              message('[SPATIAL] ADM', lv, ' intersection produced no overlaps.')
            }
            rm(lyr); gc()
          }, error=function(e){
            message('[SPATIAL] ERROR in ADM', lv, ': ', conditionMessage(e))
          })
        }
        dir_create('results')
        if(length(admin_levels)>0){
          admin_mapping <- reduce(admin_levels, full_join, by='sample') %>% left_join(pts %>% select(sample, longitude, latitude), by='sample')
          suppressWarnings(write_tsv(admin_mapping, 'results/sample_admin_mapping.tsv'))
          message('[SPATIAL] Wrote admin mapping file with ', nrow(admin_mapping), ' rows.')
          adm2_n <- if('adm2_unit' %in% names(admin_mapping)) n_distinct(admin_mapping$adm2_unit, na.rm=TRUE) else NA
          adm3_n <- if('adm3_unit' %in% names(admin_mapping)) n_distinct(admin_mapping$adm3_unit, na.rm=TRUE) else NA
          admin_summary_line <- glue('Admin coverage: ADM2={adm2_n} distinct; ADM3={adm3_n} distinct')
        } else {
          # Write empty placeholder
          admin_mapping <- pts %>% select(sample, longitude, latitude)
          suppressWarnings(write_tsv(admin_mapping, 'results/sample_admin_mapping.tsv'))
          message('[SPATIAL] Wrote placeholder admin mapping file (no overlaps). Rows=', nrow(admin_mapping))
          admin_summary_line <- 'Admin coverage: none'
          message('[SPATIAL] No administrative overlaps detected.')
        }
      } else {
        message('[SPATIAL] No usable point coordinates after merging masterlist.')
        admin_summary_line <- 'Admin coverage: no points'
      }
    } else {
      message('[SPATIAL] Could not detect latitude/longitude columns in masterlist.')
      admin_summary_line <- 'Admin coverage: coordinates missing'
    }
  }
}, silent=TRUE)
if(n_distinct(counts_long$sample) == 0){
  stop("All samples were excluded (masterlist criteria); adjust exclusion logic.")
}

# ====== 3b) Join external metadata sources (buildings, POIs, census) via authoritative ADM4 PCODEs ======
metadata_enrichment_summary <- NULL
try({
  if(exists('admin_mapping') && !is.null(admin_mapping) && 'adm4_unit' %in% names(admin_mapping)){
    # Build authoritative lookup from ADM4 shapefile (only once per run)
    adm4_lookup <- try({
      sf::st_read('phl_adm_psa_namria_20231106_shp/phl_admbnda_adm4_psa_namria_20231106.shp', quiet=TRUE) %>%
        sf::st_drop_geometry() %>%
        transmute(adm4_unit = ADM4_EN, adm4_pcode = ADM4_PCODE,
                  adm3_unit = ADM3_EN, adm3_pcode = ADM3_PCODE,
                  adm2_unit = ADM2_EN, adm2_pcode = ADM2_PCODE,
                  adm1_unit = ADM1_EN, adm1_pcode = ADM1_PCODE)
    }, silent=TRUE)
    if(inherits(adm4_lookup,'try-error')) adm4_lookup <- NULL
    # Direct name join (exact match) to get PCODEs
    admin_map_codes <- if(!is.null(adm4_lookup)) {
      joined <- admin_mapping %>% left_join(adm4_lookup %>% select(adm4_unit, adm4_pcode, adm3_pcode, adm2_pcode, adm1_pcode), by='adm4_unit')
      # If REGION codes did not map (many NA), attempt reconstruction: extract barangay numeric portion and substitute into pattern if adm3/adm2 codes available elsewhere.
      if(mean(is.na(joined$adm4_pcode)) > 0.8){
        # Attempt pattern-based generation for REGION: if adm4_unit like 'Barangay XYZ' or has trailing digits
        joined <- joined %>% mutate(
          brgy_num = stringr::str_extract(adm4_unit, '\\b[0-9]{1,3}\\b'),
          brgy_num_padded = ifelse(!is.na(brgy_num), stringr::str_pad(brgy_num, 4, pad='0'), NA_character_)
        )
        # Hypothetical REGION root codes (observed PH137401000 pattern in external datasets). We cannot guess full without authoritative table; leave placeholder if unresolvable.
        # If an adm3_pcode exists in lookup (none currently for PH137 in shapefile), could concatenate.
        # For now, retain constructed placeholder code to aid manual reconciliation.
        joined <- joined %>% mutate(adm4_pcode = ifelse(is.na(adm4_pcode) & !is.na(brgy_num_padded), paste0('PH1374010', brgy_num_padded), adm4_pcode))
      }
      joined
    } else admin_mapping
    # External datasets
    path_meta <- 'philippine_metadata'
    read_quiet <- function(path){ suppressMessages(suppressWarnings(readr::read_csv(path, show_col_types = FALSE))) }
    sx <- function(df){ if(is.null(df)) return(NULL); names(df) <- janitor::make_clean_names(names(df)); df }
    f_google <- file.path(path_meta,'google_open_buildings.csv')
    f_tm     <- file.path(path_meta,'tm_open_buildings.csv')
    f_water  <- file.path(path_meta,'osm_poi_water_body.csv')
    f_san    <- file.path(path_meta,'osm_poi_sanitation.csv')
    f_census <- file.path(path_meta,'2020-census-total-popn-brgy_adm4_new-pcode.xlsx')
    safe_csv <- function(f) if(fs::file_exists(f)) sx(read_quiet(f)) else NULL
    google_tbl <- safe_csv(f_google)
    tm_tbl     <- safe_csv(f_tm)
    water_tbl  <- safe_csv(f_water)
    san_tbl    <- safe_csv(f_san)
    census_tbl <- if(fs::file_exists(f_census) && requireNamespace('readxl', quietly=TRUE)){
      ct <- suppressMessages(readxl::read_excel(f_census, sheet=1)); names(ct) <- janitor::make_clean_names(names(ct)); ct
    } else NULL
    # Normalize census -> adm4_pcode
    if(!is.null(census_tbl)){
      pcode_col <- intersect(names(census_tbl), c('new_10_pcode','new_10_digit_psgc','bgy_code'))[1]
      if(!is.na(pcode_col)){
        census_tbl <- census_tbl %>% rename(adm4_pcode = !!pcode_col)
        # Standardize barangay name casing to match admin_map_codes
        census_tbl <- census_tbl %>% mutate(barangay_std = tolower(stringr::str_replace_all(barangay, '^barangay\\s+',''))) %>%
          select(adm4_pcode, census__barangay = barangay_std, census__pop2020 = any_of('x2020_census_popn'))
      } else census_tbl <- NULL
    }
    # Trim each dataset down to one row per adm4_pcode (while keeping primary metrics)
    dedupe <- function(df){ if(is.null(df) || !'adm4_pcode' %in% names(df)) return(NULL); df %>% filter(!is.na(adm4_pcode)) %>% group_by(adm4_pcode) %>% slice(1) %>% ungroup() }
    google_core <- dedupe(google_tbl)
    tm_core     <- dedupe(tm_tbl)
    water_core  <- dedupe(water_tbl)
    san_core    <- dedupe(san_tbl)
    # Column prefixing & selection (exclude geometry columns or large WKT)
    drop_geom_cols <- function(df){ if(is.null(df)) return(NULL); keep <- setdiff(names(df), c('geometry')); df[keep] }
    google_core <- drop_geom_cols(google_core)
    tm_core     <- drop_geom_cols(tm_core)
    # Sanitize column sets to avoid massive join if google not overlapping region
    join_pref <- function(left, right, prefix){
      if(is.null(right) || !'adm4_pcode' %in% names(right)) return(left)
      feat <- setdiff(names(right), 'adm4_pcode')
      right <- right %>% rename_with(~paste0(prefix,'__', .x), all_of(feat))
      left %>% left_join(right, by='adm4_pcode')
    }
    enriched <- admin_map_codes
    # If we lack valid adm4_pcode assignments (e.g., >80% NA), skip heavy joins to avoid all-NA columns.
    na_rate_codes <- mean(is.na(enriched$adm4_pcode))
    if(!is.na(na_rate_codes) && na_rate_codes < 0.8){
      enriched <- join_pref(enriched, google_core, 'google')
      enriched <- join_pref(enriched, tm_core, 'tm')
      enriched <- join_pref(enriched, water_core, 'water')
      enriched <- join_pref(enriched, san_core, 'san')
      if(!is.null(census_tbl)) enriched <- join_pref(enriched, census_tbl, 'census')
    } else {
      metadata_enrichment_summary <- glue('Metadata enrichment skipped: adm4_pcode coverage insufficient (NA rate={round(na_rate_codes,2)})')
    }
    # NA diagnostics
  meta_cols <- grep('^(google__|tm__|water__|san__|census__)', names(enriched), value=TRUE)
    na_fracs <- if(length(meta_cols)>0) colMeans(sapply(enriched[meta_cols], function(v) is.na(v) | v=='')) else numeric()
    mean_na <- if(length(na_fracs)>0) round(mean(na_fracs),3) else NA_real_
    dir_create('results')
    readr::write_tsv(enriched, 'results/enriched_sample_admin_metadata.tsv')
  if(is.null(metadata_enrichment_summary)) metadata_enrichment_summary <- glue('Metadata enrichment: authoritative join (adm4_pcode); datasets added cols={length(meta_cols)}; mean NA frac={mean_na}')
  }
}, silent=TRUE)

# ====== 4) Marker→ARO→Class mapping ======
# Robust map parsing: prefer numeric ARO column; fallback to old_header text.
map_sb_df_raw <- suppressMessages(read_tsv(map_sb, show_col_types = FALSE))
map_sb_df <- map_sb_df_raw %>% janitor::clean_names()

# Identify columns
marker_col <- intersect(names(map_sb_df), c("new_id","marker_id","id"))[1]
if(is.na(marker_col)) marker_col <- names(map_sb_df)[1]
aro_numeric_col <- intersect(names(map_sb_df), c("aro","aro_accession","aro_id"))[1]
old_header_col  <- intersect(names(map_sb_df), c("old_header","header","orig_header"))[1]

tmp <- map_sb_df %>% select(all_of(marker_col), any_of(c(aro_numeric_col, old_header_col)))
names(tmp)[1] <- "marker_id"

tmp <- tmp %>% mutate(
  aro = case_when(
    !is.na(aro_numeric_col) & aro_numeric_col %in% names(.) ~ {
      v_chr <- as.character(.data[[aro_numeric_col]])
      ifelse(str_detect(v_chr, "^ARO:"), v_chr, paste0("ARO:", v_chr))
    },
    !is.na(old_header_col) & old_header_col %in% names(.) ~ str_extract(.data[[old_header_col]], "ARO:\\d+"),
    TRUE ~ NA_character_
  )
) %>% select(marker_id, aro) %>% distinct()

missing_aro_frac <- mean(is.na(tmp$aro))
if(isTRUE(missing_aro_frac == 1)){
  warning("All ARO accessions failed to parse from map file. Downstream annotation will collapse to a single class.")
}

map_sb_df <- tmp

aro_index_raw <- read_tsv(aro_idx, show_col_types = FALSE)
aro_cats_raw  <- read_tsv(aro_cat, show_col_types = FALSE)

normalize_cols <- function(df){
  cn <- names(df)
  cn_norm <- cn %>%
    str_replace_all("[ /]+","_") %>%
    str_replace_all("[^A-Za-z0-9_]+","_") %>%
    str_to_lower()
  names(df) <- cn_norm
  df
}
aro_index <- aro_index_raw %>% normalize_cols()
aro_cats  <- aro_cats_raw %>% normalize_cols()

aro_key <- (intersect(names(aro_index), c("aro_accession","aro_acc","aro","aro_id")))[1]
cat_key <- (intersect(names(aro_cats),  c("aro_accession","aro_acc","aro","aro_id")))[1]
if(is.na(aro_key) || is.na(cat_key)){
  stop("Could not determine ARO key columns in index/categories tables.")
}

# Ensure ARO key columns are prefixed with "ARO:" for consistent joins
if(aro_key %in% names(aro_index)){
  aro_index[[aro_key]] <- ifelse(stringr::str_detect(as.character(aro_index[[aro_key]]), "^ARO:"),
                                 as.character(aro_index[[aro_key]]),
                                 paste0("ARO:", as.character(aro_index[[aro_key]])))
}
if(cat_key %in% names(aro_cats)){
  aro_cats[[cat_key]] <- ifelse(stringr::str_detect(as.character(aro_cats[[cat_key]]), "^ARO:"),
                                as.character(aro_cats[[cat_key]]),
                                paste0("ARO:", as.character(aro_cats[[cat_key]])))
}

class_col <- (intersect(names(aro_index), c("drug_class","drug_class_name","drug_classification","antibiotic_class")))[1]
if(is.na(class_col)) class_col <- (intersect(names(aro_index), c("resistance_mechanism","mechanism","amr_gene_family")))[1]
if(is.na(class_col)) class_col <- (intersect(names(aro_cats),  c("drug_class","drug_class_name","drug_classification","antibiotic_class")))[1]
if(is.na(class_col)) class_col <- (intersect(names(aro_cats),  c("resistance_mechanism","mechanism","amr_gene_family")))[1]
if(is.na(class_col)) class_col <- names(aro_index)[1]  # fallback

aro_annot <- map_sb_df %>%
  left_join(aro_index, by = setNames(aro_key, "aro")) %>%
  left_join(aro_cats , by = setNames(cat_key, "aro")) %>%
  rename(arg_class = !!class_col) %>%
  select(marker_id, aro, arg_class) %>% distinct()

# Diagnostics: number of distinct classes before prevalence filtering
diag_n_classes <- n_distinct(aro_annot$arg_class, na.rm = TRUE)
if(diag_n_classes <= 1){
  message("[DIAG] Only ", diag_n_classes, " ARG class detected after mapping. Check map/ARO joins.")
}

# ====== 5) Normalize to CPM / TPM and aggregate to class ======
# User directive: ignore N_eff.tsv (values are NA) and treat raw ShortBRED Count as already a comparable abundance.
fallback_neff_used <- TRUE
counts_long <- counts_long %>%
  left_join(neff, by="sample") %>%
  mutate(CPM = reads)  # direct use of raw counts as CPM surrogate

# RPK/TPM (optional but more length-aware)
counts_long <- counts_long %>%
  mutate(marker_len_kb = ifelse(!is.na(marker_len_bp), marker_len_bp/1000, NA_real_),
         RPK = ifelse(!is.na(marker_len_kb) & marker_len_kb>0, reads/marker_len_kb, NA_real_))

tpm_by_sample <- counts_long %>%
  group_by(sample) %>%
  mutate(TPM = if(all(!is.na(RPK)) && sum(RPK, na.rm=TRUE)>0) 1e6 * RPK / sum(RPK, na.rm=TRUE) else NA_real_) %>%
  ungroup()

# Aggregate to ARG class
class_cpm <- tpm_by_sample %>%
  left_join(aro_annot, by="marker_id") %>%
  mutate(arg_class = coalesce(arg_class, "Unannotated")) %>%
  group_by(sample, arg_class) %>%
  summarise(CPM = sum(CPM, na.rm=TRUE), TPM = sum(TPM, na.rm=TRUE), .groups="drop")

# Save class abundance (long and wide) for downstream analyses
try({
  dir_create("results")
  readr::write_tsv(class_cpm, "results/class_abundance_cpm.tsv")
  class_cpm_wide <- class_cpm %>%
    select(sample, arg_class, CPM) %>%
    tidyr::pivot_wider(names_from = arg_class, values_from = CPM, values_fill = 0) %>%
    arrange(sample)
  readr::write_tsv(class_cpm_wide, "results/class_abundance_cpm_wide.tsv")
  # total ARG load per sample (CPM and TPM)
  class_cpm %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(arg_total_cpm = sum(CPM, na.rm=TRUE),
                     arg_total_tpm = sum(TPM, na.rm=TRUE), .groups='drop') %>%
    readr::write_tsv("results/arg_total_by_sample.tsv")
}, silent = TRUE)

# ---- Diagnostics: ensure no silent zeroing ----
try({
  message('[DIAG] Per-sample total reads (marker level):')
  print(counts_long %>% group_by(sample) %>% summarise(total_reads = sum(reads, na.rm=TRUE)) %>% arrange(total_reads), n=Inf)
  message('[DIAG] Markers seen per sample:')
  print(counts_long %>% group_by(sample) %>% summarise(n_markers = sum(reads > 0, na.rm=TRUE)) %>% arrange(n_markers), n=Inf)
  message('[DIAG] Positive classes per sample:')
  print(class_cpm %>% group_by(sample) %>% summarise(n_classes_pos = sum(CPM > 0, na.rm=TRUE)) %>% arrange(n_classes_pos), n=Inf)
  message('[DIAG] Top classes overall:')
  print(class_cpm %>% group_by(arg_class) %>% summarise(tot = sum(CPM, na.rm=TRUE)) %>% arrange(desc(tot)) %>% head(20), n=20)
}, silent=TRUE)

# ====== 6) Class prevalence & core definition ======
n_samples <- n_distinct(class_cpm$sample)
prev_over_X <- class_cpm %>%
  group_by(arg_class) %>%
  summarise(
    n_sites_over_X = sum(CPM >= cpm_core_threshold, na.rm=TRUE),
    prevalence_over_X = n_sites_over_X / n_samples,
    mean_CPM = mean(CPM, na.rm=TRUE),
    median_CPM = median(CPM, na.rm=TRUE),
    mean_TPM = mean(TPM, na.rm=TRUE),
    .groups="drop"
  ) %>% arrange(desc(prevalence_over_X))

# TPM-based prevalence (sensitivity check)
prev_over_X_tpm <- class_cpm %>%
  group_by(arg_class) %>%
  summarise(n_sites_over_X = sum(TPM >= cpm_core_threshold, na.rm=TRUE),
            prevalence_over_X = n_sites_over_X / n_samples, .groups='drop')

core <- prev_over_X %>%
  filter(prevalence_over_X >= core_prev_threshold) %>%
  arrange(desc(prevalence_over_X))

if(nrow(core) == 0){
  core <- tibble(arg_class=character(), n_sites_over_X=integer(), prevalence_over_X=double(), mean_CPM=double(), median_CPM=double(), mean_TPM=double())
}

# ====== 7) LOD-based prevalence (reads >= lod_min_reads at any marker for a class) ======
prev_LOD <- counts_long %>%
  left_join(aro_annot, by="marker_id") %>%
  mutate(arg_class = coalesce(arg_class, "Unannotated"), detected = reads >= lod_min_reads) %>%
  group_by(sample, arg_class) %>%
  summarise(detected = any(detected), .groups="drop") %>%
  group_by(arg_class) %>%
  summarise(n_sites_detected = sum(detected), prevalence_LOD = n_sites_detected / n_samples, .groups="drop") %>%
  arrange(desc(prevalence_LOD))

# Presence-based core using LOD prevalence
core_LOD <- prev_LOD %>% filter(prevalence_LOD >= core_prev_threshold)


r2_adj <- NA_real_
X <- NULL; X_repl <- NULL; X_clr <- NULL
anova_cap <- NULL
# CLR
if(meta_available || !is.null(covariate_summary)){
  # ====== 8 & 9) CLR matrix + db-RDA ======
  mat_cpm <- class_cpm %>%
    select(sample, arg_class, CPM) %>%
    pivot_wider(names_from = arg_class, values_from = CPM, values_fill = 0) %>%
    arrange(sample) %>%
    filter(sample %in% meta$sample)
  meta_use <- meta %>% filter(sample %in% mat_cpm$sample)
  # ensure sample_base available then attach covariates from masterlist if present
  meta_use <- meta_use %>% mutate(sample_base = sub('_S[0-9]+$','', sample))
  if(exists('meta_cov') && !is.null(meta_cov) && 'sample_base' %in% names(meta_cov)){
    meta_use <- meta_use %>% left_join(meta_cov, by='sample_base')
  }
  # dynamic numeric covariates from joined meta
  numeric_covars <- setdiff(names(meta_use)[vapply(meta_use, is.numeric, logical(1))], c('sample'))
  # Restrict to physicochemical parameters only (whitelist common fields) and drop QC/positional
  physchem_whitelist <- c(
    'cond_m_s_cm','conductivity','temperature','temp_c','turbidity','turbidity_ppt','salinity','salinity_ppt','p_h','ph','do','do_mg_l'
  )
  drop_non_phys <- c('conc_ng_u_l_de_novix_ds_11','a260_280_de_novix_ds_11','a260_230_de_novix_ds_11','longitude','latitude','draw_order')
  # Map potential names in meta to whitelist keys
  nc_lower <- tolower(numeric_covars)
  keep_idx <- nc_lower %in% physchem_whitelist & !(numeric_covars %in% drop_non_phys)
  numeric_covars <- numeric_covars[keep_idx]
  if(length(numeric_covars) > 0){
    nz_vars <- vapply(numeric_covars, function(col){
      v <- meta_use[[col]]; stats::var(v, na.rm=TRUE) > 0 && sum(!is.na(v)) >= 5
    }, logical(1))
    numeric_covars <- numeric_covars[nz_vars]
  }
  if(nrow(meta_use) > 2 && length(numeric_covars) >= 1){
  rownames(meta_use) <- meta_use$sample
  X <- mat_cpm %>% select(-sample) %>% as.matrix()
  X_repl <- zCompositions::cmultRepl(X, label=0, method="CZM")
  X_clr  <- compositions::clr(compositions::acomp(X_repl)) %>% as.matrix()
  rownames(X_clr) <- mat_cpm$sample
  # build covariate frame and align
  covars <- meta_use %>% select(sample, any_of(numeric_covars))
  covars[numeric_covars] <- lapply(covars[numeric_covars], base::scale)
  # Align samples present in both X_clr and covars
  common <- intersect(rownames(X_clr), covars$sample)
  covars_aligned <- covars %>% filter(sample %in% common)
  X_clr_aligned <- X_clr[common, , drop=FALSE]
  # Drop rows with NA covariates
  cc <- stats::complete.cases(covars_aligned[numeric_covars])
  covars_cc <- covars_aligned[cc, , drop=FALSE]
  X_clr_cc <- X_clr_aligned[cc, , drop=FALSE]
    # Variable selection to avoid overfitting
    # 1) Drop highly correlated covariates (|r|>0.95)
    keep_vars <- numeric_covars
    if(length(keep_vars) > 1){
      cm <- try(stats::cor(covars_cc[keep_vars], use='pairwise.complete.obs'), silent=TRUE)
      if(!inherits(cm,'try-error')){
        to_drop <- character(0)
        ut <- which(upper.tri(cm), arr.ind=TRUE)
        for(k in seq_len(nrow(ut))){
          i <- ut[k,1]; j <- ut[k,2]
          if(abs(cm[i,j]) > 0.95){ to_drop <- unique(c(to_drop, colnames(cm)[j])) }
        }
        keep_vars <- setdiff(keep_vars, to_drop)
      }
    }
    # 2) Limit number of constraints based on samples
    max_k <- max(1, min(8, nrow(X_clr_cc) - 5))
    if(length(keep_vars) > max_k){
      # Rank covariates by correlation with first PC of response
      pc1 <- try(stats::prcomp(X_clr_cc, scale.=FALSE)$x[,1], silent=TRUE)
      if(inherits(pc1,'try-error')) pc1 <- rowMeans(X_clr_cc)
      scor <- vapply(keep_vars, function(v) abs(stats::cor(covars_cc[[v]], pc1, use='pairwise.complete.obs')), numeric(1))
      keep_vars <- names(sort(scor, decreasing=TRUE))[seq_len(max_k)]
    }
    # 3) Fit RDA on CLR (equivalent to Euclidean db-RDA)
    if(length(keep_vars) >= 1 && nrow(X_clr_cc) >= 6){
      df_cov <- covars_cc[, c('sample', keep_vars), drop=FALSE]
      rownames(df_cov) <- df_cov$sample
      df_cov <- df_cov[, keep_vars, drop=FALSE]
      form <- as.formula(paste('X_clr_cc ~', paste(keep_vars, collapse=' + ')))
      rda_fit <- try(vegan::rda(form, data = df_cov), silent=TRUE)
      if(!inherits(rda_fit,'try-error')){
        anova_cap <- try(vegan::anova.cca(rda_fit, permutations = 499), silent=TRUE)
        r2_adj <- try(vegan::RsquareAdj(rda_fit)$adj.r.squared, silent=TRUE)
        if(inherits(r2_adj,'try-error')) r2_adj <- NA_real_
      } else {
        message('RDA failed; skipping db-RDA.')
      }
    } else {
      message('Not enough samples after alignment or no covariates; skipping db-RDA.')
    }
  } else {
    message("Insufficient samples for db-RDA; skipping.")
  }
} else {
  message("Skipping CLR/db-RDA: metadata or required columns absent.")
}

# ====== 10) ILR + Bayesian regression (brms, horseshoe) ======
delta_elpd <- NA_real_
heldout_r2_mean <- NA_real_
driver_rank <- tibble(term=character(), median_abs = numeric())
if(!is.na(brms_backend) && meta_available && exists("X") && !is.null(X) && is.matrix(X)){
  keep_classes <- unique(core$arg_class)
  X_model <- if(length(keep_classes) >= 3) X[, colnames(X) %in% keep_classes, drop=FALSE] else X
  X_model_repl <- zCompositions::cmultRepl(X_model, label=0, method="CZM")
  A <- compositions::acomp(X_model_repl)
  ILR <- as.data.frame(compositions::pivotCoord(A)); ILR$sample <- rownames(X)

  # Build covariate frame consistent with earlier selection
  if (exists("covars_cc")) {
    covars_bayes <- covars_cc %>% mutate(sample = rownames(covars_cc))
    brms_covars  <- setdiff(names(covars_cc), "sample")
  } else if (exists("covars")) {
    covars_bayes <- covars
    brms_covars  <- setdiff(names(covars_bayes), "sample")
  } else {
    covars_bayes <- tibble(sample = ILR$sample)
    brms_covars  <- character(0)
  }

  has_district <- "district" %in% names(meta)
  if (has_district) {
    ddf <- meta %>% select(sample, district) %>% mutate(district = as.factor(district))
    covars_bayes <- covars_bayes %>% left_join(ddf, by="sample")
  }

  dat_model <- ILR %>% left_join(covars_bayes, by="sample")
  coords <- setdiff(names(ILR), "sample")

  if(length(coords) > 0){
    fx <- intersect(brms_covars, names(dat_model))
    rhs <- if(length(fx)>0) paste(fx, collapse=" + ") else "1"
    if (has_district && "district" %in% names(dat_model)) rhs <- paste(rhs, "+ (1|district)")
    hs_prior <- set_prior("horseshoe(1)", class="b") + set_prior("student_t(3,0,2.5)", class="Intercept")

    brms_fits <- list()
    for(co in coords){
      f <- as.formula(glue("{co} ~ {rhs}"))
      dat_co <- dat_model %>% drop_na(all_of(c(co, intersect(c(fx,"district"), names(dat_model)))))
      if(nrow(dat_co) < 10) next
      fit <- brm(formula = f, data = dat_co, family = gaussian(), prior = hs_prior,
                 backend = brms_backend, iter=2000, warmup=1000, chains=4, refresh=0)
      brms_fits[[co]] <- fit
    }

    if(length(brms_fits) > 0){
      rep_co <- names(brms_fits)[1]
      f_full <- formula(brms_fits[[rep_co]])
      null_rhs <- if (has_district && "district" %in% names(dat_model)) "1 + (1|district)" else "1"
      dat_rep <- dat_model %>% drop_na(all_of(c(rep_co, if(has_district) "district")))
      null_fit <- brm(as.formula(glue("{rep_co} ~ {null_rhs}")), data = dat_rep,
                      family = gaussian(), backend = brms_backend,
                      iter=2000, warmup=1000, chains=4, refresh=0)
      loo_full <- loo(brms_fits[[rep_co]]); loo_null <- loo(null_fit)
      delta_elpd <- loo_compare(loo_full, loo_null)[1,"elpd_diff"]

      set.seed(42)
      k <- loo_kfold
      folds <- sample(rep(1:k, length.out = nrow(dat_rep)))
      r2s <- numeric(k)
      for(i in seq_len(k)){
        tr <- dat_rep[folds!=i,]; te <- dat_rep[folds==i,]
        fk <- brm(f_full, data=tr, family=gaussian(), backend=brms_backend,
                  iter=1500, warmup=750, chains=2, refresh=0)
        pred <- posterior_epred(fk, newdata=te) %>% apply(2, mean)
        r2s[i] <- stats::cor(pred, te[[rep_co]], use="complete.obs")^2
      }
      heldout_r2_mean <- mean(r2s)

      coef_summ <- map_df(names(brms_fits), function(co){
        s <- posterior_summary(brms_fits[[co]]); rn <- rownames(s)
        tibble(coord = co, term = rn, est = s[,"Estimate"], l95 = s[,"Q2.5"], u95 = s[,"Q97.5"]) %>%
          filter(str_detect(term, "^b_")) %>% mutate(term = str_remove(term, "^b_"))
      })
      driver_rank <- coef_summ %>% group_by(term) %>% summarise(median_abs = median(abs(est), na.rm=TRUE), .groups='drop') %>% arrange(desc(median_abs)) %>% slice(1:5)
    }
  } else {
    message("No ILR coordinates available for modeling.")
  }
} else {
  message("Skipping Bayesian ILR modeling due to missing Stan backend or inputs.")
}

# ====== 11) ABSTRACT output (newline-safe + sensitivity grid) ======

n_sites <- n_distinct(class_cpm$sample)
n_classes_total <- n_distinct(class_cpm$arg_class)
n_core <- nrow(core)

top3 <- if(nrow(core) > 0){
  core %>% slice_max(prevalence_over_X, n = min(3, nrow(core))) %>%
    transmute(x = glue("{arg_class} ({fmt_pct(prevalence_over_X)})")) %>% pull(x)
} else character()
top3_str <- if(length(top3) > 0) paste(top3, collapse='; ') else "None"

# Core sensitivity grid (vary prevalence & CPM thresholds)
sense_prev <- c(0.50, 0.60, 0.70)
sense_cpm  <- c(0.01, 0.05, 0.10)
sensitivity_grid <- expand_grid(prev_thr = sense_prev, cpm_thr = sense_cpm) %>%
  mutate(core_size = purrr::pmap_int(list(prev_thr, cpm_thr), function(pv, cp){
    tmp <- class_cpm %>%
      group_by(arg_class) %>%
      summarise(n_sites_over_X = sum(CPM >= cp, na.rm=TRUE), .groups='drop') %>%
      mutate(prev = n_sites_over_X / n_sites)
    sum(tmp$prev >= pv)
  })) %>% arrange(desc(prev_thr), cpm_thr)

# Short inline sensitivity summary
sens_str <- sensitivity_grid %>%
  mutate(lbl = glue("{round(prev_thr*100)}%/{cpm_thr}=\n{core_size}")) %>%
  transmute(lbl = glue("{round(prev_thr*100)}%/{cpm_thr}={core_size}")) %>%
  pull(lbl) %>% paste(collapse='; ')

perm_p <- if(!is.null(anova_cap) && !inherits(anova_cap, "try-error")) {
  format.pval(anova_cap$`Pr(>F)`[1], digits = 3)
} else "NA"
drv_str <- driver_rank %>% mutate(term = str_replace_all(term, "_", " ")) %>% transmute(x = glue("{term}")) %>% pull(x)

abstract_lines <- c(
  "--- ABSTRACT DROP-INS ---",
  glue("Sites (after exclusions): {n_sites}"),
  if(exists("excluded_samples_summary")) glue("Excluded samples: {excluded_samples_summary}") else NULL,
  glue("ARG classes with CPM>0: {n_classes_total}"),
  glue("Mapping diagnostic (distinct raw classes pre-core): {diag_n_classes}"),
  if(fallback_neff_used) "CPM note: Using raw ShortBRED Count as CPM (N_eff ignored)" else NULL,
  if(!is.null(covariate_summary)) glue("Covariates used: {covariate_summary}") else NULL,
  if(!is.null(admin_summary_line)) admin_summary_line else NULL,
  if(!is.null(metadata_enrichment_summary)) metadata_enrichment_summary else NULL,
  glue("Core definition (primary): present in ≥{round(100*core_prev_threshold)}% of sites at ≥{cpm_core_threshold} CPM"),
  glue("Core size: {n_core} classes"),
  glue("Top core classes by prevalence: {top3_str}"),
  glue("Core sensitivity (core sizes): {sens_str}"),
  glue("db-RDA (CLR/Aitchison) adjusted R²: {ifelse(is.na(r2_adj),'NA',round(r2_adj,3))}; permutation p-value: {perm_p}"),
  glue("Bayesian ILR regression ΔELPD (rep coord) vs null: {ifelse(is.na(delta_elpd),'NA',round(delta_elpd,2))}"),
  glue("Held-out R² (mean across {loo_kfold} folds): {ifelse(is.na(heldout_r2_mean),'NA',round(heldout_r2_mean,2))}"),
  glue("Top model drivers (median |β| across balances): {if(length(drv_str)==0) 'NA' else paste(drv_str, collapse=', ')}"),
  "--- END ---",
  ""
)

writeLines(abstract_lines)

# ====== 12) OPTIONAL: save key tables ======
dir_create("results")
write_tsv(core, "results/core_classes.tsv")
write_tsv(prev_over_X, "results/class_prevalence_over_X.tsv")
write_tsv(prev_over_X_tpm, "results/class_prevalence_over_X_TPM.tsv")
write_tsv(prev_LOD, "results/class_prevalence_LOD.tsv")
write_tsv(core_LOD, "results/core_classes_LOD.tsv")
write_tsv(driver_rank, "results/top_drivers.tsv")

dir_create("figs")

# ====== BEAUTIFUL FIGURES ======
library(ggplot2)
library(patchwork)
library(viridis)
library(scales)
library(ggrepel)

# Custom theme for publication-quality figures
theme_publication <- function(base_size = 11, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 0.8),
      strip.background = element_rect(fill = "grey95", color = "black", size = 0.8),
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

# Helper function to group multi-drug resistance for visualization
simplify_arg_for_viz <- function(df, col_name = "arg_class") {
  df %>%
    mutate(
      arg_class_display = case_when(
        str_detect(!!sym(col_name), ";") ~ "Multi-drug resistance",
        TRUE ~ !!sym(col_name)
      )
    )
}

# Figure 1: ARG Class Prevalence and Abundance
message("[FIGURES] Creating Figure 1: ARG Class Prevalence...")
if(nrow(prev_over_X) > 0){
  fig1_data <- prev_over_X %>%
    simplify_arg_for_viz() %>%
    group_by(arg_class_display) %>%
    summarise(
      prevalence_over_X = max(prevalence_over_X),
      mean_CPM = sum(mean_CPM),
      .groups = "drop"
    ) %>%
    arrange(desc(prevalence_over_X)) %>%
    slice(1:20) %>%
    mutate(arg_class_display = forcats::fct_reorder(arg_class_display, prevalence_over_X),
           is_core = prevalence_over_X >= core_prev_threshold)
  
  p1a <- ggplot(fig1_data, aes(x = prevalence_over_X, y = arg_class_display, fill = is_core)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = core_prev_threshold, linetype = "dashed", 
               color = "red", size = 0.8, alpha = 0.7) +
    scale_x_continuous(labels = percent_format(), expand = c(0, 0)) +
    scale_fill_manual(values = c("TRUE" = "#2E7D32", "FALSE" = "#757575"),
                      labels = c("TRUE" = "Core", "FALSE" = "Non-core"),
                      name = "Classification") +
    labs(x = "Prevalence (%)", y = NULL,
         title = "Top 20 ARG Classes by Prevalence",
         subtitle = glue("Core threshold: ≥{round(100*core_prev_threshold)}% prevalence at ≥{cpm_core_threshold} CPM")) +
    theme_publication() +
    theme(legend.position = c(0.85, 0.15))
  
  p1b <- ggplot(fig1_data, aes(x = mean_CPM, y = arg_class_display, fill = log10(mean_CPM + 1))) +
    geom_col(width = 0.7) +
    scale_x_log10(labels = comma_format()) +
    scale_fill_viridis_c(option = "plasma", name = "log10(CPM)", guide = "none") +
    labs(x = "Mean CPM (log scale)", y = NULL,
         title = "Mean Abundance") +
    theme_publication()
  
  fig1 <- p1a + p1b + plot_layout(widths = c(2, 1))
  ggsave("figs/Fig1_ARG_prevalence_abundance.png", fig1, width = 12, height = 7, dpi = 600, bg = "white")
  ggsave("figs/Fig1_ARG_prevalence_abundance.pdf", fig1, width = 12, height = 7, bg = "white")
  message("[FIGURES] ✓ Saved Figure 1")
}

# Figure 2: Core Sensitivity Heatmap
message("[FIGURES] Creating Figure 2: Core Sensitivity...")
if(exists("sensitivity_grid") && nrow(sensitivity_grid) > 0){
  p2 <- ggplot(sensitivity_grid, aes(x = factor(cpm_thr), y = factor(prev_thr), fill = core_size)) +
    geom_tile(color = "white", size = 1.5) +
    geom_text(aes(label = core_size), color = "white", fontface = "bold", size = 6) +
    scale_fill_viridis_c(option = "inferno", name = "Core Size\n(# classes)") +
    labs(x = "CPM Threshold", y = "Prevalence Threshold",
         title = "ARG Core Definition Sensitivity Analysis",
         subtitle = "Number of core classes under varying thresholds",
         caption = glue("Primary definition: ≥{round(100*core_prev_threshold)}% prevalence, ≥{cpm_core_threshold} CPM")) +
    theme_publication() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 11, face = "bold"))
  
  ggsave("figs/Fig2_core_sensitivity.png", p2, width = 8, height = 6, dpi = 600, bg = "white")
  ggsave("figs/Fig2_core_sensitivity.pdf", p2, width = 8, height = 6, bg = "white")
  message("[FIGURES] ✓ Saved Figure 2")
}

# Figure 3: ARG Composition Stacked Bar
message("[FIGURES] Creating Figure 3: ARG Composition by Sample...")
if(exists("class_cpm") && nrow(class_cpm) > 0){
  # Group multi-drug resistance first
  comp_data_prep <- class_cpm %>%
    simplify_arg_for_viz() %>%
    group_by(sample, arg_class_display) %>%
    summarise(CPM = sum(CPM), .groups = "drop")
  
  # Get top classes for visualization (after grouping)
  top_classes <- comp_data_prep %>%
    group_by(arg_class_display) %>%
    summarise(total = sum(CPM), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice(1:15) %>%
    pull(arg_class_display)
  
  comp_data <- comp_data_prep %>%
    mutate(arg_class_grouped = ifelse(arg_class_display %in% top_classes, arg_class_display, "Other")) %>%
    group_by(sample, arg_class_grouped) %>%
    summarise(CPM = sum(CPM), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(proportion = CPM / sum(CPM))
  
  # Order samples by total ARG abundance
  sample_order <- comp_data %>%
    group_by(sample) %>%
    summarise(total = sum(CPM), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(sample)
  
  comp_data <- comp_data %>%
    mutate(sample = factor(sample, levels = sample_order))
  
  # Create color palette
  n_colors <- length(unique(comp_data$arg_class_grouped))
  colors <- c(viridis::turbo(n_colors - 1), "grey70")
  names(colors) <- c(setdiff(unique(comp_data$arg_class_grouped), "Other"), "Other")
  
  p3 <- ggplot(comp_data, aes(x = sample, y = proportion, fill = arg_class_grouped)) +
    geom_bar(stat = "identity", width = 0.9) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    scale_fill_manual(values = colors, name = "ARG Class") +
    labs(x = "Sample", y = "Relative Abundance (%)",
         title = "ARG Class Composition Across Samples",
         subtitle = glue("Top 15 classes shown, multi-drug resistances grouped (n={n_sites} samples)")) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          legend.position = "right",
          legend.text = element_text(size = 8))
  
  ggsave("figs/Fig3_ARG_composition.png", p3, width = 14, height = 6, dpi = 600, bg = "white")
  ggsave("figs/Fig3_ARG_composition.pdf", p3, width = 14, height = 6, bg = "white")
  message("[FIGURES] ✓ Saved Figure 3")
}

# Figure 4: Spatial Distribution (if admin mapping exists)
message("[FIGURES] Creating Figure 4: Spatial Patterns...")
if(exists("admin_mapping") && !is.null(admin_mapping) && "longitude" %in% names(admin_mapping)){
  # Join with ARG totals
  arg_totals <- class_cpm %>%
    group_by(sample) %>%
    summarise(total_arg_cpm = sum(CPM), 
              richness = sum(CPM > 0), .groups = "drop")
  
  spatial_data <- admin_mapping %>%
    left_join(arg_totals, by = "sample") %>%
    filter(!is.na(longitude), !is.na(latitude), !is.na(total_arg_cpm))
  
  if(nrow(spatial_data) > 0){
    p4a <- ggplot(spatial_data, aes(x = longitude, y = latitude)) +
      geom_point(aes(size = total_arg_cpm, color = log10(total_arg_cpm + 1)), alpha = 0.7) +
      scale_size_continuous(range = c(2, 12), name = "Total ARG\nCPM") +
      scale_color_viridis_c(option = "magma", name = "log10(CPM)") +
      labs(x = "Longitude", y = "Latitude",
           title = "Spatial Distribution of ARG Abundance",
           subtitle = "Point size and color represent total ARG load") +
      theme_publication() +
      theme(legend.position = "right")
    
    p4b <- ggplot(spatial_data, aes(x = longitude, y = latitude)) +
      geom_point(aes(size = richness, color = richness), alpha = 0.7) +
      scale_size_continuous(range = c(2, 12), name = "ARG\nRichness") +
      scale_color_viridis_c(option = "viridis", name = "# Classes") +
      labs(x = "Longitude", y = "Latitude",
           title = "ARG Richness Distribution") +
      theme_publication() +
      theme(legend.position = "right")
    
    fig4 <- p4a / p4b
    ggsave("figs/Fig4_spatial_patterns.png", fig4, width = 10, height = 12, dpi = 600, bg = "white")
    ggsave("figs/Fig4_spatial_patterns.pdf", fig4, width = 10, height = 12, bg = "white")
    message("[FIGURES] ✓ Saved Figure 4")
  }
}

# Figure 5: db-RDA Ordination (if available)
message("[FIGURES] Creating Figure 5: db-RDA Ordination...")
if(exists("rda_fit") && !inherits(rda_fit, "try-error") && !is.null(rda_fit)){
  # Extract site scores
  site_scores <- try(vegan::scores(rda_fit, display = "sites", choices = 1:2), silent = TRUE)
  
  if(!inherits(site_scores, "try-error")){
    ord_data <- as.data.frame(site_scores)
    ord_data$sample <- rownames(ord_data)
    
    # Join with total ARG for coloring
    ord_data <- ord_data %>%
      left_join(class_cpm %>% group_by(sample) %>% 
                  summarise(total_cpm = sum(CPM), .groups = "drop"), 
                by = "sample")
    
    # Get variance explained
    var_exp <- try(vegan::eigenvals(rda_fit) / sum(vegan::eigenvals(rda_fit)) * 100, silent = TRUE)
    if(!inherits(var_exp, "try-error") && length(var_exp) >= 2){
      xlab <- glue("RDA1 ({round(var_exp[1], 1)}%)")
      ylab <- glue("RDA2 ({round(var_exp[2], 1)}%)")
    } else {
      xlab <- "RDA1"; ylab <- "RDA2"
    }
    
    p5 <- ggplot(ord_data, aes(x = RDA1, y = RDA2)) +
      geom_point(aes(color = log10(total_cpm + 1), size = total_cpm), alpha = 0.7) +
      scale_color_viridis_c(option = "plasma", name = "log10(Total CPM)") +
      scale_size_continuous(range = c(2, 8), name = "Total CPM") +
      labs(x = xlab, y = ylab,
           title = "Distance-based RDA of ARG Composition",
           subtitle = glue("Environmental drivers explain {ifelse(is.na(r2_adj), 'NA', round(100*r2_adj, 1))}% of variance (adjusted R²)"),
           caption = glue("Permutation p-value: {perm_p}")) +
      theme_publication() +
      theme(legend.position = "right")
    
    ggsave("figs/Fig5_dbRDA_ordination.png", p5, width = 10, height = 7, dpi = 600, bg = "white")
    ggsave("figs/Fig5_dbRDA_ordination.pdf", p5, width = 10, height = 7, bg = "white")
    message("[FIGURES] ✓ Saved Figure 5")
  }
}

# Figure 6: Bayesian Model Coefficients (if available)
message("[FIGURES] Creating Figure 6: Bayesian Driver Effects...")
if(exists("brms_fits") && length(brms_fits) > 0){
  # Extract and combine coefficients
  coef_all <- map_df(names(brms_fits), function(co){
    s <- posterior_summary(brms_fits[[co]])
    rn <- rownames(s)
    tibble(coord = co, term = rn, 
           est = s[,"Estimate"], 
           l95 = s[,"Q2.5"], 
           u95 = s[,"Q97.5"]) %>%
      filter(str_detect(term, "^b_"), !str_detect(term, "Intercept"))
  })
  
  if(nrow(coef_all) > 0){
    # Clean term names
    coef_all <- coef_all %>%
      mutate(term = str_remove(term, "^b_"),
             term = str_replace_all(term, "_", " "),
             term = str_to_title(term),
             sig = ifelse(sign(l95) == sign(u95), "Significant", "Not significant"))
    
    # Get top effects
    top_terms <- coef_all %>%
      group_by(term) %>%
      summarise(mean_abs = mean(abs(est)), .groups = "drop") %>%
      arrange(desc(mean_abs)) %>%
      slice(1:10) %>%
      pull(term)
    
    plot_data <- coef_all %>%
      filter(term %in% top_terms) %>%
      mutate(term = forcats::fct_reorder(term, abs(est), .fun = median))
    
    p6 <- ggplot(plot_data, aes(x = est, y = term, color = sig)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_linerange(aes(xmin = l95, xmax = u95), size = 0.8, alpha = 0.7) +
      geom_point(size = 2.5, alpha = 0.9) +
      scale_color_manual(values = c("Significant" = "#D32F2F", "Not significant" = "#757575"),
                        name = "95% CI") +
      labs(x = "Standardized Coefficient", y = NULL,
           title = "Environmental Driver Effects on ARG Composition",
           subtitle = "Bayesian horseshoe regression coefficients (ILR coordinates)",
           caption = glue("ΔELPD vs null: {ifelse(is.na(delta_elpd), 'NA', round(delta_elpd, 1))} | Held-out R²: {ifelse(is.na(heldout_r2_mean), 'NA', round(heldout_r2_mean, 2))}")) +
      theme_publication() +
      theme(legend.position = c(0.85, 0.15))
    
    ggsave("figs/Fig6_bayesian_drivers.png", p6, width = 10, height = 7, dpi = 600, bg = "white")
    ggsave("figs/Fig6_bayesian_drivers.pdf", p6, width = 10, height = 7, bg = "white")
    message("[FIGURES] ✓ Saved Figure 6")
  }
}

# Figure 7: Prevalence Comparison (CPM vs LOD)
message("[FIGURES] Creating Figure 7: Prevalence Methods Comparison...")
if(exists("prev_over_X") && exists("prev_LOD")){
  # Group multi-drug resistance
  prev_grouped <- prev_over_X %>%
    simplify_arg_for_viz() %>%
    group_by(arg_class_display) %>%
    summarise(prevalence_over_X = max(prevalence_over_X), .groups = "drop")
  
  prev_lod_grouped <- prev_LOD %>%
    simplify_arg_for_viz() %>%
    group_by(arg_class_display) %>%
    summarise(prevalence_LOD = max(prevalence_LOD), .groups = "drop")
  
  comp_prev <- prev_grouped %>%
    left_join(prev_lod_grouped, by = "arg_class_display") %>%
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
  
  # Label high-prevalence and core classes
  labels_to_show <- comp_prev %>%
    filter(is_core_CPM | is_core_LOD | prev_CPM >= 0.5 | prev_LOD >= 0.5) %>%
    distinct(arg_class_display, .keep_all = TRUE)
  
  p7 <- ggplot(comp_prev, aes(x = prev_CPM, y = prev_LOD, color = category)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = core_prev_threshold, linetype = "dotted", color = "red", alpha = 0.5) +
    geom_vline(xintercept = core_prev_threshold, linetype = "dotted", color = "red", alpha = 0.5) +
    geom_point(size = 3, alpha = 0.7) +
    ggrepel::geom_label_repel(
      data = labels_to_show,
      aes(label = arg_class_display),
      size = 3,
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
         caption = glue("Core threshold: {round(100*core_prev_threshold)}% | LOD: ≥{lod_min_reads} reads | n={nrow(comp_prev)} classes")) +
    theme_publication() +
    theme(legend.position = c(0.2, 0.75),
          legend.background = element_rect(fill = "white", color = "black"))
  
  ggsave("figs/Fig7_prevalence_methods.png", p7, width = 11, height = 9, dpi = 600, bg = "white")
  ggsave("figs/Fig7_prevalence_methods.pdf", p7, width = 11, height = 9, bg = "white")
  message("[FIGURES] ✓ Saved Figure 7")
}

message("[FIGURES] ═══════════════════════════════════════")
message("[FIGURES] All publication-quality figures created!")
message("[FIGURES] Location: figs/ directory")
message("[FIGURES] Formats: PNG (300 dpi) and PDF")
message("[FIGURES] ═══════════════════════════════════════")

# ====== 13) Pathogen–ARG analysis (integrated) ======
try({
  merged_path <- "diversity_results/metaphlan_bracken_merged_both.tsv"
  pathogen_list <- "clean_pathogens.txt"
  if(fs::file_exists(merged_path) && fs::file_exists(pathogen_list)){
    message("[PATHOGEN] Building pathogen tables and CLR correlations…")
    merged <- suppressMessages(readr::read_tsv(merged_path, show_col_types = FALSE)) %>% janitor::clean_names()
    if(!all(c("sample","taxon","rank") %in% names(merged))){
      stop("Merged taxonomy file missing required columns: sample, taxon, rank")
    }
    if(!"bracken" %in% names(merged)) merged$bracken <- NA_real_
    if(!"metaphlan" %in% names(merged)) merged$metaphlan <- NA_real_
    merged <- merged %>% mutate(abund = dplyr::coalesce(bracken, metaphlan, 0))

    pat_lines <- readLines(pathogen_list, warn = FALSE)
    pat_genera <- tibble(line = pat_lines) %>%
      filter(nchar(trimws(line)) > 0) %>%
      mutate(genus = stringr::str_extract(line, "^[A-Za-z]+")) %>%
      filter(!is.na(genus)) %>% distinct(genus)
    message("[PATHOGEN] Unique pathogen genera in list: ", nrow(pat_genera))

    genus_tbl <- merged %>%
      mutate(genus = dplyr::case_when(
        stringr::str_detect(taxon, '^g__') ~ stringr::str_remove(taxon, '^g__'),
        tolower(rank) == 'genus' ~ as.character(taxon),
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(genus)) %>%
      transmute(sample, genus, value = abund) %>%
      group_by(sample, genus) %>%
      summarise(value = max(replace_na(value, 0), na.rm=TRUE), .groups='drop')

    P_long <- genus_tbl %>% filter(genus %in% pat_genera$genus) %>%
      mutate(present = value > 0) %>%
      group_by(sample, genus) %>% summarise(present = any(present),
                                            value = max(value, na.rm=TRUE), .groups='drop')

    metrics <- P_long %>% group_by(sample) %>% summarise(
      pathogen_pct = sum(value, na.rm=TRUE),
      pathogen_richness = sum(present, na.rm=TRUE), .groups='drop') %>% arrange(sample)

    pat_wide <- genus_tbl %>% filter(genus %in% pat_genera$genus) %>%
      tidyr::pivot_wider(names_from = genus, values_from = value, values_fill = 0) %>%
      arrange(sample)

    readr::write_tsv(metrics, "results/pathogen_sample_metrics.tsv")
    readr::write_tsv(pat_wide, "results/pathogen_genus_wide.tsv")
    message("[PATHOGEN] Wrote pathogen tables: metrics + genus wide.")

    # Join with ARG totals (already computed above) and wide class CPM
    arg_total_df <- class_cpm %>% group_by(sample) %>% summarise(
      arg_total_cpm = sum(CPM, na.rm=TRUE), arg_total_tpm = sum(TPM, na.rm=TRUE), .groups='drop')
    arg_wide_df <- class_cpm %>% select(sample, arg_class, CPM) %>%
      tidyr::pivot_wider(names_from = arg_class, values_from = CPM, values_fill = 0) %>%
      arrange(sample)

    df_merge <- metrics %>% left_join(arg_total_df, by='sample') %>% left_join(arg_wide_df, by='sample')
    if(nrow(df_merge) >= 5){
      # Sample-level associations
      sp1 <- try(stats::cor.test(df_merge$pathogen_pct, df_merge$arg_total_cpm, method='spearman'), silent=TRUE)
      sp2 <- try(stats::cor.test(df_merge$pathogen_richness, df_merge$arg_total_cpm, method='spearman'), silent=TRUE)
      if(!inherits(sp1,'try-error')) message("[PATHOGEN] Spearman pathogen_pct ~ ARG_total_cpm: rho=", round(sp1$estimate,3), " p=", signif(sp1$p.value,3))
      if(!inherits(sp2,'try-error')) message("[PATHOGEN] Spearman pathogen_richness ~ ARG_total_cpm: rho=", round(sp2$estimate,3), " p=", signif(sp2$p.value,3))
    }

    # CLR co-structure (Aitchison distances) and pairwise CLR correlations
    common <- intersect(pat_wide$sample, arg_wide_df$sample)
    P_mat <- pat_wide %>% filter(sample %in% common) %>% column_to_rownames('sample') %>% as.matrix()
    A_mat <- arg_wide_df %>% filter(sample %in% common) %>% column_to_rownames('sample') %>% as.matrix()

    # Drop zero-variance columns
    A_mat <- A_mat[, apply(A_mat, 2, function(v) stats::var(v) > 0), drop=FALSE]
    P_mat <- P_mat[, apply(P_mat, 2, function(v) stats::var(v) > 0), drop=FALSE]

    if(nrow(A_mat) >= 4 && nrow(P_mat) >= 4 && nrow(A_mat) == nrow(P_mat)){
      A_repl <- zCompositions::cmultRepl(A_mat, label=0, method='CZM')
      P_repl <- zCompositions::cmultRepl(P_mat, label=0, method='CZM')
      A_clr  <- compositions::clr(compositions::acomp(A_repl)) %>% as.matrix()
      P_clr  <- compositions::clr(compositions::acomp(P_repl)) %>% as.matrix()

      # Mantel-like via vegan::mantel (Spearman on Euclidean distances of CLR)
      DA <- dist(A_clr)
      DP <- dist(P_clr)
      mt <- try(vegan::mantel(DA, DP, method='spearman', permutations = 999), silent=TRUE)
      if(!inherits(mt,'try-error')) message('[PATHOGEN] CLR Mantel (Spearman) ARG vs pathogen: r=', round(mt$statistic,3), ' p=', signif(mt$signif,3))

      # Pairwise CLR correlations
      arg_cols <- colnames(A_clr); gen_cols <- colnames(P_clr)
      res_list <- list(); idx <- 1L
      for(g in gen_cols){
        for(cn in arg_cols){
          x <- P_clr[, g]; y <- A_clr[, cn]
          if(stats::var(x) == 0 || stats::var(y) == 0) next
          ct <- try(stats::cor.test(x, y, method='pearson'), silent=TRUE)
          if(inherits(ct,'try-error')) next
          res_list[[idx]] <- tibble(genus = g, arg_class = cn, r_clr = unname(ct$estimate), p = ct$p.value)
          idx <- idx + 1L
        }
      }
      if(length(res_list) > 0){
        pairs <- dplyr::bind_rows(res_list)
        pairs <- pairs %>% mutate(q = p.adjust(p, method='BH')) %>% arrange(q, desc(r_clr))
        readr::write_tsv(pairs %>% filter(q <= 0.1), 'results/pathogen_genus__arg_class_clr_pairs.tsv')
        message('[PATHOGEN] Pairwise CLR correlations written: ', sum(pairs$q <= 0.1), ' significant pairs (q<=0.1).')
      } else {
        message('[PATHOGEN] No pairwise CLR correlations computed (empty column sets).')
      }
    } else {
      message('[PATHOGEN] Skipping CLR analysis: insufficient or misaligned samples/columns.')
    }
  } else {
    message('[PATHOGEN] Skipped: required inputs missing (', merged_path, ' or ', pathogen_list, ').')
  }
}, silent=TRUE)

# ====== PATHOGEN-ARG FIGURES ======
message("[FIGURES] Creating pathogen-ARG correlation figures...")
try({
  if(fs::file_exists("results/pathogen_sample_metrics.tsv") && 
     fs::file_exists("results/arg_total_by_sample.tsv")){
    
    pat_metrics <- readr::read_tsv("results/pathogen_sample_metrics.tsv", show_col_types = FALSE)
    arg_metrics <- readr::read_tsv("results/arg_total_by_sample.tsv", show_col_types = FALSE)
    
    combined <- pat_metrics %>% 
      left_join(arg_metrics, by = "sample") %>%
      filter(!is.na(arg_total_cpm))
    
    if(nrow(combined) >= 5){
      # Figure 8: Pathogen vs ARG correlations
      p8a <- ggplot(combined, aes(x = pathogen_pct, y = arg_total_cpm)) +
        geom_point(aes(color = pathogen_richness), size = 3, alpha = 0.7) +
        geom_smooth(method = "lm", color = "red", linetype = "dashed", se = TRUE, alpha = 0.2) +
        scale_color_viridis_c(option = "mako", name = "Pathogen\nRichness") +
        scale_y_log10(labels = comma_format()) +
        labs(x = "Pathogen Abundance (%)", y = "Total ARG CPM (log scale)",
             title = "Pathogen Abundance vs. ARG Load",
             subtitle = "Potential co-occurrence patterns") +
        theme_publication()
      
      p8b <- ggplot(combined, aes(x = pathogen_richness, y = arg_total_cpm)) +
        geom_point(aes(color = pathogen_pct), size = 3, alpha = 0.7) +
        geom_smooth(method = "lm", color = "blue", linetype = "dashed", se = TRUE, alpha = 0.2) +
        scale_color_viridis_c(option = "rocket", name = "Pathogen\nAbund (%)") +
        scale_y_log10(labels = comma_format()) +
        labs(x = "Pathogen Richness (# genera)", y = "Total ARG CPM (log scale)",
             title = "Pathogen Richness vs. ARG Load") +
        theme_publication()
      
      fig8 <- p8a + p8b + plot_layout(guides = "collect")
      ggsave("figs/Fig8_pathogen_ARG_correlation.png", fig8, width = 14, height = 6, dpi = 600, bg = "white")
      ggsave("figs/Fig8_pathogen_ARG_correlation.pdf", fig8, width = 14, height = 6, bg = "white")
      message("[FIGURES] ✓ Saved Figure 8")
    }
  }
  
  # Figure 9: Top pathogen-ARG CLR pairs
  if(fs::file_exists("results/pathogen_genus__arg_class_clr_pairs.tsv")){
    clr_pairs <- readr::read_tsv("results/pathogen_genus__arg_class_clr_pairs.tsv", show_col_types = FALSE)
    
    if(nrow(clr_pairs) > 0){
      top_pairs <- clr_pairs %>%
        arrange(q) %>%
        slice(1:min(30, nrow(clr_pairs))) %>%
        mutate(pair = glue("{genus} - {arg_class}"),
               pair = forcats::fct_reorder(pair, r_clr))
      
      p9 <- ggplot(top_pairs, aes(x = r_clr, y = pair, fill = r_clr)) +
        geom_col(width = 0.7) +
        geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                            midpoint = 0, name = "CLR\nCorrelation") +
        labs(x = "CLR Correlation Coefficient", y = NULL,
             title = "Top Pathogen Genus - ARG Class Associations",
             subtitle = "Compositionally-aware (CLR) correlations",
             caption = glue("Showing top {nrow(top_pairs)} pairs with q ≤ 0.1 (FDR-corrected)")) +
        theme_publication() +
        theme(axis.text.y = element_text(size = 8))
      
      ggsave("figs/Fig9_pathogen_ARG_pairs.png", p9, width = 10, height = 10, dpi = 600, bg = "white")
      ggsave("figs/Fig9_pathogen_ARG_pairs.pdf", p9, width = 10, height = 10, bg = "white")
      message("[FIGURES] ✓ Saved Figure 9")
    }
  }
}, silent = TRUE)

# Figure 10: Summary Dashboard
message("[FIGURES] Creating summary dashboard...")
try({
  # Create a multi-panel summary figure
  
  # Panel A: Core vs Non-core
  if(exists("core") && nrow(core) > 0){
    summary_data <- prev_over_X %>%
      mutate(status = ifelse(arg_class %in% core$arg_class, "Core", "Non-core"))
    
    pA <- ggplot(summary_data, aes(x = status, fill = status)) +
      geom_bar(stat = "count", width = 0.6) +
      geom_text(stat = "count", aes(label = after_stat(count)), 
                vjust = -0.5, fontface = "bold", size = 5) +
      scale_fill_manual(values = c("Core" = "#388E3C", "Non-core" = "#9E9E9E"), 
                       guide = "none") +
      labs(x = NULL, y = "Number of Classes",
           title = "A. Core Classification",
           subtitle = glue("{nrow(core)} core / {nrow(prev_over_X)} total classes")) +
      theme_publication() +
      theme(axis.text.x = element_text(size = 12, face = "bold"))
  } else {
    pA <- ggplot() + theme_void()
  }
  
  # Panel B: Sample distribution
  if(exists("class_cpm")){
    sample_totals <- class_cpm %>%
      group_by(sample) %>%
      summarise(total = sum(CPM), .groups = "drop")
    
    pB <- ggplot(sample_totals, aes(x = total)) +
      geom_histogram(bins = 30, fill = "#1976D2", color = "white", alpha = 0.8) +
      geom_vline(aes(xintercept = median(total)), color = "red", 
                linetype = "dashed", size = 1) +
      scale_x_log10(labels = comma_format()) +
      labs(x = "Total ARG CPM (log scale)", y = "Number of Samples",
           title = "B. ARG Load Distribution",
           subtitle = glue("Median: {round(median(sample_totals$total), 1)} CPM")) +
      theme_publication()
  } else {
    pB <- ggplot() + theme_void()
  }
  
  # Panel C: Richness
  if(exists("class_cpm")){
    richness <- class_cpm %>%
      filter(CPM > 0) %>%
      group_by(sample) %>%
      summarise(richness = n(), .groups = "drop")
    
    pC <- ggplot(richness, aes(x = richness)) +
      geom_histogram(bins = 20, fill = "#7B1FA2", color = "white", alpha = 0.8) +
      geom_vline(aes(xintercept = mean(richness)), color = "orange", 
                linetype = "dashed", size = 1) +
      labs(x = "ARG Class Richness", y = "Number of Samples",
           title = "C. ARG Richness Distribution",
           subtitle = glue("Mean: {round(mean(richness$richness), 1)} classes/sample")) +
      theme_publication()
  } else {
    pC <- ggplot() + theme_void()
  }
  
  # Panel D: Model performance
  if(!is.na(r2_adj) || !is.na(heldout_r2_mean)){
    model_data <- tibble(
      metric = c("db-RDA\nAdj. R²", "Bayesian\nHeld-out R²"),
      value = c(ifelse(is.na(r2_adj), 0, r2_adj), 
                ifelse(is.na(heldout_r2_mean), 0, heldout_r2_mean)),
      available = c(!is.na(r2_adj), !is.na(heldout_r2_mean))
    )
    
    pD <- ggplot(model_data %>% filter(available), 
                 aes(x = metric, y = value, fill = metric)) +
      geom_col(width = 0.6, alpha = 0.8) +
      geom_text(aes(label = sprintf("%.2f", value)), 
                vjust = -0.5, fontface = "bold", size = 5) +
      scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
      scale_fill_manual(values = c("#FF6F00", "#00897B"), guide = "none") +
      labs(x = NULL, y = "Variance Explained",
           title = "D. Model Performance",
           subtitle = "Compositional analysis results") +
      theme_publication() +
      theme(axis.text.x = element_text(size = 10, face = "bold"))
  } else {
    pD <- ggplot() + theme_void()
  }
  
  # Combine panels
  fig10 <- (pA + pB) / (pC + pD) +
    plot_annotation(
      title = "ARG Analysis Summary Dashboard",
      subtitle = glue("Dataset: {n_sites} samples | {n_classes_total} ARG classes detected"),
      theme = theme_publication()
    )
  
  ggsave("figs/Fig10_summary_dashboard.png", fig10, width = 14, height = 10, dpi = 600, bg = "white")
  ggsave("figs/Fig10_summary_dashboard.pdf", fig10, width = 14, height = 10, bg = "white")
  message("[FIGURES] ✓ Saved Figure 10 (Summary Dashboard)")
}, silent = TRUE)

message("\n")
message("╔════════════════════════════════════════════════════════════╗")
message("║  ✨ PUBLICATION-QUALITY FIGURES COMPLETE ✨                 ║")
message("╠════════════════════════════════════════════════════════════╣")
message("║  📊 Up to 10 figures generated in figs/ directory          ║")
message("║  📁 Formats: High-res PNG (300 dpi) + Vector PDF           ║")
message("║                                                            ║")
message("║  Figure 1: ARG Class Prevalence & Abundance                ║")
message("║  Figure 2: Core Definition Sensitivity Heatmap             ║")
message("║  Figure 3: Sample Composition Stacked Bars                 ║")
message("║  Figure 4: Spatial Distribution Maps                       ║")
message("║  Figure 5: db-RDA Ordination Plot                          ║")
message("║  Figure 6: Bayesian Driver Effects                         ║")
message("║  Figure 7: Prevalence Method Comparison                    ║")
message("║  Figure 8: Pathogen-ARG Correlations                       ║")
message("║  Figure 9: Top Pathogen-ARG Pairs (CLR)                    ║")
message("║  Figure 10: Summary Dashboard (4-panel)                    ║")
message("╚════════════════════════════════════════════════════════════╝")
message("\n")
