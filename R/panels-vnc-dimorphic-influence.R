###########################################################
### Sensory → Dimorphic/Specific Influence Heatmaps
###
### Rows:    sensory cell_function × body_part_sensory
###          (only groups present in both BANC & maleCNS)
### Columns: dimorphic + specific cell types per super_class
###          (split: cross-mapped dimorphic | other-dimorphic | specific)
###          Cross-mapped: same cell type names in both heatmaps
###          (uses malecns_cell_type convention for BANC)
###
### Min-max normalised across both heatmaps (per super_class set)
###
### Three sets of paired heatmaps (BANC female / maleCNS male):
###   1. ventral_nerve_cord_intrinsic
###   2. ascending
###   3. descending
###
### Inspired by connectome_influence_maps/07_publication_heatmaps.R
###########################################################

.this_script <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) normalizePath(f) else stop("Cannot determine script path")
})
.this_dir <- dirname(.this_script)
.this_repo <- normalizePath(file.path(.this_dir, ".."))

source(file.path(.this_dir, "startup.R"))

suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dendextend)
})

# ── Paths ─────────────────────────────────────────────────
output_dir <- file.path(.this_repo, "figs/figure_dimorphic/links")
cache_dir  <- file.path(.this_repo, "data/banc")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

banc_cache_file   <- file.path(cache_dir, "banc_sensory_all_dim_influence.csv")
malecns_cache_file <- file.path(cache_dir, "malecns_sensory_all_dim_influence.csv")
force_recompute <- TRUE  # always recalculate — influence parameters likely to change

# ── Influence palette (blue → mauve → red) ───────────────
influence_palette <- colorRampPalette(
  c("#1f4e79", "#4a90a4", "#7ba7bc", "#a67c8a", "#c4967d", "#b22222")
)

# ═══════════════════════════════════════════════════════════
# 1. LOAD METADATA
# ═══════════════════════════════════════════════════════════
message("=== Loading metadata ===")

# BANC metadata (banc.meta loaded by startup)
# Ensure sexually_dimorphic + malecns_cell_type are present
banc_cols_needed <- c("sexually_dimorphic", "malecns_cell_type")
missing_cols <- setdiff(banc_cols_needed, names(banc.meta))
if (length(missing_cols) > 0) {
  message("  Fetching missing columns from SeaTable: ", paste(missing_cols, collapse = ", "))
  banc_extra <- banctable_query(paste0(
    "SELECT root_id, ", paste(banc_cols_needed, collapse = ", "), " FROM banc_meta"
  )) %>%
    dplyr::mutate(root_id = as.character(root_id)) %>%
    dplyr::distinct(root_id, .keep_all = TRUE)
  banc.meta <- banc.meta %>%
    dplyr::left_join(banc_extra %>% dplyr::select(root_id, dplyr::all_of(missing_cols)),
                     by = "root_id")
}

# maleCNS metadata already loaded by startup.R (malecns.meta, malecns.edgelist.simple)

# maleCNS sexually_dimorphic from SeaTable
message("  Querying maleCNS sexually_dimorphic from SeaTable...")
malecns_sd <- tryCatch({
  banctable_query(
    "SELECT malecns_09_id, sexually_dimorphic FROM malecns",
    base = "cns_meta"
  ) %>%
    dplyr::filter(!is.na(sexually_dimorphic), sexually_dimorphic != "") %>%
    dplyr::mutate(malecns_09_id = as.character(malecns_09_id)) %>%
    dplyr::distinct(malecns_09_id, .keep_all = TRUE)
}, error = function(e) {
  message("  Could not query maleCNS SeaTable: ", e$message)
  NULL
})

if (!is.null(malecns_sd) && nrow(malecns_sd) > 0) {
  malecns.meta <- malecns.meta %>%
    dplyr::left_join(malecns_sd, by = "malecns_09_id")
  message(sprintf("  Joined %d maleCNS sexually_dimorphic labels", nrow(malecns_sd)))
}

message(sprintf("  BANC: %d neurons | maleCNS: %d neurons",
                nrow(banc.meta), nrow(malecns.meta)))

# ═══════════════════════════════════════════════════════════
# 2. BUILD TYPE-LEVEL DIMORPHISM LOOKUP + TYPE BRIDGE
# ═══════════════════════════════════════════════════════════
message("=== Identifying dimorphic/sex-specific targets (all super_classes) ===")

# -- BANC dimorphic/sex-specific types (all super_classes) --
banc_dim_all <- banc.meta %>%
  dplyr::filter(sexually_dimorphic %in% c("dimorphic", "female-specific", "male-specific"),
                !is.na(cell_type), cell_type != "",
                !grepl("^AN_", cell_type))  # drop improperly matched AN_ types

banc_type_dim <- banc_dim_all %>%
  dplyr::distinct(cell_type, sexually_dimorphic)

message(sprintf("  BANC: %d dimorphic/specific neurons, %d unique types",
                nrow(banc_dim_all), nrow(banc_type_dim)))

# -- Type bridge: BANC cell_type <-> maleCNS cell_type --
banc_bridge <- banc.meta %>%
  dplyr::filter(!is.na(malecns_cell_type), malecns_cell_type != "") %>%
  dplyr::distinct(cell_type, malecns_cell_type)

# -- maleCNS dimorphic/sex-specific VNC intrinsic targets --
has_mc_dim <- "sexually_dimorphic" %in% names(malecns.meta) &&
  sum(!is.na(malecns.meta$sexually_dimorphic) & malecns.meta$sexually_dimorphic != "") > 0

if (has_mc_dim) {
  message("  Using maleCNS sexually_dimorphic labels from SeaTable")
  malecns_dim_all <- malecns.meta %>%
    dplyr::filter(sexually_dimorphic %in% c("dimorphic", "female-specific", "male-specific"),
                  !is.na(cell_type), cell_type != "")
} else {
  message("  Deriving maleCNS dimorphism from BANC cell_type matching")
  type_map_for_derive <- banc_type_dim %>%
    dplyr::left_join(banc_bridge, by = "cell_type") %>%
    dplyr::mutate(mc_type = dplyr::coalesce(malecns_cell_type, cell_type)) %>%
    dplyr::select(mc_type, sexually_dimorphic)
  malecns_dim_all <- malecns.meta %>%
    dplyr::filter(cell_type %in% type_map_for_derive$mc_type,
                  !is.na(cell_type), cell_type != "") %>%
    dplyr::left_join(type_map_for_derive, by = c("cell_type" = "mc_type")) %>%
    dplyr::filter(!is.na(sexually_dimorphic))
}

malecns_type_dim <- malecns_dim_all %>%
  dplyr::distinct(cell_type, sexually_dimorphic)

message(sprintf("  maleCNS: %d dimorphic/specific neurons, %d unique types",
                nrow(malecns_dim_all), nrow(malecns_type_dim)))

# -- Determine cross-mapping between datasets --
# For each BANC dimorphic type, find the maleCNS equivalent
banc_dim_only <- banc_type_dim %>%
  dplyr::filter(sexually_dimorphic == "dimorphic") %>%
  dplyr::left_join(banc_bridge, by = "cell_type") %>%
  dplyr::mutate(mc_type = dplyr::coalesce(malecns_cell_type, cell_type))

malecns_dim_type_set <- malecns_type_dim %>%
  dplyr::filter(sexually_dimorphic == "dimorphic") %>%
  dplyr::pull(cell_type) %>% unique()

# Cross-mapped: BANC dimorphic type whose maleCNS equivalent is also dimorphic
cross_mapped <- banc_dim_only %>%
  dplyr::filter(mc_type %in% malecns_dim_type_set) %>%
  dplyr::distinct(banc_type = cell_type, mc_type)

message(sprintf("  Cross-mapped dimorphic types: %d", nrow(cross_mapped)))

# -- Assign 3-level categories --
# BANC: cross-mapped dimorphic uses malecns_cell_type name for consistency
banc_dim_all <- banc_dim_all %>%
  dplyr::left_join(cross_mapped, by = c("cell_type" = "banc_type")) %>%
  dplyr::mutate(
    dim_category = dplyr::case_when(
      sexually_dimorphic %in% c("female-specific", "male-specific") ~ "specific",
      !is.na(mc_type) ~ "cross-mapped dimorphic",
      TRUE ~ "other-dimorphic"
    ),
    target_type = dplyr::if_else(dim_category == "cross-mapped dimorphic",
                                  mc_type, cell_type)
  )

# maleCNS: cross-mapped types already use maleCNS cell_type names
malecns_dim_all <- malecns_dim_all %>%
  dplyr::mutate(
    dim_category = dplyr::case_when(
      sexually_dimorphic %in% c("female-specific", "male-specific") ~ "specific",
      cell_type %in% cross_mapped$mc_type ~ "cross-mapped dimorphic",
      TRUE ~ "other-dimorphic"
    ),
    target_type = cell_type
  )

# Summary
banc_cats <- table(banc_dim_all %>% dplyr::distinct(target_type, dim_category) %>%
                     dplyr::pull(dim_category))
mc_cats <- table(malecns_dim_all %>% dplyr::distinct(target_type, dim_category) %>%
                   dplyr::pull(dim_category))
message(sprintf("  BANC types: %s",
                paste(names(banc_cats), banc_cats, sep = "=", collapse = ", ")))
message(sprintf("  maleCNS types: %s",
                paste(names(mc_cats), mc_cats, sep = "=", collapse = ", ")))

# -- Super_class breakdown --
message("  BANC super_class breakdown:")
for (sc in sort(unique(banc_dim_all$super_class))) {
  n <- sum(banc_dim_all$super_class == sc, na.rm = TRUE)
  message(sprintf("    %s: %d neurons", sc, n))
}

# ═══════════════════════════════════════════════════════════
# 3. IDENTIFY SENSORY NEURONS + SEX-SPECIFIC GROUPS
# ═══════════════════════════════════════════════════════════
message("=== Identifying sensory neuron groups ===")

make_sensory_groups <- function(meta, id_col) {
  meta %>%
    dplyr::filter(grepl("sensory", super_class, ignore.case = TRUE)) %>%
    dplyr::filter(!grepl("unknown|orphan", cell_type, ignore.case = TRUE)) %>%
    dplyr::filter(!grepl("^JO", cell_type)) %>%  # drop Johnston's organ (poorly reconstructed)
    dplyr::filter(!grepl(",| ", cell_type)) %>%
    dplyr::filter(!is.na(body_part_sensory), body_part_sensory != "") %>%
    dplyr::mutate(
      function_label = dplyr::coalesce(cell_function_detailed, cell_function),
      source_group = paste0(function_label, " [", body_part_sensory, "]")
    ) %>%
    dplyr::filter(!is.na(function_label)) %>%
    dplyr::select(id = dplyr::all_of(id_col), cell_type, source_group,
                  function_label, body_part_sensory)
}

banc_sensory   <- make_sensory_groups(banc.meta, "root_id")
malecns_sensory <- make_sensory_groups(malecns.meta, "malecns_09_id")

# Track all groups (no filtering — compute influence for all)
all_sensory_groups <- union(unique(banc_sensory$source_group),
                            unique(malecns_sensory$source_group))
common_groups <- intersect(unique(banc_sensory$source_group),
                           unique(malecns_sensory$source_group))

# Identify sex-specific sensory groups (from BANC metadata)
sex_specific_groups <- banc.meta %>%
  dplyr::filter(grepl("sensory", super_class, ignore.case = TRUE),
                sexually_dimorphic %in% c("female-specific", "male-specific"),
                !is.na(body_part_sensory), body_part_sensory != "") %>%
  dplyr::mutate(
    function_label = dplyr::coalesce(cell_function_detailed, cell_function),
    source_group = paste0(function_label, " [", body_part_sensory, "]")
  ) %>%
  dplyr::filter(!is.na(function_label)) %>%
  dplyr::distinct(source_group) %>%
  dplyr::pull(source_group)

message(sprintf("  All sensory groups: %d (common: %d, sex-specific: %d)",
                length(all_sensory_groups), length(common_groups),
                length(sex_specific_groups)))
message(sprintf("  BANC: %d neurons | maleCNS: %d neurons",
                nrow(banc_sensory), nrow(malecns_sensory)))

# ═══════════════════════════════════════════════════════════
# 4. COMPUTE BANC INFLUENCE (via influence_calculator_py)
# ═══════════════════════════════════════════════════════════

if (!file.exists(banc_cache_file) || force_recompute) {
  message("=== Computing BANC influence ===")

  # Use BANC edgelist from startup (already loaded by banc-edgelist.R)
  if (!exists("banc.edgelist.simple") || is.null(banc.edgelist.simple)) {
    stop("banc.edgelist.simple not loaded — ensure banc-startup.R ran successfully")
  }
  message(sprintf("  Using BANC edgelist: %d rows", nrow(banc.edgelist.simple)))

  # Build influence calculator
  message("  Building BANC influence calculator...")
  ic_banc <- influence_calculator_py(
    edgelist_simple = banc.edgelist.simple %>% dplyr::filter(count >= 3),
    meta = banc.meta
  )
  message("  Influence calculator ready")

  # Source sensory combos
  banc_combos <- banc_sensory %>%
    dplyr::distinct(cell_type, body_part_sensory, source_group)

  # Target IDs and type mapping
  banc_target_map <- banc_dim_all %>%
    dplyr::distinct(root_id, target_type, dim_category, super_class)
  banc_target_ids <- unique(banc_target_map$root_id)

  message(sprintf("  Computing influence for %d sensory combos → %d targets...",
                  nrow(banc_combos), length(banc_target_ids)))
  banc_results <- list()

  for (i in seq_len(nrow(banc_combos))) {
    ct <- banc_combos$cell_type[i]
    bp <- banc_combos$body_part_sensory[i]
    sg <- banc_combos$source_group[i]

    source_ids <- banc_sensory$id[
      banc_sensory$cell_type == ct &
        banc_sensory$body_part_sensory == bp]

    if (length(source_ids) == 0) next

    if (i %% 25 == 0 || i == 1 || i == nrow(banc_combos)) {
      message(sprintf("    [%d/%d] %s (%d seeds)",
                      i, nrow(banc_combos), sg, length(source_ids)))
    }

    tryCatch({
      inf <- calculate_influence_py(ic_banc, source_ids)
      inf_filtered <- inf %>%
        dplyr::filter(id %in% banc_target_ids) %>%
        dplyr::mutate(
          source_group = sg,
          body_part_sensory = bp,
          raw_influence = `Influence_score_(unsigned)`
        ) %>%
        dplyr::select(id, source_group, body_part_sensory, raw_influence)

      if (nrow(inf_filtered) > 0) {
        banc_results[[length(banc_results) + 1]] <- inf_filtered
      }
    }, error = function(e) {
      message(sprintf("    Error for %s: %s", sg, e$message))
    })
  }

  banc_inf_raw <- dplyr::bind_rows(banc_results)
  message(sprintf("  Collected %d source-target pairs", nrow(banc_inf_raw)))

  # Aggregate by target_type: median influence per source_group × target_type
  banc_agg <- banc_inf_raw %>%
    dplyr::inner_join(banc_target_map, by = c("id" = "root_id")) %>%
    dplyr::group_by(source_group, body_part_sensory, target_type, dim_category, super_class) %>%
    dplyr::summarise(
      median_influence = median(raw_influence, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      adjusted_influence = pmax(0, log(median_influence) + 24),
      adjusted_influence = dplyr::if_else(
        is.infinite(adjusted_influence) | is.nan(adjusted_influence), 0, adjusted_influence)
    )

  readr::write_csv(banc_agg, banc_cache_file)
  message(sprintf("  Cached BANC influence: %d rows → %s", nrow(banc_agg), banc_cache_file))

  rm(ic_banc, banc_results, banc_inf_raw, banc.edgelist.simple); gc()
} else {
  message("=== Loading cached BANC influence ===")
  banc_agg <- readr::read_csv(banc_cache_file, show_col_types = FALSE)
  message(sprintf("  Loaded %d rows from cache", nrow(banc_agg)))
}

# Drop improperly matched AN_ types from BANC aggregation
banc_agg <- banc_agg %>% dplyr::filter(!grepl("^AN_", target_type))

# ═══════════════════════════════════════════════════════════
# 5. COMPUTE maleCNS INFLUENCE (influence_calculator_py)
# ═══════════════════════════════════════════════════════════

if (!file.exists(malecns_cache_file) || force_recompute) {
  message("=== Computing maleCNS influence (this may take a while) ===")

  # Load maleCNS edgelist
  message("  Loading maleCNS edgelist from GCS...")
  malecns.edgelist.simple <- read_feather_gcs(
    file.path(malecns.gcs.path, "malecns_09_simple_edgelist.feather"))

  # Build influence calculator
  message("  Building maleCNS influence calculator...")
  ic_malecns <- influence_calculator_py(
    edgelist_simple = malecns.edgelist.simple %>% dplyr::filter(count > 0),
    meta = malecns.meta %>% dplyr::rename(root_id = malecns_09_id)
  )
  message("  Influence calculator ready")

  # Unique sensory combos
  malecns_combos <- malecns_sensory %>%
    dplyr::distinct(cell_type, body_part_sensory, source_group)

  # Target IDs
  malecns_target_ids <- unique(malecns_dim_all$malecns_09_id)

  # Target type map: malecns_09_id → (target_type, dim_category, super_class)
  mc_target_map <- malecns_dim_all %>%
    dplyr::distinct(malecns_09_id, target_type, dim_category, super_class)

  # Compute influence for each sensory combo
  message(sprintf("  Computing influence for %d sensory combos...", nrow(malecns_combos)))
  malecns_results <- list()

  for (i in seq_len(nrow(malecns_combos))) {
    ct <- malecns_combos$cell_type[i]
    bp <- malecns_combos$body_part_sensory[i]
    sg <- malecns_combos$source_group[i]

    source_ids <- malecns_sensory$id[
      malecns_sensory$cell_type == ct &
        malecns_sensory$body_part_sensory == bp]

    if (length(source_ids) == 0) next

    if (i %% 25 == 0 || i == 1 || i == nrow(malecns_combos)) {
      message(sprintf("    [%d/%d] %s (%d seeds)",
                      i, nrow(malecns_combos), sg, length(source_ids)))
    }

    tryCatch({
      inf <- calculate_influence_py(ic_malecns, source_ids)
      inf_filtered <- inf %>%
        dplyr::filter(id %in% malecns_target_ids) %>%
        dplyr::mutate(
          source_group = sg,
          body_part_sensory = bp,
          raw_influence = `Influence_score_(unsigned)`
        ) %>%
        dplyr::select(id, source_group, body_part_sensory, raw_influence)

      if (nrow(inf_filtered) > 0) {
        malecns_results[[length(malecns_results) + 1]] <- inf_filtered
      }
    }, error = function(e) {
      message(sprintf("    Error for %s: %s", sg, e$message))
    })
  }

  malecns_inf_raw <- dplyr::bind_rows(malecns_results)
  message(sprintf("  Collected %d source-target pairs", nrow(malecns_inf_raw)))

  # Aggregate by target_type: median influence per source_group × target_type
  malecns_agg <- malecns_inf_raw %>%
    dplyr::inner_join(mc_target_map, by = c("id" = "malecns_09_id")) %>%
    dplyr::group_by(source_group, body_part_sensory, target_type, dim_category, super_class) %>%
    dplyr::summarise(
      median_influence = median(raw_influence, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      adjusted_influence = pmax(0, log(median_influence) + 24),
      adjusted_influence = dplyr::if_else(
        is.infinite(adjusted_influence) | is.nan(adjusted_influence), 0, adjusted_influence)
    )

  readr::write_csv(malecns_agg, malecns_cache_file)
  message(sprintf("  Cached maleCNS influence: %d rows → %s",
                  nrow(malecns_agg), malecns_cache_file))

  rm(ic_malecns, malecns_results, malecns_inf_raw, malecns.edgelist.simple); gc()
} else {
  message("=== Loading cached maleCNS influence ===")
  malecns_agg <- readr::read_csv(malecns_cache_file, show_col_types = FALSE)
  message(sprintf("  Loaded %d rows from cache", nrow(malecns_agg)))
}

# ═══════════════════════════════════════════════════════════
# 6-9. BUILD MATRICES + HEATMAPS (per super_class group)
# ═══════════════════════════════════════════════════════════

# --- Shared helper functions ---

build_matrix <- function(agg_df, target_types) {
  agg_df %>%
    dplyr::filter(target_type %in% target_types) %>%
    dplyr::select(source_group, target_type, adjusted_influence) %>%
    dplyr::group_by(source_group, target_type) %>%
    dplyr::summarise(adjusted_influence = mean(adjusted_influence, na.rm = TRUE),
                     .groups = "drop") %>%
    tidyr::pivot_wider(names_from = target_type, values_from = adjusted_influence,
                       values_fill = list(adjusted_influence = 0)) %>%
    tibble::column_to_rownames("source_group") %>%
    as.matrix()
}

pad_matrix <- function(mat, rows, cols = colnames(mat)) {
  out <- matrix(0, nrow = length(rows), ncol = length(cols),
                dimnames = list(rows, cols))
  shared_r <- intersect(rows, rownames(mat))
  shared_c <- intersect(cols, colnames(mat))
  if (length(shared_r) > 0 && length(shared_c) > 0) {
    out[shared_r, shared_c] <- mat[shared_r, shared_c]
  }
  out
}

cluster_cols <- function(mat) {
  if (ncol(mat) <= 1) colnames(mat)
  else {
    # Replace NA with 0 for clustering only
    mat_cl <- mat; mat_cl[is.na(mat_cl)] <- 0
    hc <- hclust(dist(t(mat_cl)), method = "ward.D2")
    colnames(mat)[hc$order]
  }
}

extract_body_part <- function(sg) {
  bp <- gsub(".*\\[(.*)\\]", "\\1", sg)
  ifelse(bp == "NA" | bp == sg, NA_character_, bp)
}

minmax_norm <- function(mat, lo, hi) {
  if (hi == lo) { mat[!is.na(mat)] <- 0; return(mat) }
  (mat - lo) / (hi - lo)  # NA stays NA
}

colwise_minmax_norm <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    vals <- mat[, j]
    non_na <- !is.na(vals)
    if (sum(non_na) == 0) next
    lo <- quantile(vals[non_na], 0.01)
    hi <- quantile(vals[non_na], 0.99)
    if (hi > lo) {
      mat[non_na, j] <- pmin(pmax((vals[non_na] - lo) / (hi - lo), 0), 1)
    } else {
      mat[non_na, j] <- 0
    }
  }
  mat
}

# Category colors (2 levels — other-dimorphic dropped)
cat_levels <- c("cross-mapped dimorphic", "specific")
dim_colors <- c("cross-mapped dimorphic" = "#984EA3",
                "specific" = "#E41A1C")

# Fixed cell size (mm) for cross-plot comparability
cell_w_mm <- 3
cell_h_mm <- 4

# Normalisation modes (raw removed — only normalised plots needed)
norm_modes <- list(
  list(tag = "minmax", title_suffix = "matrix min-max normalised",
       legend_title = "Influence\n(min-max)"),
  list(tag = "colnorm", title_suffix = "column min-max normalised",
       legend_title = "Influence\n(col-norm)")
)

# --- Super_class groups to plot ---
sc_groups <- list(
  list(sc_filter = "ventral_nerve_cord_intrinsic",
       label = "VNC Intrinsic",
       file_tag = "vnci"),
  list(sc_filter = "ascending",
       label = "Ascending",
       file_tag = "ascending"),
  list(sc_filter = "descending",
       label = "Descending",
       file_tag = "descending")
)

for (scg in sc_groups) {
  message(sprintf("\n=== Heatmaps: %s ===", scg$label))

  # Filter aggregated data to this super_class (exclude other-dimorphic)
  banc_agg_sc <- banc_agg %>%
    dplyr::filter(super_class == scg$sc_filter, dim_category != "other-dimorphic")
  malecns_agg_sc <- malecns_agg %>%
    dplyr::filter(super_class == scg$sc_filter, dim_category != "other-dimorphic")

  # Filter dim targets to this super_class (exclude other-dimorphic)
  banc_dim_sc <- banc_dim_all %>%
    dplyr::filter(super_class == scg$sc_filter, dim_category != "other-dimorphic")
  malecns_dim_sc <- malecns_dim_all %>%
    dplyr::filter(super_class == scg$sc_filter, dim_category != "other-dimorphic")

  if (nrow(banc_dim_sc) == 0 && nrow(malecns_dim_sc) == 0) {
    message(sprintf("  No dimorphic/specific %s neurons in either dataset; skipping", scg$label))
    next
  }

  # Cross-mapped types for this super_class
  banc_cross_sc <- banc_dim_sc %>%
    dplyr::filter(dim_category == "cross-mapped dimorphic") %>%
    dplyr::distinct(target_type) %>% dplyr::pull(target_type)
  mc_cross_sc <- malecns_dim_sc %>%
    dplyr::filter(dim_category == "cross-mapped dimorphic") %>%
    dplyr::distinct(target_type) %>% dplyr::pull(target_type)
  cross_types_sc_all <- union(banc_cross_sc, mc_cross_sc)
  banc_has_data_sc <- unique(banc_agg_sc$target_type)
  mc_has_data_sc <- unique(malecns_agg_sc$target_type)
  # Keep only dimorphic types present in BOTH datasets
  cross_types_sc <- cross_types_sc_all[cross_types_sc_all %in% banc_has_data_sc &
                                        cross_types_sc_all %in% mc_has_data_sc]
  cross_banc_only_sc <- cross_types_sc_all[cross_types_sc_all %in% banc_has_data_sc &
                                            !cross_types_sc_all %in% mc_has_data_sc]
  cross_mc_only_sc <- cross_types_sc_all[!cross_types_sc_all %in% banc_has_data_sc &
                                          cross_types_sc_all %in% mc_has_data_sc]

  if (length(cross_banc_only_sc) > 0 || length(cross_mc_only_sc) > 0) {
    disc_df <- dplyr::bind_rows(
      if (length(cross_banc_only_sc) > 0)
        data.frame(target_type = cross_banc_only_sc, dataset = "BANC_only",
                   super_class = scg$sc_filter, stringsAsFactors = FALSE),
      if (length(cross_mc_only_sc) > 0)
        data.frame(target_type = cross_mc_only_sc, dataset = "maleCNS_only",
                   super_class = scg$sc_filter, stringsAsFactors = FALSE)
    )
    disc_file <- file.path(output_dir, sprintf("dimorphic_discrepancies_%s.csv", scg$file_tag))
    readr::write_csv(disc_df, disc_file)
    message(sprintf("  Discrepancies: %d BANC-only, %d maleCNS-only → %s",
                    length(cross_banc_only_sc), length(cross_mc_only_sc), basename(disc_file)))
  }

  banc_spec_sc <- banc_dim_sc %>%
    dplyr::filter(dim_category == "specific") %>%
    dplyr::distinct(target_type) %>% dplyr::pull(target_type) %>%
    intersect(unique(banc_agg_sc$target_type))

  mc_spec_sc <- malecns_dim_sc %>%
    dplyr::filter(dim_category == "specific") %>%
    dplyr::distinct(target_type) %>% dplyr::pull(target_type) %>%
    intersect(unique(malecns_agg_sc$target_type))

  n_total_cols <- length(cross_types_sc) + length(banc_spec_sc) + length(mc_spec_sc)
  if (n_total_cols == 0) {
    message(sprintf("  No target types with influence data for %s; skipping", scg$label))
    next
  }

  message(sprintf("  Columns — cross-mapped: %d | BANC specific: %d | maleCNS specific: %d",
                  length(cross_types_sc), length(banc_spec_sc), length(mc_spec_sc)))

  # Build sub-matrices
  banc_mat_cross  <- build_matrix(banc_agg_sc, cross_types_sc)
  banc_mat_spec   <- build_matrix(banc_agg_sc, banc_spec_sc)
  mc_mat_cross    <- build_matrix(malecns_agg_sc, cross_types_sc)
  mc_mat_spec     <- build_matrix(malecns_agg_sc, mc_spec_sc)

  # All rows from both datasets
  all_rows <- Reduce(union, lapply(
    list(banc_mat_cross, banc_mat_spec, mc_mat_cross, mc_mat_spec),
    rownames))
  all_rows <- intersect(all_rows, all_sensory_groups)

  banc_mat_cross <- pad_matrix(banc_mat_cross, all_rows, cross_types_sc)
  banc_mat_spec  <- pad_matrix(banc_mat_spec, all_rows, banc_spec_sc)
  mc_mat_cross   <- pad_matrix(mc_mat_cross, all_rows, cross_types_sc)
  mc_mat_spec    <- pad_matrix(mc_mat_spec, all_rows, mc_spec_sc)

  # Filter weak rows (50th percentile threshold)
  combined_maxes <- pmax(
    apply(cbind(banc_mat_cross, banc_mat_spec), 1, max, na.rm = TRUE),
    apply(cbind(mc_mat_cross, mc_mat_spec), 1, max, na.rm = TRUE)
  )
  threshold <- quantile(combined_maxes, 0.75, na.rm = TRUE)
  keep_rows <- names(combined_maxes)[combined_maxes >= threshold]

  if (length(keep_rows) == 0) {
    message(sprintf("  No rows above threshold for %s; skipping", scg$label))
    next
  }

  banc_mat_cross <- banc_mat_cross[keep_rows, , drop = FALSE]
  banc_mat_spec  <- banc_mat_spec[keep_rows, , drop = FALSE]
  mc_mat_cross   <- mc_mat_cross[keep_rows, , drop = FALSE]
  mc_mat_spec    <- mc_mat_spec[keep_rows, , drop = FALSE]

  message(sprintf("  After threshold filter: %d rows (threshold = %.2f)",
                  length(keep_rows), threshold))

  # Column ordering
  cross_order      <- cluster_cols(banc_mat_cross)
  banc_spec_order  <- cluster_cols(banc_mat_spec)
  mc_spec_order    <- cluster_cols(mc_mat_spec)

  banc_mat_full <- cbind(
    banc_mat_cross[, cross_order, drop = FALSE],
    banc_mat_spec[, banc_spec_order, drop = FALSE]
  )
  mc_mat_full <- cbind(
    mc_mat_cross[, cross_order, drop = FALSE],
    mc_mat_spec[, mc_spec_order, drop = FALSE]
  )

  banc_col_split <- factor(
    c(rep("cross-mapped dimorphic", length(cross_order)),
      rep("specific", length(banc_spec_order))),
    levels = cat_levels)
  mc_col_split <- factor(
    c(rep("cross-mapped dimorphic", length(cross_order)),
      rep("specific", length(mc_spec_order))),
    levels = cat_levels)

  # Row ordering (BANC clustering within body_part)
  body_parts <- sapply(rownames(banc_mat_full), extract_body_part)
  body_parts[is.na(body_parts)] <- "unknown"
  bp_counts <- sort(table(body_parts), decreasing = TRUE)
  row_split <- factor(body_parts, levels = names(bp_counts))

  row_order_final <- character(0)
  for (bp in levels(row_split)) {
    bp_rows <- names(row_split)[row_split == bp]
    if (length(bp_rows) <= 1) {
      row_order_final <- c(row_order_final, bp_rows)
    } else {
      bp_mat <- banc_mat_full[bp_rows, , drop = FALSE]
      hc_bp <- hclust(dist(bp_mat), method = "ward.D2")
      row_order_final <- c(row_order_final, bp_rows[hc_bp$order])
    }
  }

  banc_mat_full <- banc_mat_full[row_order_final, , drop = FALSE]
  mc_mat_full   <- mc_mat_full[row_order_final, , drop = FALSE]
  row_split     <- row_split[row_order_final]

  # --- Handle absent rows ---
  # Rows all-zero in one heatmap: drop unless sex-specific (→ grey out)
  banc_row_zero <- rownames(banc_mat_full)[rowSums(banc_mat_full) == 0]
  mc_row_zero   <- rownames(mc_mat_full)[rowSums(mc_mat_full) == 0]
  one_sided <- union(
    setdiff(banc_row_zero, mc_row_zero),
    setdiff(mc_row_zero, banc_row_zero)
  )
  drop_rows <- setdiff(one_sided, sex_specific_groups)
  sex_spec_grey <- intersect(one_sided, sex_specific_groups)

  if (length(drop_rows) > 0) {
    keep <- setdiff(rownames(banc_mat_full), drop_rows)
    banc_mat_full <- banc_mat_full[keep, , drop = FALSE]
    mc_mat_full   <- mc_mat_full[keep, , drop = FALSE]
    row_split     <- row_split[keep]
    message(sprintf("  Dropped %d non-sex-specific absent rows", length(drop_rows)))
  }
  if (length(sex_spec_grey) > 0) {
    for (r in sex_spec_grey) {
      if (r %in% rownames(banc_mat_full) && sum(banc_mat_full[r, ]) == 0)
        banc_mat_full[r, ] <- NA
      if (r %in% rownames(mc_mat_full) && sum(mc_mat_full[r, ]) == 0)
        mc_mat_full[r, ] <- NA
    }
    message(sprintf("  Greyed out %d sex-specific absent rows", length(sex_spec_grey)))
  }

  # --- Grey out all-zero columns (missing data → NA) ---
  for (j in seq_len(ncol(banc_mat_full))) {
    if (all(banc_mat_full[, j] == 0, na.rm = TRUE)) banc_mat_full[, j] <- NA
  }
  for (j in seq_len(ncol(mc_mat_full))) {
    if (all(mc_mat_full[, j] == 0, na.rm = TRUE)) mc_mat_full[, j] <- NA
  }

  n_rows <- nrow(banc_mat_full)
  message(sprintf("  BANC matrix: %d × %d | maleCNS matrix: %d × %d",
                  nrow(banc_mat_full), ncol(banc_mat_full),
                  nrow(mc_mat_full), ncol(mc_mat_full)))

  # Row annotation (shared across normalisation modes)
  unique_bps <- levels(row_split)
  bp_colors <- setNames(rainbow(length(unique_bps), s = 0.6, v = 0.8), unique_bps)

  # --- Loop over normalisation modes ---
  for (nm in norm_modes) {
    message(sprintf("  -- %s --", nm$title_suffix))

    if (nm$tag == "minmax") {
      all_vals <- c(as.vector(banc_mat_full[!is.na(banc_mat_full)]),
                    as.vector(mc_mat_full[!is.na(mc_mat_full)]))
      lo <- quantile(all_vals, 0.01, na.rm = TRUE)
      hi <- quantile(all_vals, 0.99, na.rm = TRUE)
      message(sprintf("    Min-max range (1st-99th pctl): %.4f to %.4f", lo, hi))
      banc_plot <- minmax_norm(banc_mat_full, lo, hi)
      mc_plot   <- minmax_norm(mc_mat_full, lo, hi)
      banc_plot[!is.na(banc_plot)] <- pmin(pmax(banc_plot[!is.na(banc_plot)], 0), 1)
      mc_plot[!is.na(mc_plot)]     <- pmin(pmax(mc_plot[!is.na(mc_plot)], 0), 1)
      col_fun <- colorRamp2(seq(0, 1, length.out = 200), influence_palette(200))
    } else {  # colnorm
      banc_plot <- colwise_minmax_norm(banc_mat_full)
      mc_plot   <- colwise_minmax_norm(mc_mat_full)
      col_fun <- colorRamp2(seq(0, 1, length.out = 200), influence_palette(200))
    }

    # --- Individual heatmaps (as before) ---
    for (ds in list(
      list(mat = banc_plot, cs = banc_col_split, ds_label = "BANC"),
      list(mat = mc_plot, cs = mc_col_split, ds_label = "maleCNS")
    )) {
      n_cols <- ncol(ds$mat)
      fontsize_col <- max(5, min(9, 160 / n_cols))
      fontsize_row <- max(5, min(9, 160 / n_rows))

      col_ha <- HeatmapAnnotation(
        Category = as.character(ds$cs),
        col = list(Category = dim_colors),
        show_legend = TRUE,
        annotation_legend_param = list(Category = list(title = "Category")),
        show_annotation_name = FALSE)

      row_ha <- rowAnnotation(
        body_part = body_parts[rownames(ds$mat)],
        col = list(body_part = bp_colors),
        show_annotation_name = FALSE,
        show_legend = FALSE,
        na_col = "grey90")

      ht <- Heatmap(ds$mat,
        name = nm$legend_title,
        col = col_fun,
        na_col = "grey90",
        width = unit(n_cols * cell_w_mm, "mm"),
        height = unit(n_rows * cell_h_mm, "mm"),
        cluster_columns = FALSE, cluster_rows = FALSE,
        row_split = row_split,
        row_gap = unit(1.5, "mm"),
        row_title_gp = gpar(fontsize = 7, fontface = "bold"),
        row_title_rot = 0,
        column_split = ds$cs,
        column_gap = unit(3, "mm"),
        column_title_gp = gpar(fontsize = 9, fontface = "bold"),
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = fontsize_col),
        row_names_gp = gpar(fontsize = fontsize_row),
        row_names_max_width = unit(12, "cm"),
        top_annotation = col_ha,
        left_annotation = row_ha,
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 8, fontface = "bold"),
          labels_gp = gpar(fontsize = 7)),
        border = FALSE,
        rect_gp = gpar(col = NA))

      body_w_in <- n_cols * cell_w_mm / 25.4
      body_h_in <- n_rows * cell_h_mm / 25.4
      w <- body_w_in + 7
      h <- body_h_in + 4

      title_str <- sprintf("Sensory → Dimorphic/Specific %s (%s, %s)",
                            scg$label, ds$ds_label, nm$title_suffix)
      ds_tag <- tolower(ds$ds_label)

      for (ext in c("png", "pdf")) {
        fname <- file.path(output_dir,
          sprintf("sensory_to_dimorphic_influence_%s_%s_%s.%s",
                  nm$tag, scg$file_tag, ds_tag, ext))
        if (ext == "png") {
          png(fname, width = w, height = h, units = "in", res = 300)
        } else {
          pdf(fname, width = w, height = h)
        }
        draw(ht,
             column_title = title_str,
             column_title_gp = gpar(fontsize = 12, fontface = "bold"),
             padding = unit(c(2, 2, 2, 15), "mm"),
             merge_legend = TRUE)
        dev.off()
      }
      message(sprintf("    Saved: %s (%.1f × %.1f in)", basename(fname), w, h))
    }

    # --- Side-by-side paired heatmap (cross-mapped columns only) ---
    # This enables direct visual comparison of female vs male influence patterns
    if (length(cross_order) > 0) {
      message("    Building side-by-side paired heatmap (cross-mapped columns)...")

      banc_cross_plot <- banc_plot[, seq_len(length(cross_order)), drop = FALSE]
      mc_cross_plot   <- mc_plot[, seq_len(length(cross_order)), drop = FALSE]

      n_cross <- length(cross_order)
      fontsize_col_paired <- max(5, min(9, 160 / n_cross))
      fontsize_row_paired <- max(5, min(9, 160 / n_rows))

      # Difference matrix (BANC - maleCNS) for cross-mapped columns
      diff_mat <- banc_cross_plot - mc_cross_plot
      # Where either is NA, diff is NA
      diff_mat[is.na(banc_cross_plot) | is.na(mc_cross_plot)] <- NA
      diff_max <- max(abs(diff_mat), na.rm = TRUE)
      if (diff_max == 0) diff_max <- 1
      diff_col_fun <- colorRamp2(
        c(-diff_max, 0, diff_max),
        c("#2166AC", "white", "#B2182B")  # blue = male-biased, red = female-biased
      )

      row_ha_paired <- rowAnnotation(
        body_part = body_parts[rownames(banc_cross_plot)],
        col = list(body_part = bp_colors),
        show_annotation_name = FALSE,
        show_legend = TRUE,
        na_col = "grey90")

      ht_banc <- Heatmap(banc_cross_plot,
        name = paste0(nm$legend_title, "\n(shared)"),
        col = col_fun, na_col = "grey90",
        width = unit(n_cross * cell_w_mm, "mm"),
        height = unit(n_rows * cell_h_mm, "mm"),
        cluster_columns = FALSE, cluster_rows = FALSE,
        row_split = row_split,
        row_gap = unit(1.5, "mm"),
        row_title_gp = gpar(fontsize = 7, fontface = "bold"),
        row_title_rot = 0,
        column_title = "BANC (female)",
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = fontsize_col_paired),
        row_names_gp = gpar(fontsize = fontsize_row_paired),
        row_names_max_width = unit(10, "cm"),
        left_annotation = row_ha_paired,
        show_row_names = TRUE,
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 8, fontface = "bold"),
          labels_gp = gpar(fontsize = 7)),
        border = FALSE,
        rect_gp = gpar(col = NA))

      ht_mc <- Heatmap(mc_cross_plot,
        name = " ",  # suppress duplicate legend
        col = col_fun, na_col = "grey90",
        width = unit(n_cross * cell_w_mm, "mm"),
        height = unit(n_rows * cell_h_mm, "mm"),
        cluster_columns = FALSE, cluster_rows = FALSE,
        row_split = row_split,
        row_gap = unit(1.5, "mm"),
        column_title = "maleCNS (male)",
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = fontsize_col_paired),
        show_row_names = FALSE,
        show_heatmap_legend = FALSE,
        border = FALSE,
        rect_gp = gpar(col = NA))

      ht_diff <- Heatmap(diff_mat,
        name = "Difference\n(F - M)",
        col = diff_col_fun, na_col = "grey90",
        width = unit(n_cross * cell_w_mm, "mm"),
        height = unit(n_rows * cell_h_mm, "mm"),
        cluster_columns = FALSE, cluster_rows = FALSE,
        row_split = row_split,
        row_gap = unit(1.5, "mm"),
        column_title = "Difference (F - M)",
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = fontsize_col_paired),
        show_row_names = FALSE,
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 8, fontface = "bold"),
          labels_gp = gpar(fontsize = 7)),
        border = FALSE,
        rect_gp = gpar(col = NA))

      ht_list <- ht_banc + ht_mc + ht_diff

      body_w_paired <- 3 * n_cross * cell_w_mm / 25.4
      body_h_paired <- n_rows * cell_h_mm / 25.4
      w_paired <- body_w_paired + 12  # extra margin for row labels + legends + gaps
      h_paired <- body_h_paired + 4

      title_paired <- sprintf("Sensory → Cross-mapped Dimorphic %s (female vs male, %s)",
                               scg$label, nm$title_suffix)

      for (ext in c("png", "pdf")) {
        fname_paired <- file.path(output_dir,
          sprintf("sensory_dimorphic_paired_%s_%s.%s",
                  nm$tag, scg$file_tag, ext))
        if (ext == "png") {
          png(fname_paired, width = w_paired, height = h_paired, units = "in", res = 300)
        } else {
          pdf(fname_paired, width = w_paired, height = h_paired)
        }
        draw(ht_list,
             column_title = title_paired,
             column_title_gp = gpar(fontsize = 13, fontface = "bold"),
             padding = unit(c(2, 2, 2, 15), "mm"),
             merge_legend = FALSE)
        dev.off()
      }
      message(sprintf("    Saved paired: %s (%.1f × %.1f in)",
                      basename(fname_paired), w_paired, h_paired))
    }
  }
}

# ═══════════════════════════════════════════════════════════
# 10. CONDENSED PAIRED HEATMAPS
# ═══════════════════════════════════════════════════════════
# Fixes:
# - Row names always shown (source_group = function [body_part])
# - Rows included if they pass percentile filter in EITHER dataset
# - Only columns (cell types) present in BOTH datasets are shown;
#   discrepancies (types in one but not the other) saved to CSV
# - VNC intrinsic: select columns by most divergent sensory profiles
#   between BANC and maleCNS (not merged clusters)
# - Two percentile thresholds for VNC intrinsic: 50th and 85th
# - Column-normalised, side-by-side BANC | maleCNS | Difference

message("\n=== Condensed paired heatmaps ===")

# Helper: build and draw a condensed paired heatmap
draw_condensed_paired <- function(banc_mat, mc_mat, row_split_c, body_parts_c,
                                   scg_label, pctl_tag, output_dir) {
  # Column-normalise
  banc_plot <- colwise_minmax_norm(banc_mat)
  mc_plot   <- colwise_minmax_norm(mc_mat)

  # Difference
  diff_mat <- banc_plot - mc_plot
  diff_mat[is.na(banc_plot) | is.na(mc_plot)] <- NA
  diff_max <- max(abs(diff_mat), na.rm = TRUE)
  if (diff_max == 0) diff_max <- 1

  col_fun_c <- colorRamp2(seq(0, 1, length.out = 200), influence_palette(200))
  diff_col_fun_c <- colorRamp2(c(-diff_max, 0, diff_max),
                                c("#2166AC", "white", "#B2182B"))

  n_cols_c <- ncol(banc_plot)
  n_rows_c <- nrow(banc_plot)
  fontsize_col_c <- max(5, min(9, 200 / n_cols_c))
  fontsize_row_c <- max(5, min(9, 200 / n_rows_c))

  unique_bps_c <- levels(row_split_c)
  bp_colors_c <- setNames(rainbow(length(unique_bps_c), s = 0.6, v = 0.8), unique_bps_c)

  row_ha_c <- rowAnnotation(
    body_part = body_parts_c[rownames(banc_plot)],
    col = list(body_part = bp_colors_c),
    show_annotation_name = FALSE,
    show_legend = TRUE,
    na_col = "grey90")

  ht_banc_c <- Heatmap(banc_plot,
    name = "Influence\n(col-norm)",
    col = col_fun_c, na_col = "grey90",
    width = unit(n_cols_c * cell_w_mm, "mm"),
    height = unit(n_rows_c * cell_h_mm, "mm"),
    cluster_columns = FALSE, cluster_rows = FALSE,
    row_split = row_split_c,
    row_gap = unit(1.5, "mm"),
    row_title_gp = gpar(fontsize = 7, fontface = "bold"),
    row_title_rot = 0,
    column_title = "BANC (female)",
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = fontsize_col_c),
    row_names_gp = gpar(fontsize = fontsize_row_c),
    row_names_max_width = unit(12, "cm"),
    left_annotation = row_ha_c,
    show_row_names = TRUE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)),
    border = FALSE,
    rect_gp = gpar(col = NA))

  ht_mc_c <- Heatmap(mc_plot,
    name = " ",
    col = col_fun_c, na_col = "grey90",
    width = unit(n_cols_c * cell_w_mm, "mm"),
    height = unit(n_rows_c * cell_h_mm, "mm"),
    cluster_columns = FALSE, cluster_rows = FALSE,
    row_split = row_split_c,
    row_gap = unit(1.5, "mm"),
    column_title = "maleCNS (male)",
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = fontsize_col_c),
    show_row_names = FALSE,
    show_heatmap_legend = FALSE,
    border = FALSE,
    rect_gp = gpar(col = NA))

  ht_diff_c <- Heatmap(diff_mat,
    name = "Difference\n(F - M)",
    col = diff_col_fun_c, na_col = "grey90",
    width = unit(n_cols_c * cell_w_mm, "mm"),
    height = unit(n_rows_c * cell_h_mm, "mm"),
    cluster_columns = FALSE, cluster_rows = FALSE,
    row_split = row_split_c,
    row_gap = unit(1.5, "mm"),
    column_title = "Difference (F - M)",
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = fontsize_col_c),
    show_row_names = FALSE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)),
    border = FALSE,
    rect_gp = gpar(col = NA))

  ht_list_c <- ht_banc_c + ht_mc_c + ht_diff_c

  body_w_c <- 3 * n_cols_c * cell_w_mm / 25.4
  body_h_c <- n_rows_c * cell_h_mm / 25.4
  w_c <- body_w_c + 14
  h_c <- body_h_c + 4

  title_c <- sprintf("Sensory → Cross-mapped Dimorphic %s (female vs male, %s)",
                      scg_label, pctl_tag)

  for (ext in c("png", "pdf")) {
    fname_c <- file.path(output_dir,
      sprintf("sensory_dimorphic_condensed_%s_%s.%s",
              gsub(" ", "_", pctl_tag), gsub(" ", "_", tolower(scg_label)), ext))
    if (ext == "png") {
      png(fname_c, width = w_c, height = h_c, units = "in", res = 300)
    } else {
      pdf(fname_c, width = w_c, height = h_c)
    }
    draw(ht_list_c,
         column_title = title_c,
         column_title_gp = gpar(fontsize = 13, fontface = "bold"),
         padding = unit(c(2, 2, 2, 15), "mm"),
         merge_legend = FALSE)
    dev.off()
  }
  message(sprintf("  Saved: %s (%.1f × %.1f in, %d rows × %d cols)",
                  basename(fname_c), w_c, h_c, n_rows_c, n_cols_c))
}

# Helper: order rows by body_part clustering
order_rows_by_bp <- function(mat) {
  body_parts_c <- sapply(rownames(mat), extract_body_part)
  body_parts_c[is.na(body_parts_c)] <- "unknown"
  bp_counts_c <- sort(table(body_parts_c), decreasing = TRUE)
  row_split_c <- factor(body_parts_c, levels = names(bp_counts_c))

  row_order_c <- character(0)
  for (bp in levels(row_split_c)) {
    bp_rows <- names(row_split_c)[row_split_c == bp]
    if (length(bp_rows) <= 1) {
      row_order_c <- c(row_order_c, bp_rows)
    } else {
      bp_m <- mat[bp_rows, , drop = FALSE]
      hc_bp <- hclust(dist(bp_m), method = "ward.D2")
      row_order_c <- c(row_order_c, bp_rows[hc_bp$order])
    }
  }
  list(order = row_order_c, split = row_split_c[row_order_c],
       body_parts = body_parts_c)
}

for (scg in sc_groups) {
  message(sprintf("\n--- Condensed: %s ---", scg$label))

  banc_agg_sc <- banc_agg %>%
    dplyr::filter(super_class == scg$sc_filter, dim_category != "other-dimorphic")
  malecns_agg_sc <- malecns_agg %>%
    dplyr::filter(super_class == scg$sc_filter, dim_category != "other-dimorphic")

  banc_dim_sc <- banc_dim_all %>%
    dplyr::filter(super_class == scg$sc_filter, dim_category != "other-dimorphic")
  malecns_dim_sc <- malecns_dim_all %>%
    dplyr::filter(super_class == scg$sc_filter, dim_category != "other-dimorphic")

  # Cross-mapped types: only those present in BOTH datasets
  banc_cross_sc <- banc_dim_sc %>%
    dplyr::filter(dim_category == "cross-mapped dimorphic") %>%
    dplyr::distinct(target_type) %>% dplyr::pull(target_type)
  mc_cross_sc <- malecns_dim_sc %>%
    dplyr::filter(dim_category == "cross-mapped dimorphic") %>%
    dplyr::distinct(target_type) %>% dplyr::pull(target_type)

  # Types with influence data in each dataset
  banc_has_data <- unique(banc_agg_sc$target_type)
  mc_has_data   <- unique(malecns_agg_sc$target_type)

  # Shared: present in both
  shared_types <- intersect(
    intersect(union(banc_cross_sc, mc_cross_sc), banc_has_data),
    mc_has_data
  )
  # Discrepancies: in one but not both
  banc_only_types <- setdiff(intersect(banc_cross_sc, banc_has_data), mc_has_data)
  mc_only_types   <- setdiff(intersect(mc_cross_sc, mc_has_data), banc_has_data)

  if (length(banc_only_types) > 0 || length(mc_only_types) > 0) {
    discrepancy_df <- dplyr::bind_rows(
      if (length(banc_only_types) > 0) {
        banc_dim_sc %>%
          dplyr::filter(target_type %in% banc_only_types) %>%
          dplyr::transmute(target_type, dataset = "BANC_only",
                           neuron_id = root_id, super_class, dim_category)
      },
      if (length(mc_only_types) > 0) {
        malecns_dim_sc %>%
          dplyr::filter(target_type %in% mc_only_types) %>%
          dplyr::transmute(target_type, dataset = "maleCNS_only",
                           neuron_id = malecns_09_id, super_class, dim_category)
      }
    )
    disc_file <- file.path(output_dir,
      sprintf("dimorphic_type_discrepancies_%s.csv", scg$file_tag))
    readr::write_csv(discrepancy_df, disc_file)
    message(sprintf("  Discrepancies: %d BANC-only types, %d maleCNS-only types → %s",
                    length(banc_only_types), length(mc_only_types), basename(disc_file)))
  }

  if (length(shared_types) == 0) {
    message("  No shared cross-mapped types; skipping")
    next
  }

  # Build matrices with shared types only
  banc_mat <- build_matrix(banc_agg_sc, shared_types)
  mc_mat   <- build_matrix(malecns_agg_sc, shared_types)

  # Rows must have data in BOTH datasets for a meaningful comparison
  all_rows <- intersect(rownames(banc_mat), rownames(mc_mat))
  all_rows <- intersect(all_rows, all_sensory_groups)
  banc_mat <- pad_matrix(banc_mat, all_rows, shared_types)
  mc_mat   <- pad_matrix(mc_mat, all_rows, shared_types)

  # Row thresholding: keep if max in EITHER dataset passes threshold
  banc_row_max <- apply(banc_mat, 1, max, na.rm = TRUE)
  mc_row_max   <- apply(mc_mat, 1, max, na.rm = TRUE)
  combined_maxes <- pmax(banc_row_max, mc_row_max)

  # Drop all-zero columns
  keep_cols <- (colSums(banc_mat, na.rm = TRUE) > 0) |
               (colSums(mc_mat, na.rm = TRUE) > 0)
  banc_mat <- banc_mat[, keep_cols, drop = FALSE]
  mc_mat   <- mc_mat[, keep_cols, drop = FALSE]

  # --- VNC intrinsic: select most divergent columns, two thresholds ---
  if (scg$sc_filter == "ventral_nerve_cord_intrinsic") {
    message("  VNC intrinsic: selecting most divergent columns...")

    # Column divergence = mean absolute difference in column-normalised profiles
    banc_cn <- colwise_minmax_norm(banc_mat)
    mc_cn   <- colwise_minmax_norm(mc_mat)
    col_divergence <- colMeans(abs(banc_cn - mc_cn), na.rm = TRUE)

    for (pctl_thresh in 0.75) {
      pctl_label <- sprintf("p%02d", as.integer(pctl_thresh * 100))

      # Row filter
      row_thresh <- quantile(combined_maxes, pctl_thresh, na.rm = TRUE)
      keep_r <- names(combined_maxes)[combined_maxes >= row_thresh]
      if (length(keep_r) < 3) {
        message(sprintf("  %s: fewer than 3 rows; skipping", pctl_label))
        next
      }

      # Column filter: top 50% most divergent columns
      col_thresh <- quantile(col_divergence, 0.50, na.rm = TRUE)
      keep_c <- names(col_divergence)[col_divergence >= col_thresh]
      if (length(keep_c) < 3) keep_c <- names(sort(col_divergence, decreasing = TRUE))[1:min(20, length(col_divergence))]

      b_sub <- banc_mat[keep_r, keep_c, drop = FALSE]
      m_sub <- mc_mat[keep_r, keep_c, drop = FALSE]

      # Order columns
      c_ord <- cluster_cols(b_sub)
      b_sub <- b_sub[, c_ord, drop = FALSE]
      m_sub <- m_sub[, c_ord, drop = FALSE]

      # Order rows
      ro <- order_rows_by_bp(b_sub)
      b_sub <- b_sub[ro$order, , drop = FALSE]
      m_sub <- m_sub[ro$order, , drop = FALSE]

      message(sprintf("  %s: %d rows × %d cols (most divergent)", pctl_label,
                      nrow(b_sub), ncol(b_sub)))

      draw_condensed_paired(b_sub, m_sub, ro$split, ro$body_parts,
                             scg$label, sprintf("%s most_divergent", pctl_label),
                             output_dir)
    }

  } else {
    # --- Ascending / Descending: single 75th percentile plot ---
    threshold_75 <- quantile(combined_maxes, 0.75, na.rm = TRUE)
    keep_rows <- names(combined_maxes)[combined_maxes >= threshold_75]
    if (length(keep_rows) < 3) {
      message("  Fewer than 3 rows above 75th pctl; skipping")
      next
    }

    banc_sub <- banc_mat[keep_rows, , drop = FALSE]
    mc_sub   <- mc_mat[keep_rows, , drop = FALSE]

    # Order columns and rows
    col_order <- cluster_cols(banc_sub)
    banc_sub <- banc_sub[, col_order, drop = FALSE]
    mc_sub   <- mc_sub[, col_order, drop = FALSE]

    ro <- order_rows_by_bp(banc_sub)
    banc_sub <- banc_sub[ro$order, , drop = FALSE]
    mc_sub   <- mc_sub[ro$order, , drop = FALSE]

    message(sprintf("  75th pctl: %d rows × %d cols", nrow(banc_sub), ncol(banc_sub)))

    draw_condensed_paired(banc_sub, mc_sub, ro$split, ro$body_parts,
                           scg$label, "p75", output_dir)
  }
}

message("=== Done ===")
