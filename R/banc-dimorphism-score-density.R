#' banc-dimorphism-score-density.R — VNC dimorphism-score weighted density maps
#'
#' Visualises continuous dimorphism scores (DimScore, computed by
#' murthylab/banc_connectivity_analysis) across the VNC in JRCVNC2018U
#' space, in two complementary projections:
#'   1. Synapse maps (pre + post) — 2D KDEs weighted by DimScore so that
#'      high-scoring neurons dominate the density.
#'   2. Root-position maps — soma anchor coordinates plotted as points
#'      coloured/sized by DimScore.
#'
#' Two scoring directions are visualised side-by-side:
#'   * banc_vs_manc — score per BANC root (v626 IDs in source CSV).
#'   * manc_vs_banc — score per MANC body ID.
#'
#' The script is the continuous-score complement to banc-dimorphic-density.R
#' (which uses the categorical sexually_dimorphic SeaTable column).
#'
#' @section Reads:
#'   * data/dimorphism_scores_banc_vs_manc.csv    — BANC root_626 IDs + DimScore
#'   * data/dimorphism_scores_manc_vs_banc.csv    — MANC manc_121_id + DimScore
#'   * banc.meta via startup.R                    — for super_class, root_888, root_position_nm
#'   * GCS banc_888_synapses_v2_enriched.parquet  — BANC synapses (nm)
#'   * GCS manc_121_meta.feather                  — MANC neuron metadata
#'   * GCS manc_121_synapses.parquet              — MANC synapses (µm)
#'   * neuPrint (via malevnc::manc_neuprint_meta) — MANC rootLocation (voxel × 8nm)
#'
#' @section Writes:
#'   * figs/figure_dimorphic/links/dim_scores/<prefix>_<view>.pdf
#'   * figs/figure_dimorphic/links/dim_scores/dimscore_synapses_weighted.txt (plain-language legend)
#'   * figs/figure_dimorphic/links/dim_scores/dimscore_*_jrcvnc2018u.feather (cached point tables)
#'
#' @section Paper:
#'   * BANC paper 2 (Murthy lab) dimorphism scores; doc id
#'     1MUOX8YmrFuuWjsmmGci5Kq9HxQTjkg0KoB2Pt4m_6Xo.
#'
#' @section Notes:
#'   * MANC root positions come from neuPrint via malevnc; if the call
#'     fails (auth / network) the MANC root panel is skipped with a warning.
#'   * BANC IDs in the source CSV are v626; we join through banc.meta to
#'     get the v888 root ID used by the v2-enriched synapse parquet.
#'   * First run downloads two large parquets (BANC ~10 GiB, MANC ~3 GiB)
#'     into the OS temp dir via gsutil; cached for the session.
#'
#' @section Reproduce: BANC_NCORES=1 Rscript R/banc-dimorphism-score-density.R

.this_script <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) normalizePath(f) else stop("Cannot determine script path")
})
.this_dir  <- dirname(.this_script)
.this_repo <- normalizePath(file.path(.this_dir, ".."))

# Deliberately do NOT source startup.R: the dim-score panels only need
# banc.meta (for super_class + root_position_nm + root_626 join), not the
# edgelist or maleCNS bits. We load banc.meta inline from the same GCS
# bucket startup.R uses, so behaviour stays consistent.

suppressMessages({
  library(bancr)
  library(ggplot2)
  library(dplyr)
  library(arrow)
  library(readr)
  library(nat)
  library(nat.flybrains)
  library(nat.templatebrains)
  library(nat.jrcbrains)
  library(nat.ggplot)
  library(rgl)
})

banc.version <- Sys.getenv("BANC_VERSION", unset = "banc_888")
gcs.bucket   <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data"
cache_dir    <- file.path(.this_repo, "data", "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Environment / data-source resolution --------------------------------
# On O2 the large inputs already live on disk under $DATA/connectomes; read
# them in place (NEVER copy multi-GiB parquets into data/cache). Off O2 (e.g.
# a laptop) the same files are pulled from GCS and cached to data/cache. This
# is the only place machine-specific paths live.
on_o2  <- dir.exists("/n/data1/hms/neurobio/wilson")
o2_data <- "/n/data1/hms/neurobio/wilson"                       # == $DATA
o2_connectomes <- file.path(o2_data, "connectomes")

# Per-input registry: logical name -> list(o2 = <local path>, gcs = <gs:// path>).
.inputs <- list(
  banc_meta = list(
    o2  = file.path(o2_connectomes, "banc", banc.version,
                    sprintf("%s_meta.feather", banc.version)),
    gcs = sprintf("%s/%s/%s_meta.feather", gcs.bucket, banc.version, banc.version)),
  banc_synapses = list(
    o2  = file.path(o2_connectomes, "banc", banc.version,
                    sprintf("%s_synapses_v2_enriched.parquet", banc.version)),
    gcs = sprintf("%s/%s/%s_synapses_v2_enriched.parquet",
                  gcs.bucket, banc.version, banc.version)),
  manc_meta = list(
    o2  = file.path(o2_connectomes, "manc", "manc_121_meta.feather"),
    gcs = sprintf("%s/manc_121/manc_121_meta.feather", gcs.bucket)),
  manc_synapses = list(
    o2  = file.path(o2_connectomes, "manc", "manc_121_synapses.parquet"),
    gcs = sprintf("%s/manc_121/manc_121_synapses.parquet", gcs.bucket))
)

# Per (dataset x role) synapse cap applied BEFORE the expensive coordinate
# transform, so draft maps render quickly. DimScore weights are preserved, so
# the weighted density is an unbiased subsample. Set to Inf / 0 to disable.
max_syn <- suppressWarnings(as.numeric(Sys.getenv("BANC_DIM_MAXSYN", unset = "200000")))
if (is.na(max_syn) || max_syn <= 0) max_syn <- Inf

message("### dimorphism-score weighted density maps (VNC) ###")
t_start <- Sys.time()

tryCatch(nat.jrcbrains::register_saalfeldlab_registrations(),
         error = function(e) NULL)

# Persistent cache for the large synapse parquets so reruns don't pay
# the multi-GiB download cost. Lives under data/cache/ (gitignored).
syn_cache_dir <- cache_dir

output_dir <- file.path(.this_repo, "figs", "figure_dimorphic", "links", "dim_scores")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cached point tables (written at the end of the data stage). With
# BANC_DIM_REPLOT=1 we skip the multi-GiB data pull + coordinate transforms and
# re-plot straight from these — seconds instead of ~25 min. Useful when only the
# plotting code changed.
cache_syn_path   <- file.path(output_dir, "dimscore_synapses_jrcvnc2018u.feather")
cache_roots_path <- file.path(output_dir, "dimscore_roots_jrcvnc2018u.feather")
replot_only <- Sys.getenv("BANC_DIM_REPLOT", "0") == "1" &&
  file.exists(cache_syn_path) && file.exists(cache_roots_path)

# Template surface — needed by both the data stage (point-in-VNC test) and the
# plot stage (mesh outline), so define it unconditionally.
jrc_surf <- nat.flybrains::JRCVNC2018U.surf

# ---- 1. Helpers -----------------------------------------------------------

gcs_cache <- function(gcs_path, cache_dir = syn_cache_dir) {
  local_file <- file.path(cache_dir, basename(gcs_path))
  if (!file.exists(local_file)) {
    message(sprintf("  Caching %s ...", basename(gcs_path)))
    out <- system2("gsutil", c("cp", gcs_path, local_file),
                   stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
      unlink(local_file)
      stop(sprintf("gsutil cp failed for %s:\n%s", basename(gcs_path),
                   paste(out, collapse = "\n")))
    }
  } else {
    message(sprintf("  Using cached %s", basename(local_file)))
  }
  local_file
}

# Resolve a logical input name to a readable local path. On O2 returns the
# in-place file under $DATA/connectomes (no copy); otherwise downloads + caches
# the GCS object. Falls back to GCS if the expected O2 file is missing.
resolve_input <- function(name) {
  spec <- .inputs[[name]]
  if (is.null(spec)) stop("Unknown input: ", name)
  if (on_o2 && file.exists(spec$o2)) {
    message(sprintf("  [O2 local] %s", spec$o2))
    return(spec$o2)
  }
  if (on_o2)
    message(sprintf("  [O2] %s missing; falling back to GCS", basename(spec$o2)))
  gcs_cache(spec$gcs)
}

# Collect a filtered arrow synapse query, capping at `cap` rows BEFORE pulling
# into R (the full result can be tens of millions of rows). Subsampling is done
# server-side via deterministic modulo on the integer-cast string id column, so
# memory stays bounded and the sample is spatially unbiased (uniform over
# synapses => DimScore-weighted density is unbiased). Requires `key_col` to be a
# numeric-as-string id (arrow casts it to int64).
collect_capped <- function(query, key_col, label, cap = max_syn) {
  if (is.infinite(cap)) return(dplyr::collect(query))
  n <- query %>% dplyr::count() %>% dplyr::collect() %>% dplyr::pull(n)
  if (n <= cap) {
    message(sprintf("    %s: %d rows (<= cap)", label, n))
    return(dplyr::collect(query))
  }
  k <- as.integer(ceiling(n / cap))
  message(sprintf("    %s: %d rows -> ~%d (every %dth: int64(%s) %%%% %d == 0)",
                  label, n, n %/% k, k, key_col, k))
  mod_expr <- rlang::expr(cast(!!rlang::sym(key_col), int64()) %% !!k == 0L)
  query %>% dplyr::filter(!!mod_expr) %>% dplyr::collect()
}

# `open_dataset` is the predicate-pushdown path; fall back to `read_parquet`
# on HPC mounts where the URI dispatcher misbehaves.
open_parquet <- function(path) {
  tryCatch(arrow::open_dataset(path, format = "parquet"),
           error = function(e) {
             if (grepl("filesystem", e$message, ignore.case = TRUE)) {
               message("  open_dataset failed; using read_parquet fallback")
               arrow::read_parquet(path, as_data_frame = FALSE)
             } else stop(e)
           })
}

project_points <- function(xyz, rotation_matrix) {
  projected <- t(rotation_matrix[1:3, 1:3] %*% t(as.matrix(xyz)))
  data.frame(X = projected[, 1], Y = projected[, 2], Z = projected[, 3])
}

# Parse a "x,y,z" or "x, y, z" string column into a numeric matrix.
parse_xyz_string <- function(s) {
  m <- matrix(NA_real_, nrow = length(s), ncol = 3)
  ok <- !is.na(s) & nchar(s) > 0
  if (!any(ok)) return(m)
  parts <- strsplit(gsub("\\s+", "", s[ok]), ",", fixed = TRUE)
  lens  <- lengths(parts)
  good  <- which(ok)[lens == 3]
  if (length(good)) {
    vals <- suppressWarnings(as.numeric(unlist(parts[lens == 3])))
    m[good, ] <- matrix(vals, ncol = 3, byrow = TRUE)
  }
  m
}

# VNC dorsal: 270° rotation (long axis vertical)
jrc_vnc_mat <- matrix(c( 0,-1, 0, 0,
                         1, 0, 0, 0,
                         0, 0, 1, 0,
                         0, 0, 0, 1), 4, 4, byrow = TRUE)

# VNC lateral: Z->display plane, 270° rotation
jrc_vnc_side_mat <- matrix(c( 0,-1, 0, 0,
                              0, 0, 1, 0,
                              1, 0, 0, 0,
                              0, 0, 0, 1), 4, 4, byrow = TRUE)

# ---- 2. Coordinate transforms ---------------------------------------------

banc_to_jrcvnc2018u <- function(xyz) {
  xyz_jrc <- bancr::banc_to_JRC2018F(xyz, region = "vnc",
                                       method = "tpsreg", banc.units = "nm",
                                       inverse = FALSE)
  nat.templatebrains::xform_brain(xyz_jrc, sample = "JRCVNC2018F",
                                   reference = "JRCVNC2018U")
}

manc_to_jrcvnc2018u <- function(xyz_um) {
  xyz_jrc <- nat.templatebrains::xform_brain(xyz_um, sample = "MANC",
                                              reference = "JRCVNC2018F")
  nat.templatebrains::xform_brain(xyz_jrc, sample = "JRCVNC2018F",
                                   reference = "JRCVNC2018U")
}

# ---- 3a. Load BANC metadata ----------------------------------------------

if (replot_only) {
  message("\n### Replot mode (BANC_DIM_REPLOT=1): loading cached point tables ###")
  syn_all   <- arrow::read_feather(cache_syn_path)
  roots_all <- arrow::read_feather(cache_roots_path)
} else {

message("\n=== Loading BANC metadata ===")
banc_meta_file <- resolve_input("banc_meta")
banc.meta <- arrow::read_feather(banc_meta_file) %>%
  dplyr::mutate(across(where(bit64::is.integer64), as.character))
.id_col <- grep("^banc_[0-9]+_id$", colnames(banc.meta), value = TRUE)
if (length(.id_col) == 1 && .id_col != "root_id") {
  banc.meta$root_id <- as.character(banc.meta[[.id_col]])
}
banc.meta$root_id  <- as.character(banc.meta$root_id)
banc.meta$root_626 <- as.character(banc.meta$root_626)
message(sprintf("  banc.meta: %d rows, %d columns", nrow(banc.meta), ncol(banc.meta)))

# ---- 3. Load dimorphism-score tables -------------------------------------

message("\n=== Loading dimorphism-score CSVs ===")
banc_scores <- readr::read_csv(
  file.path(.this_repo, "data", "dimorphism_scores_banc_vs_manc.csv"),
  col_types = "cd", show_col_types = FALSE) %>%
  dplyr::rename(root_626 = `BANC ID`, DimScore = DimScore) %>%
  dplyr::mutate(root_626 = as.character(root_626)) %>%
  dplyr::distinct(root_626, .keep_all = TRUE)
message(sprintf("  banc_vs_manc: %d rows, DimScore range [%.2f, %.2f]",
                nrow(banc_scores), min(banc_scores$DimScore),
                max(banc_scores$DimScore)))

manc_scores <- readr::read_csv(
  file.path(.this_repo, "data", "dimorphism_scores_manc_vs_banc.csv"),
  col_types = "cd", show_col_types = FALSE) %>%
  dplyr::rename(manc_id = `MANC ID`, DimScore = DimScore) %>%
  dplyr::mutate(manc_id = as.character(manc_id)) %>%
  dplyr::distinct(manc_id, .keep_all = TRUE)
message(sprintf("  manc_vs_banc: %d rows, DimScore range [%.2f, %.2f]",
                nrow(manc_scores), min(manc_scores$DimScore),
                max(manc_scores$DimScore)))

# ---- 4. Join scores to per-neuron metadata --------------------------------

message("\n=== Joining scores to neuron metadata ===")

bc <- banc_scores %>%
  dplyr::inner_join(
    banc.meta %>%
      dplyr::select(root_626, root_id, super_class, cell_type, root_position_nm) %>%
      dplyr::distinct(root_626, .keep_all = TRUE),
    by = "root_626")
message(sprintf("  BANC scores joined to banc.meta: %d / %d",
                nrow(bc), nrow(banc_scores)))

# MANC metadata (manc_121): local on O2, else GCS.
manc_meta_local <- resolve_input("manc_meta")
manc.meta <- arrow::read_feather(manc_meta_local) %>%
  dplyr::mutate(manc_id = as.character(manc_121_id))

mc <- manc_scores %>%
  dplyr::inner_join(
    manc.meta %>%
      dplyr::select(manc_id, super_class, cell_type) %>%
      dplyr::distinct(manc_id, .keep_all = TRUE),
    by = "manc_id")
message(sprintf("  MANC scores joined to manc_121 meta: %d / %d",
                nrow(mc), nrow(manc_scores)))

# ---- 5. BANC root positions ----------------------------------------------

message("\n=== BANC root positions ===")
banc_roots_nm <- parse_xyz_string(bc$root_position_nm)
banc_root_valid <- complete.cases(banc_roots_nm)
message(sprintf("  %d / %d BANC neurons have a parseable root_position_nm",
                sum(banc_root_valid), nrow(bc)))

banc_roots_u <- matrix(NA_real_, nrow = nrow(bc), ncol = 3)
if (any(banc_root_valid)) {
  banc_roots_u[banc_root_valid, ] <-
    banc_to_jrcvnc2018u(banc_roots_nm[banc_root_valid, , drop = FALSE])
}

# ---- 6. MANC root positions (neuPrint via malevnc) -----------------------

message("\n=== MANC root positions (neuPrint) ===")
manc_roots_df <- tryCatch({
  suppressMessages({
    library(malevnc)
    nm <- malevnc::manc_neuprint_meta(mc$manc_id)
  })
  bid_col <- intersect(c("bodyid", "bodyId", "manc_id"), names(nm))[1]
  if (is.na(bid_col))
    stop("neuPrint meta lacks bodyid column")
  nm$.bid <- as.character(nm[[bid_col]])
  # Per neuron prefer somaLocation, then tosomaLocation, then rootLocation.
  # somaLocation is the true cell body; the others are skeleton roots.
  pick_loc <- function(...) {
    cols <- c(...)
    cols <- cols[cols %in% names(nm)]
    if (!length(cols)) return(rep(NA_character_, nrow(nm)))
    out <- rep(NA_character_, nrow(nm))
    for (cc in cols) {
      v <- nm[[cc]]; v[is.na(v) | nchar(v) == 0] <- NA_character_
      out[is.na(out)] <- v[is.na(out)]
    }
    out
  }
  nm$.loc <- pick_loc("somaLocation", "tosomaLocation", "rootLocation")
  message(sprintf(
    "  neuPrint location coverage — soma: %d, tosoma: %d, root: %d, any: %d / %d",
    sum(!is.na(nm$somaLocation) & nchar(nm$somaLocation) > 0),
    sum(!is.na(nm$tosomaLocation) & nchar(nm$tosomaLocation) > 0),
    sum(!is.na(nm$rootLocation) & nchar(nm$rootLocation) > 0),
    sum(!is.na(nm$.loc)), nrow(nm)))
  dplyr::left_join(data.frame(manc_id = mc$manc_id, stringsAsFactors = FALSE),
                   nm[, c(".bid", ".loc")], by = c("manc_id" = ".bid")) %>%
    dplyr::rename(rootLocation = .loc)
}, error = function(e) {
  warning("MANC neuPrint fetch failed (", e$message,
          "); skipping MANC root panel.")
  NULL
})

manc_roots_u <- matrix(NA_real_, nrow = nrow(mc), ncol = 3)
manc_root_valid <- rep(FALSE, nrow(mc))
if (!is.null(manc_roots_df)) {
  manc_roots_vox <- parse_xyz_string(manc_roots_df$rootLocation)
  manc_root_valid <- complete.cases(manc_roots_vox)
  message(sprintf("  %d / %d MANC neurons have a neuPrint rootLocation",
                  sum(manc_root_valid), nrow(mc)))
  if (any(manc_root_valid)) {
    # neuPrint voxel coords at 8 nm/voxel → MANC native µm (×0.008).
    manc_roots_um <- manc_roots_vox[manc_root_valid, , drop = FALSE] * 0.008
    manc_roots_u[manc_root_valid, ] <- manc_to_jrcvnc2018u(manc_roots_um)
  }
}

# ---- 7. BANC synapses ----------------------------------------------------

message("\n=== BANC synapses ===")
set.seed(42)  # reproducible subsampling
banc_syn_file <- resolve_input("banc_synapses")

# Drop neurons that did not survive the join (no v888 root_id available).
bc_with_888 <- bc %>%
  dplyr::filter(!is.na(root_id), nchar(root_id) > 0) %>%
  dplyr::distinct(root_id, .keep_all = TRUE)
banc_888_ids <- bc_with_888$root_id
message(sprintf("  %d BANC neurons with v888 root_id (of %d scored)",
                length(banc_888_ids), nrow(bc)))

banc_syn_ds <- open_parquet(banc_syn_file)
banc_pre <- collect_capped(
  banc_syn_ds %>%
    dplyr::select(id, pre_root_id, X, Y, Z) %>%
    dplyr::filter(pre_root_id %in% banc_888_ids),
  key_col = "id", label = "BANC pre") %>%
  dplyr::inner_join(bc_with_888 %>% dplyr::select(root_id, DimScore, super_class),
                    by = c("pre_root_id" = "root_id"))
gc()
message(sprintf("  presynapses from scored BANC neurons (collected): %d", nrow(banc_pre)))

banc_post <- collect_capped(
  banc_syn_ds %>%
    dplyr::select(id, post_root_id, X, Y, Z) %>%
    dplyr::filter(post_root_id %in% banc_888_ids),
  key_col = "id", label = "BANC post") %>%
  dplyr::inner_join(bc_with_888 %>% dplyr::select(root_id, DimScore, super_class),
                    by = c("post_root_id" = "root_id"))
gc()
message(sprintf("  postsynapses to scored BANC neurons (collected): %d", nrow(banc_post)))

# Transform + restrict to JRCVNC2018U volume.
jrc_surf <- nat.flybrains::JRCVNC2018U.surf
banc_pre_in <- (function() {
  if (!nrow(banc_pre)) return(NULL)
  u <- banc_to_jrcvnc2018u(nat::xyzmatrix(banc_pre[, c("X", "Y", "Z")]))
  inside <- nat::pointsinside(u, jrc_surf); inside[is.na(inside)] <- FALSE
  message(sprintf("  BANC presynapses within JRCVNC2018U: %d", sum(inside)))
  data.frame(X = u[inside, 1], Y = u[inside, 2], Z = u[inside, 3],
             DimScore = banc_pre$DimScore[inside],
             super_class = banc_pre$super_class[inside],
             dataset = "BANC", role = "pre")
})()
banc_post_in <- (function() {
  if (!nrow(banc_post)) return(NULL)
  u <- banc_to_jrcvnc2018u(nat::xyzmatrix(banc_post[, c("X", "Y", "Z")]))
  inside <- nat::pointsinside(u, jrc_surf); inside[is.na(inside)] <- FALSE
  message(sprintf("  BANC postsynapses within JRCVNC2018U: %d", sum(inside)))
  data.frame(X = u[inside, 1], Y = u[inside, 2], Z = u[inside, 3],
             DimScore = banc_post$DimScore[inside],
             super_class = banc_post$super_class[inside],
             dataset = "BANC", role = "post")
})()
rm(banc_pre, banc_post); gc()

# ---- 8. MANC synapses ----------------------------------------------------

message("\n=== MANC synapses ===")
manc_syn_file <- resolve_input("manc_synapses")
manc_syn_ds <- open_parquet(manc_syn_file)
manc_ids_str <- as.character(mc$manc_id)

manc_pre <- collect_capped(
  manc_syn_ds %>%
    dplyr::filter(prepost == 0L, pre %in% manc_ids_str) %>%
    dplyr::select(connector_id, pre, x, y, z),
  key_col = "connector_id", label = "MANC pre") %>%
  dplyr::inner_join(mc %>% dplyr::select(manc_id, DimScore, super_class),
                    by = c("pre" = "manc_id"))
gc()
message(sprintf("  presynapses from scored MANC neurons (collected): %d", nrow(manc_pre)))

manc_post <- collect_capped(
  manc_syn_ds %>%
    dplyr::filter(prepost == 1L, post %in% manc_ids_str) %>%
    dplyr::select(connector_id, post, x, y, z),
  key_col = "connector_id", label = "MANC post") %>%
  dplyr::inner_join(mc %>% dplyr::select(manc_id, DimScore, super_class),
                    by = c("post" = "manc_id"))
gc()
message(sprintf("  postsynapses to scored MANC neurons (collected): %d", nrow(manc_post)))

manc_pre_in <- (function() {
  if (!nrow(manc_pre)) return(NULL)
  u <- manc_to_jrcvnc2018u(nat::xyzmatrix(manc_pre[, c("x", "y", "z")]))
  inside <- nat::pointsinside(u, jrc_surf); inside[is.na(inside)] <- FALSE
  message(sprintf("  MANC presynapses within JRCVNC2018U: %d", sum(inside)))
  data.frame(X = u[inside, 1], Y = u[inside, 2], Z = u[inside, 3],
             DimScore = manc_pre$DimScore[inside],
             super_class = manc_pre$super_class[inside],
             dataset = "MANC", role = "pre")
})()
manc_post_in <- (function() {
  if (!nrow(manc_post)) return(NULL)
  u <- manc_to_jrcvnc2018u(nat::xyzmatrix(manc_post[, c("x", "y", "z")]))
  inside <- nat::pointsinside(u, jrc_surf); inside[is.na(inside)] <- FALSE
  message(sprintf("  MANC postsynapses within JRCVNC2018U: %d", sum(inside)))
  data.frame(X = u[inside, 1], Y = u[inside, 2], Z = u[inside, 3],
             DimScore = manc_post$DimScore[inside],
             super_class = manc_post$super_class[inside],
             dataset = "MANC", role = "post")
})()
rm(manc_pre, manc_post); gc()

# ---- 9. Combine + persist point tables -----------------------------------

message("\n=== Combining tables ===")
syn_all <- dplyr::bind_rows(banc_pre_in, banc_post_in,
                             manc_pre_in, manc_post_in)
syn_all$dataset <- factor(syn_all$dataset, levels = c("BANC", "MANC"))
syn_all$role    <- factor(syn_all$role,    levels = c("pre", "post"))
arrow::write_feather(syn_all,
  file.path(output_dir, "dimscore_synapses_jrcvnc2018u.feather"))
message(sprintf("  %d total synapse points across both datasets",
                nrow(syn_all)))

roots_banc <- data.frame(
  X = banc_roots_u[banc_root_valid, 1],
  Y = banc_roots_u[banc_root_valid, 2],
  Z = banc_roots_u[banc_root_valid, 3],
  DimScore = bc$DimScore[banc_root_valid],
  super_class = bc$super_class[banc_root_valid],
  dataset = "BANC")
roots_manc <- data.frame(
  X = manc_roots_u[manc_root_valid, 1],
  Y = manc_roots_u[manc_root_valid, 2],
  Z = manc_roots_u[manc_root_valid, 3],
  DimScore = mc$DimScore[manc_root_valid],
  super_class = mc$super_class[manc_root_valid],
  dataset = "MANC")
roots_all <- dplyr::bind_rows(roots_banc, roots_manc)
roots_all$dataset <- factor(roots_all$dataset, levels = c("BANC", "MANC"))
# Restrict roots to inside the JRCVNC2018U template; somas outside the
# VNC (e.g. brain ANs / DNs in BANC) would otherwise dominate the plot
# extent and squash the meaningful range.
roots_inside <- nat::pointsinside(as.matrix(roots_all[, c("X", "Y", "Z")]),
                                   jrc_surf)
roots_inside[is.na(roots_inside)] <- FALSE
roots_all <- roots_all[roots_inside, , drop = FALSE]
arrow::write_feather(roots_all,
  file.path(output_dir, "dimscore_roots_jrcvnc2018u.feather"))
message(sprintf("  %d root points within JRCVNC2018U (BANC + MANC)",
                nrow(roots_all)))

}  # end if(!replot_only)

# Ensure consistent factor levels regardless of branch (read_feather can return
# these columns as plain character).
syn_all$dataset   <- factor(syn_all$dataset, levels = c("BANC", "MANC"))
syn_all$role      <- factor(syn_all$role,    levels = c("pre", "post"))
roots_all$dataset <- factor(roots_all$dataset, levels = c("BANC", "MANC"))

# ---- 10. Plot helpers ----------------------------------------------------

# Projected mesh triangles (one group per face) for the VNC neuropil surface.
# Rendered as a *solid filled silhouette* (fill == stroke colour, hairline width)
# so adjacent triangles merge into one clean shape with no visible wireframe —
# the previous per-triangle stroke read as a wire mesh.
template_2d_vnc  <- nat.ggplot::ggplot2_neuron_path(
  rgl::as.mesh3d(jrc_surf), rotation_matrix = jrc_vnc_mat)
template_2d_side <- nat.ggplot::ggplot2_neuron_path(
  rgl::as.mesh3d(jrc_surf), rotation_matrix = jrc_vnc_side_mat)

# Replicate the projected surface across facets (group offset keeps polygons
# from bleeding between panels). `facet_cols` is a named list of factor columns
# to attach (e.g. list(dataset=..., role=...)).
template_per_facet <- function(template_2d, facets) {
  facets <- as.data.frame(facets, stringsAsFactors = FALSE)
  base_max <- max(template_2d$group, na.rm = TRUE)
  do.call(rbind, lapply(seq_len(nrow(facets)), function(i) {
    d <- template_2d
    for (nm in names(facets)) d[[nm]] <- facets[[nm]][i]
    d$group <- d$group + (i - 1) * base_max
    d
  }))
}

# Clean silhouette layer: fill and stroke share a colour so the triangulation
# disappears, leaving a single smooth neuropil outline.
geom_surface <- function(template_all, fill = "grey92") {
  geom_polygon(data = template_all,
               aes(x = X, y = Y, group = group),
               fill = fill, colour = fill, linewidth = 0.05)
}

# Row-normalised 1D Gaussian smoothing matrix: n bins of width `bw` (µm),
# bandwidth `sigma` (µm). Edges handled by truncation (rows renormalised).
gauss_smooth_mat <- function(n, bw, sigma) {
  idx <- seq_len(n)
  W <- exp(-(outer(idx, idx, "-") * bw)^2 / (2 * sigma^2))
  W / rowSums(W)
}

# Build a 1D bin grid covering `vals` at resolution `bw` (µm).
make_grid <- function(vals, bw) {
  lo <- min(vals); hi <- max(vals)
  n  <- max(1L, as.integer(ceiling((hi - lo) / bw)))
  edges <- lo + (0:n) * bw
  list(n = n, edges = edges, centers = edges[-1] - bw / 2, bw = bw)
}

# Depth max-intensity projection (MIP) of a DimScore-weighted point cloud.
# Points are binned into a 3D grid (display-x, display-y, depth) with DimScore
# as the weight, smoothed with a separable 3D Gaussian, then the MAXIMUM along
# the depth axis is taken for each (x, y) pixel. This surfaces depth-localised
# hotspots (e.g. a single neuropil layer) that a flat through-depth sum would
# bury under thicker, busier regions. Returns a long data.frame (px, py, value)
# with value normalised to [0, 1]; bins below `thr` are set NA (transparent).
mip_density <- function(px, py, depth, w, gx, gy, gz, Wx, Wy, Wz, thr = 0.03) {
  ix <- findInterval(px,    gx$edges, rightmost.closed = TRUE)
  iy <- findInterval(py,    gy$edges, rightmost.closed = TRUE)
  iz <- findInterval(depth, gz$edges, rightmost.closed = TRUE)
  ok <- ix >= 1 & ix <= gx$n & iy >= 1 & iy <= gy$n & iz >= 1 & iz <= gz$n
  A  <- array(0, dim = c(gx$n, gy$n, gz$n))
  if (any(ok)) {
    lin <- ix[ok] + (iy[ok] - 1) * gx$n + (iz[ok] - 1) * gx$n * gy$n
    rs  <- rowsum(pmax(w[ok], 0), group = lin)          # fast weighted accumulate
    A[as.integer(rownames(rs))] <- rs[, 1]
  }
  nx <- gx$n; ny <- gy$n; nz <- gz$n
  # Separable Gaussian smoothing, one axis at a time, via matrix products.
  A <- array(Wx %*% matrix(A, nrow = nx), dim = c(nx, ny, nz))
  A <- aperm(A, c(2, 1, 3))
  A <- array(Wy %*% matrix(A, nrow = ny), dim = c(ny, nx, nz))
  A <- aperm(A, c(2, 1, 3))
  A <- aperm(A, c(3, 1, 2))
  A <- array(Wz %*% matrix(A, nrow = nz), dim = c(nz, nx, ny))
  A <- aperm(A, c(2, 3, 1))
  M <- apply(A, c(1, 2), max)                            # max-intensity projection
  mx <- max(M, na.rm = TRUE)
  if (is.finite(mx) && mx > 0) M <- M / mx               # per-panel ndensity 0..1
  out <- expand.grid(px = gx$centers, py = gy$centers)
  out$value <- as.vector(M)
  out$value[out$value < thr] <- NA_real_
  out
}

# Density panel: one facet per dataset × {pre, post}. Each panel is a DimScore-
# weighted depth max-intensity projection (see mip_density), rendered as a raster
# over the clean neuropil silhouette.
plot_weighted_density <- function(df, rotation_matrix, template_2d,
                                  bw = 2.5, sigma = 5) {
  facets <- expand.grid(role = levels(df$role),
                         dataset = levels(df$dataset),
                         stringsAsFactors = FALSE)
  template_all <- template_per_facet(template_2d, list(
    role    = factor(facets$role,    levels = levels(df$role)),
    dataset = factor(facets$dataset, levels = levels(df$dataset))))

  proj <- project_points(df[, c("X", "Y", "Z")], rotation_matrix)
  tproj <- project_points(template_2d[, c("X", "Y", "Z")], rotation_matrix)

  # Common (px, py) grid across panels (spans points + silhouette so the raster
  # registers with the surface); shared depth grid for the projection axis.
  gx <- make_grid(c(proj$X, tproj$X), bw)
  gy <- make_grid(c(proj$Y, tproj$Y), bw)
  gz <- make_grid(proj$Z, bw)
  Wx <- gauss_smooth_mat(gx$n, bw, sigma)
  Wy <- gauss_smooth_mat(gy$n, bw, sigma)
  Wz <- gauss_smooth_mat(gz$n, bw, sigma)

  pts <- data.frame(px = proj$X, py = proj$Y, depth = proj$Z,
                     DimScore = df$DimScore,
                     dataset = df$dataset, role = df$role)
  rast <- do.call(rbind, lapply(
    split(pts, list(pts$dataset, pts$role), drop = TRUE),
    function(d) {
      if (nrow(d) < 10) return(NULL)
      m <- mip_density(d$px, d$py, d$depth, d$DimScore, gx, gy, gz, Wx, Wy, Wz)
      m <- m[!is.na(m$value), , drop = FALSE]
      if (!nrow(m)) return(NULL)
      m$dataset <- factor(d$dataset[1], levels = levels(df$dataset))
      m$role    <- factor(d$role[1],    levels = levels(df$role))
      m
    }))

  ggplot() +
    geom_surface(template_all, fill = "grey92") +
    geom_raster(data = rast, aes(x = px, y = py, fill = value),
                interpolate = TRUE) +
    scale_fill_viridis_c(option = "magma", na.value = "transparent",
                         name = "DimScore-weighted\ndensity\n(depth max-proj,\nper-panel 0–1)") +
    facet_grid(dataset ~ role) +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(strip.text = element_text(size = 11, face = "bold"),
          panel.spacing = unit(0.5, "cm"),
          legend.position = "right")
}

# Root scatter panel: one facet per dataset, points sized + coloured by DimScore.
plot_roots <- function(df, rotation_matrix, template_2d) {
  template_all <- template_per_facet(template_2d, list(
    dataset = factor(levels(df$dataset), levels = levels(df$dataset))))

  proj <- project_points(df[, c("X", "Y", "Z")], rotation_matrix)
  pts <- data.frame(px = proj$X, py = proj$Y,
                     DimScore = df$DimScore,
                     dataset  = df$dataset)
  # Plot lowest-scoring roots first so the top scorers sit on top.
  pts <- pts[order(pts$DimScore), , drop = FALSE]

  ggplot() +
    geom_surface(template_all, fill = "grey90") +
    geom_point(data = pts,
                aes(x = px, y = py, colour = DimScore, size = DimScore),
                alpha = 0.7) +
    scale_colour_viridis_c(option = "magma", name = "DimScore",
                            limits = range(pts$DimScore, na.rm = TRUE)) +
    scale_size_continuous(range = c(0.4, 3), name = "DimScore") +
    guides(colour = guide_colourbar(), size = "none") +
    facet_wrap(~ dataset, nrow = 1) +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(strip.text = element_text(size = 11, face = "bold"),
          panel.spacing = unit(0.5, "cm"),
          legend.position = "right")
}

# ---- 10b. Figure legend (synapse panel) ----------------------------------
# Plain-language caption written alongside the synapse figure so reviewers can
# read it without the script. Kept here (not a static file) so it always tracks
# the current method.
write_synapse_legend <- function(path) {
  writeLines(c(
"Dimorphism-score-weighted synaptic density across the ventral nerve cord (VNC)",
"==============================================================================",
"",
"Files: dimscore_synapses_weighted_vnc.{png,pdf}        (dorsal / top-down view)",
"       dimscore_synapses_weighted_vnc_side.{png,pdf}   (lateral / side view)",
"",
"What the figure shows",
"---------------------",
"Where in the VNC the synapses of sexually dimorphic neurons are concentrated.",
"Every synapse is weighted by its neuron's dimorphism score (DimScore): a higher",
"score means the neuron differs more between the male and female nerve cord, so",
"bright regions are places rich in strongly dimorphic synapses. Brightness is",
"NOT a raw synapse count.",
"",
"Layout: four panels.",
"  Rows    - dataset: BANC (top) and MANC (bottom).",
"  Columns - synapse side: 'pre' (outputs, where the neuron talks to others)",
"            and 'post' (inputs, where it is talked to).",
"All synapses are warped into a common template nerve cord (JRCVNC2018U) so the",
"two datasets sit in the same anatomical space and can be compared panel to panel.",
"",
"How depth is handled (maximum-intensity projection)",
"---------------------------------------------------",
"The nerve cord is 3D but the figure is 2D. For each pixel we look straight",
"through the full depth of the cord and keep the single brightest value along",
"that line of sight (a maximum-intensity projection, as in confocal microscopy).",
"Nothing is hidden behind a front surface, and a hotspot confined to one layer",
"(e.g. a wing-related region) stays visible instead of being averaged out by the",
"thicker, busier leg regions in front of or behind it.",
"",
"Colour scale",
"------------",
"Magma colour = relative density, scaled independently within each panel from 0",
"(dark / transparent) to 1 (bright). Because each panel is scaled to its own",
"maximum, colours show WHERE the dimorphic synapses sit within a panel, not how",
"the four panels compare in absolute magnitude.",
"",
"Caveats",
"-------",
"- Draft: synapses are sub-sampled per panel for speed (see BANC_DIM_MAXSYN);",
"  the weighting keeps the sample unbiased but final figures should raise the cap.",
"- Grey shape = the template VNC neuropil surface, for spatial reference only.",
"",
"Source: R/banc-dimorphism-score-density.R (murthylab/banc_connectivity_analysis)."),
    con = path)
}

# ---- 11. Plots -----------------------------------------------------------

message("\n=== Plotting ===")
for (view in c("vnc", "vnc_side")) {
  rot_mat <- if (view == "vnc") jrc_vnc_mat else jrc_vnc_side_mat
  tmpl    <- if (view == "vnc") template_2d_vnc else template_2d_side

  if (nrow(syn_all) > 0) {
    stem <- sprintf("dimscore_synapses_weighted_%s", view)
    message(sprintf("  Plotting %s ...", stem))
    g <- plot_weighted_density(syn_all, rot_mat, tmpl)
    ggsave(file.path(output_dir, paste0(stem, ".pdf")), g,
           width = 12, height = 12, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, paste0(stem, ".png")), g,
           width = 12, height = 12, dpi = 200, bg = "white")
    write_synapse_legend(file.path(output_dir, "dimscore_synapses_weighted.txt"))
  }

  if (nrow(roots_all) > 0) {
    stem <- sprintf("dimscore_roots_%s", view)
    message(sprintf("  Plotting %s ...", stem))
    g <- plot_roots(roots_all, rot_mat, tmpl)
    ggsave(file.path(output_dir, paste0(stem, ".pdf")), g,
           width = 12, height = 8, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, paste0(stem, ".png")), g,
           width = 12, height = 8, dpi = 200, bg = "white")
  }
}

message(sprintf("\n### dim-score density plots complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))
