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

template_2d_vnc  <- nat.ggplot::ggplot2_neuron_path(
  rgl::as.mesh3d(jrc_surf), rotation_matrix = jrc_vnc_mat)
template_2d_side <- nat.ggplot::ggplot2_neuron_path(
  rgl::as.mesh3d(jrc_surf), rotation_matrix = jrc_vnc_side_mat)

# Weighted-density panel: one facet per dataset × {pre, post}, weighted by DimScore.
plot_weighted_density <- function(df, rotation_matrix, template_2d) {
  facets <- expand.grid(role = levels(df$role),
                         dataset = levels(df$dataset),
                         stringsAsFactors = FALSE)
  template_all <- do.call(rbind, lapply(seq_len(nrow(facets)), function(i) {
    d <- template_2d
    d$role    <- factor(facets$role[i],    levels = levels(df$role))
    d$dataset <- factor(facets$dataset[i], levels = levels(df$dataset))
    d$group <- d$group + (i - 1) * max(template_2d$group, na.rm = TRUE)
    d
  }))

  proj <- project_points(df[, c("X", "Y", "Z")], rotation_matrix)
  pts <- data.frame(px = proj$X, py = proj$Y,
                     DimScore = df$DimScore,
                     dataset  = df$dataset,
                     role     = df$role)

  # stat_density_2d_filled silently DROPS the `weight` aesthetic (its MASS::kde2d
  # backend is unweighted), so passing weight = DimScore would just give a raw
  # synapse-count density. To get a genuine DimScore-weighted density we instead
  # resample points within each facet with probability proportional to DimScore,
  # then estimate an ordinary (unweighted) KDE on that resample.
  pts <- do.call(rbind, lapply(
    split(pts, list(pts$dataset, pts$role), drop = TRUE),
    function(d) {
      if (nrow(d) < 10) return(d)
      w <- pmax(d$DimScore, 0)
      if (sum(w) <= 0) return(d)
      d[sample.int(nrow(d), nrow(d), replace = TRUE, prob = w), , drop = FALSE]
    }))

  ggplot() +
    geom_polygon(data = template_all,
                 aes(x = X, y = Y, group = group),
                 fill = "grey80", colour = NA, alpha = 0.1) +
    stat_density_2d_filled(data = pts,
                            aes(x = px, y = py),
                            h = c(5, 5), bins = 20, n = 200,
                            contour_var = "ndensity",
                            breaks = seq(0.001, 1, length.out = 20),
                            alpha = 0.85) +
    scale_fill_viridis_d(option = "magma",
                          name = "DimScore-\nweighted density",
                          guide = guide_legend(override.aes = list(alpha = 1))) +
    facet_grid(dataset ~ role) +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(strip.text = element_text(size = 11, face = "bold"),
          panel.spacing = unit(0.5, "cm"),
          legend.position = "right")
}

# Root scatter panel: one facet per dataset, points sized + coloured by DimScore.
plot_roots <- function(df, rotation_matrix, template_2d) {
  facets <- data.frame(dataset = levels(df$dataset), stringsAsFactors = FALSE)
  template_all <- do.call(rbind, lapply(seq_len(nrow(facets)), function(i) {
    d <- template_2d
    d$dataset <- factor(facets$dataset[i], levels = levels(df$dataset))
    d$group <- d$group + (i - 1) * max(template_2d$group, na.rm = TRUE)
    d
  }))

  proj <- project_points(df[, c("X", "Y", "Z")], rotation_matrix)
  pts <- data.frame(px = proj$X, py = proj$Y,
                     DimScore = df$DimScore,
                     dataset  = df$dataset)
  # Plot lowest-scoring roots first so the top scorers sit on top.
  pts <- pts[order(pts$DimScore), , drop = FALSE]

  ggplot() +
    geom_polygon(data = template_all,
                 aes(x = X, y = Y, group = group),
                 fill = "grey90", colour = "grey70",
                 linewidth = 0.3, alpha = 0.4) +
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

# ---- 11. Plots -----------------------------------------------------------

message("\n=== Plotting ===")
for (view in c("vnc", "vnc_side")) {
  rot_mat <- if (view == "vnc") jrc_vnc_mat else jrc_vnc_side_mat
  tmpl    <- if (view == "vnc") template_2d_vnc else template_2d_side

  if (nrow(syn_all) > 0) {
    fname <- sprintf("dimscore_synapses_weighted_%s.pdf", view)
    message(sprintf("  Plotting %s ...", fname))
    g <- plot_weighted_density(syn_all, rot_mat, tmpl)
    ggsave(file.path(output_dir, fname), g,
           width = 12, height = 12, dpi = 300, bg = "white")
  }

  if (nrow(roots_all) > 0) {
    fname <- sprintf("dimscore_roots_%s.pdf", view)
    message(sprintf("  Plotting %s ...", fname))
    g <- plot_roots(roots_all, rot_mat, tmpl)
    ggsave(file.path(output_dir, fname), g,
           width = 12, height = 8, dpi = 300, bg = "white")
  }
}

message(sprintf("\n### dim-score density plots complete [%s] ###",
                format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))
