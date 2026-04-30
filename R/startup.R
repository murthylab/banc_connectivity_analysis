###############################################
### banc_connectivity_analysis — startup.R  ###
###############################################
# Minimal startup for BANC Paper 2 analysis scripts.
# Loads BANC + maleCNS metadata and edgelists from GCS.
# Independent of the BANC-project repo.

# ── Configuration ──────────────────────────────────────────
banc.version <- Sys.getenv("BANC_VERSION", unset = "banc_888")
.version_num <- as.integer(sub("^banc_", "", banc.version))

gcs.bucket <- "gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data"
gcs.bucket.legacy <- "gs://brain-and-nerve-cord_exports"

.repo_root <- tryCatch(
  normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE),
  error = function(e) getwd()
)
cache_dir <- file.path(.repo_root, "data", "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ── Libraries ──────────────────────────────────────────────
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(arrow)
  library(reticulate)
  library(bancr)
  library(influencer)
})

# ── GCS helpers ────────────────────────────────────────────
# Read a feather file from GCS, caching locally.
read_feather_gcs <- function(gcs_path, cache = TRUE) {
  fname <- basename(gcs_path)
  local_path <- file.path(cache_dir, fname)
  if (cache && file.exists(local_path)) {
    message("Loading from cache: ", local_path)
    return(arrow::read_feather(local_path))
  }
  message("Downloading from GCS: ", gcs_path)
  tmp <- tempfile(fileext = paste0(".", tools::file_ext(fname)))
  system2("gsutil", c("cp", gcs_path, tmp), stdout = TRUE, stderr = TRUE)
  df <- arrow::read_feather(tmp)
  if (cache) {
    file.copy(tmp, local_path, overwrite = TRUE)
    message("Cached to: ", local_path)
  }
  unlink(tmp)
  df
}

read_parquet_gcs <- function(gcs_path, cache = TRUE) {
  fname <- basename(gcs_path)
  local_path <- file.path(cache_dir, fname)
  if (cache && file.exists(local_path)) {
    message("Loading from cache: ", local_path)
    return(arrow::read_parquet(local_path))
  }
  message("Downloading from GCS: ", gcs_path)
  tmp <- tempfile(fileext = ".parquet")
  system2("gsutil", c("cp", gcs_path, tmp), stdout = TRUE, stderr = TRUE)
  df <- arrow::read_parquet(tmp)
  if (cache) {
    file.copy(tmp, local_path, overwrite = TRUE)
    message("Cached to: ", local_path)
  }
  unlink(tmp)
  df
}

# ── Paper colours ──────────────────────────────────────────
.colours_file <- file.path(.repo_root, "settings", "paper_colours_lacroix.csv")
if (!file.exists(.colours_file)) {
  warning("paper_colours_lacroix.csv not found at ", .colours_file)
}
if (file.exists(.colours_file)) {
  .pcdf <- read.csv(.colours_file, stringsAsFactors = FALSE)
  paper.cols <- setNames(.pcdf$hex, .pcdf$label)
  message("Loaded ", length(paper.cols), " paper colours")
} else {
  paper.cols <- character(0)
  message("WARNING: paper_colours_lacroix.csv not found")
}

# ── BANC metadata ──────────────────────────────────────────
message("Loading BANC metadata...")
.banc_meta_gcs <- sprintf("%s/%s/%s_meta.feather", gcs.bucket, banc.version, banc.version)
banc.meta <- read_feather_gcs(.banc_meta_gcs) %>%
  dplyr::mutate(across(where(is.integer64), as.character))

# Normalise root_id column
.id_col <- grep("^banc_[0-9]+_id$", colnames(banc.meta), value = TRUE)
if (length(.id_col) == 1 && .id_col != "root_id") {
  banc.meta$root_id <- as.character(banc.meta[[.id_col]])
} else if ("root_id" %in% colnames(banc.meta)) {
  banc.meta$root_id <- as.character(banc.meta$root_id)
}

# Enrich from SeaTable if reachable (for sexually_dimorphic, malecns_cell_type, etc.)
tryCatch({
  .st <- banctable_query("SELECT root_id, sexually_dimorphic, malecns_cell_type FROM banc_meta")
  .st$root_id <- as.character(.st$root_id)
  .new_cols <- setdiff(colnames(.st), c("root_id", colnames(banc.meta)))
  .update_cols <- intersect(c("sexually_dimorphic", "malecns_cell_type"), colnames(.st))
  if (length(.update_cols) > 0 || length(.new_cols) > 0) {
    banc.meta <- banc.meta %>%
      dplyr::left_join(.st %>% dplyr::select(root_id, dplyr::all_of(c(.update_cols, .new_cols))),
                       by = "root_id", suffix = c("", ".st")) %>%
      dplyr::mutate(across(ends_with(".st"), ~ NULL))
    # Coalesce SeaTable values over GCS where SeaTable is non-NA
    for (col in .update_cols) {
      st_col <- paste0(col, ".st")
      if (st_col %in% colnames(banc.meta)) {
        banc.meta[[col]] <- dplyr::coalesce(banc.meta[[st_col]], banc.meta[[col]])
        banc.meta[[st_col]] <- NULL
      }
    }
    message("Enriched banc.meta from SeaTable (", length(.update_cols), " columns)")
  }
}, error = function(e) {
  message("SeaTable enrichment failed: ", e$message, " — using GCS-only metadata")
})

cat(sprintf("banc.meta: %d rows, %d columns\n", nrow(banc.meta), ncol(banc.meta)))

# ── BANC edgelist ──────────────────────────────────────────
message("Loading BANC edgelist...")
.edgelist_suffix <- if (.version_num >= 888) "_edgelist_simple_v2.feather" else "_edgelist_simple.feather"
.banc_el_gcs <- sprintf("%s/%s/%s%s", gcs.bucket, banc.version, banc.version, .edgelist_suffix)
banc.edgelist.simple <- read_feather_gcs(.banc_el_gcs) %>%
  dplyr::mutate(pre = as.character(pre), post = as.character(post))
# Proofread filter
.proofread_ids <- banc.meta %>%
  dplyr::filter(as.logical(proofread) %in% TRUE |
                  as.logical(roughly_proofread) %in% TRUE) %>%
  dplyr::pull(root_id)
.pre_filter <- nrow(banc.edgelist.simple)
banc.edgelist.simple <- banc.edgelist.simple %>%
  dplyr::filter(pre %in% .proofread_ids, post %in% .proofread_ids)
cat(sprintf("Edgelist: %d connections (dropped %d unproofread)\n",
            nrow(banc.edgelist.simple), .pre_filter - nrow(banc.edgelist.simple)))

# ── maleCNS data ───────────────────────────────────────────
message("Loading maleCNS data...")
.malecns_gcs <- sprintf("%s/malecns_09", gcs.bucket)
malecns.meta <- read_feather_gcs(file.path(.malecns_gcs, "malecns_09_meta.feather")) %>%
  dplyr::mutate(across(where(is.integer64), as.character))
malecns.edgelist.simple <- read_feather_gcs(file.path(.malecns_gcs, "malecns_09_simple_edgelist.feather")) %>%
  dplyr::mutate(pre = as.character(pre), post = as.character(post))
cat(sprintf("malecns.meta: %d rows, malecns.edgelist: %d connections\n",
            nrow(malecns.meta), nrow(malecns.edgelist.simple)))

# ── NBLAST data (GCS paths for Script 2) ──────────────────
nblast.gcs.path <- sprintf("%s/nblast", gsub("/compiled_data$", "", gcs.bucket))

message("=== Startup complete ===")
