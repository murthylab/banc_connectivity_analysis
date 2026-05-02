###########################################################
### Compare new vs old VNC types
###
### Part A: VNC intrinsic neurons
### Reads banc_vnc_dimorphism CSV, joins to SeaTable
### by banc ID -> root_626 to get old cell_type.
### Computes type-changed status, hemilineage changes.
###
### Part B: Effector neurons (motor + visceral_circulatory)
### Reads banc_to_manc_effectors CSV, joins manc match ID
### to franken_meta to find "new" type, assigns dimorphism
### from median dimorph score per type.
###
### Adapted from bancpipeline/banc/update/banc-vnc-type-changes.R
### Standalone — does not require bancpipeline.
###
### Inputs:  data/codex/  (drop CSVs here; missing inputs skipped)
### Outputs: data/        (CSVs)
###
### SeaTable update blocks remain commented out — uncomment to
### push back to SeaTable (requires bancr authentication).
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
  library(bancr)
  library(fafbseg)
})

message("### VNC type change analysis ###")

# --- Paths ----------------------------------------------------------------

data.dir   <- file.path(.this_repo, "data")
codex.dir  <- file.path(data.dir, "codex")  # input CSVs not in repo go here
dir.create(codex.dir, recursive = TRUE, showWarnings = FALSE)

# --- Inline helpers (replace bancpipeline / hemibrainr deps) -------------

# Status string helper: append a comma-separated tag, deduplicating + sorting.
append_status <- function(status, update) {
  vapply(status, function(s) {
    s <- paste(c(s, update), collapse = ",")
    s <- paste(sort(unique(unlist(strsplit(s, split = ",|, ")))), collapse = ",")
    gsub("^,| ", "", s)
  }, character(1), USE.NAMES = FALSE)
}

# Convenience wrapper to read CSVs only when they exist.
read_csv_if <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, ...)
}

###########################
### Read data           ###
###########################

# VNC intrinsic CSV (BANC ID + new type + dimorphism). Optional input.
vnc_csv_path <- file.path(codex.dir, "banc_vnc_cell_type_and_dimorphism.csv")
vnc <- read_csv_if(vnc_csv_path,
                    col_types = readr::cols(ID = "c"))
if (is.null(vnc)) {
  message(sprintf("  VNC dimorphism CSV not found: %s", vnc_csv_path))
  message("  Skipping VNC intrinsic analysis.")
} else {
  vnc <- vnc %>%
    dplyr::rename(banc_id = ID, new_type = Type, dimorphism = Dimorphism) %>%
    dplyr::mutate(dimorphism = gsub("^known_", "", dimorphism),
                  dimorphism = dplyr::case_when(
                    dimorphism == "Homologous" ~ "isomorphic",
                    dimorphism == "Sex Specific" ~ "female-specific",
                    dimorphism == "Dimorphic" ~ "dimorphic",
                    dimorphism == "female_specific" ~ "female-specific",
                    dimorphism == "male_specific" ~ "male-specific",
                    TRUE ~ dimorphism
                  ))
}

# BANC SeaTable metadata
bc <- bancr::banctable_query(paste0(
  "SELECT _id, root_id, root_626, supervoxel_id, super_class, flow, cell_type, side, region, ",
  "manc_cell_type, malecns_cell_type, manc_png_match, manc_nblast_match, manc_match, ",
  "malecns_match, status, sexually_dimorphic, body_part_effector, hemilineage, banc_match, ",
  "banc_match_supervoxel_id FROM banc_meta")) %>%
  dplyr::filter(!is.na(root_626), root_626 != "") %>%
  dplyr::mutate(root_626 = as.character(root_626),
                root_id = as.character(root_id),
                supervoxel_id = as.character(supervoxel_id),
                manc_png_match = as.character(manc_png_match)) %>%
  dplyr::distinct(root_626, .keep_all = TRUE)

# Recover old cell_type via manc_png_match -> franken_meta manc_id -> cell_type
fm_types <- bancr::franken_meta("SELECT manc_id, cell_type FROM franken_meta") %>%
  dplyr::filter(!is.na(manc_id), !is.na(cell_type), cell_type != "") %>%
  dplyr::mutate(manc_id = as.character(manc_id)) %>%
  dplyr::distinct(manc_id, .keep_all = TRUE)

bc <- bc %>%
  dplyr::left_join(fm_types %>% dplyr::rename(seatable_type = cell_type),
                   by = c("manc_png_match" = "manc_id"))

if (!is.null(vnc)) {
  df <- vnc %>%
    dplyr::left_join(bc, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(super_class == "ventral_nerve_cord_intrinsic") %>%
    dplyr::mutate(dimorphism = dplyr::if_else(is.na(dimorphism) | dimorphism == "",
                                               "isomorphic", dimorphism))

  # Propagate dimorphism within cell types (priority: female-specific > male-specific > dimorphic > isomorphic)
  .dim_rank <- c("female-specific" = 4, "male-specific" = 3,
                 "dimorphic" = 2, "isomorphic" = 1)
  type_level_dim <- df %>%
    dplyr::filter(!is.na(new_type), new_type != "") %>%
    dplyr::mutate(.rank = .dim_rank[dimorphism]) %>%
    dplyr::group_by(new_type) %>%
    dplyr::slice_max(.rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(new_type, dimorphism_type = dimorphism)
  n_propagated <- sum(df$new_type %in% type_level_dim$new_type &
    df$dimorphism != type_level_dim$dimorphism_type[match(df$new_type, type_level_dim$new_type)],
    na.rm = TRUE)
  df <- df %>%
    dplyr::left_join(type_level_dim, by = "new_type") %>%
    dplyr::mutate(dimorphism = dplyr::if_else(!is.na(dimorphism_type),
                                               dimorphism_type, dimorphism)) %>%
    dplyr::select(-dimorphism_type)
  if (n_propagated > 0)
    message(sprintf("  Propagated dimorphism labels to %d neurons via cell type", n_propagated))

  df <- df %>%
    dplyr::mutate(type_changed = dplyr::if_else(
      !is.na(seatable_type) & !is.na(new_type) & seatable_type != new_type,
      "type changed", "type did not change"
    ),
    effective_type = new_type)

  message(sprintf("  VNC intrinsic: %d neurons total, %d matched to seatable",
                  nrow(df), sum(!is.na(df$seatable_type))))
  message(sprintf("  Type changes: %d changed, %d unchanged",
                  sum(df$type_changed == "type changed"),
                  sum(df$type_changed == "type did not change")))
}

#########################################
### Effector neurons (motor/endocrine) ##
#########################################

eff_csv_path <- file.path(codex.dir, "id_id_matches", "banc_to_manc_effectors.csv")
eff <- read_csv_if(eff_csv_path,
                    col_types = readr::cols(`banc ID` = "c",
                                            `twin ID` = "c",
                                            `manc match ID` = "c"))
if (is.null(eff)) {
  message(sprintf("  Effector CSV not found: %s — skipping effector analysis", eff_csv_path))
} else {
  eff <- eff %>%
    dplyr::rename(banc_id = `banc ID`,
                  csv_type = type,
                  dimorph_score = `dimorph score`,
                  twin_id = `twin ID`,
                  twin_score = `twin score`,
                  manc_match_id = `manc match ID`,
                  manc_match_score = `manc match score`) %>%
    dplyr::left_join(fm_types, by = c("manc_match_id" = "manc_id")) %>%
    dplyr::rename(new_type = cell_type)

  eff_dimorph <- eff %>%
    dplyr::filter(!is.na(csv_type), csv_type != "") %>%
    dplyr::group_by(csv_type) %>%
    dplyr::summarise(median_dimorph = stats::median(dimorph_score, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(dimorphism = dplyr::case_when(
      median_dimorph > 4 ~ "sex-specific",
      median_dimorph > 2 ~ "dimorphic",
      TRUE ~ "isomorphic"
    ))

  eff <- eff %>%
    dplyr::left_join(eff_dimorph %>% dplyr::select(csv_type, dimorphism),
                     by = "csv_type") %>%
    dplyr::mutate(dimorphism = dplyr::if_else(is.na(dimorphism), "isomorphic", dimorphism)) %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, root_id, supervoxel_id, super_class,
                                           manc_nblast_match, seatable_type),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::mutate(
      effective_type = dplyr::coalesce(new_type, csv_type),
      type_changed = dplyr::case_when(
        is.na(seatable_type) & !is.na(new_type) ~ "type changed",
        !is.na(seatable_type) & !is.na(new_type) & seatable_type != new_type ~ "type changed",
        TRUE ~ "type did not change"
      )
    )

  .dim_rank <- c("female-specific" = 4, "male-specific" = 3,
                 "dimorphic" = 2, "isomorphic" = 1, "sex-specific" = 3)
  eff_type_level_dim <- eff %>%
    dplyr::filter(!is.na(effective_type), effective_type != "") %>%
    dplyr::mutate(.rank = .dim_rank[dimorphism]) %>%
    dplyr::group_by(effective_type) %>%
    dplyr::slice_max(.rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(effective_type, dimorphism_type = dimorphism)
  eff <- eff %>%
    dplyr::left_join(eff_type_level_dim, by = "effective_type") %>%
    dplyr::mutate(dimorphism = dplyr::if_else(!is.na(dimorphism_type),
                                               dimorphism_type, dimorphism)) %>%
    dplyr::select(-dimorphism_type)

  message(sprintf("  Effectors: %d neurons, %d with MANC match type, %d type changes",
                  nrow(eff),
                  sum(!is.na(eff$new_type)),
                  sum(eff$type_changed == "type changed")))
}

###########################
### Hemilineage lookup  ###
###########################

# Extract hemilineage encoded in VNC type names (e.g. IN19B057 -> 19B).
extract_hemilineage <- function(type_name) {
  result <- rep(NA_character_, length(type_name))
  valid <- !is.na(type_name)
  is_vnc <- valid & grepl("^IN\\d+[A-Z]", type_name, perl = TRUE)
  hemi <- sub("\\d+$", "", type_name[is_vnc])
  result[is_vnc] <- sub("^IN", "", hemi)
  result
}

fm <- bancr::franken_meta()

hemi_lookup <- fm %>%
  dplyr::filter(!is.na(cell_type), cell_type != "",
                !is.na(hemilineage), hemilineage != "") %>%
  dplyr::distinct(cell_type, .keep_all = TRUE) %>%
  dplyr::select(cell_type, hemilineage)

if (exists("df") && is.data.frame(df)) {
  df <- df %>%
    dplyr::mutate(
      old_hemilineage = extract_hemilineage(seatable_type),
      new_hemilineage = extract_hemilineage(new_type)
    ) %>%
    dplyr::left_join(hemi_lookup %>% dplyr::rename(old_hemilineage_fm = hemilineage),
                     by = c("seatable_type" = "cell_type")) %>%
    dplyr::left_join(hemi_lookup %>% dplyr::rename(new_hemilineage_fm = hemilineage),
                     by = c("new_type" = "cell_type")) %>%
    dplyr::mutate(
      old_hemilineage = dplyr::coalesce(old_hemilineage, old_hemilineage_fm),
      new_hemilineage = dplyr::coalesce(new_hemilineage, new_hemilineage_fm)
    ) %>%
    dplyr::select(-old_hemilineage_fm, -new_hemilineage_fm) %>%
    dplyr::mutate(hemi_changed = dplyr::case_when(
      is.na(old_hemilineage) | is.na(new_hemilineage) ~ "no hemilineage",
      old_hemilineage == new_hemilineage ~ "hemilineage did not change",
      TRUE ~ "hemilineage changed"
    ))

  message(sprintf("  Hemilineage changes: %d changed, %d unchanged, %d unknown",
                  sum(df$hemi_changed == "hemilineage changed"),
                  sum(df$hemi_changed == "hemilineage did not change"),
                  sum(df$hemi_changed == "no hemilineage")))
}

################################################
### NBLAST scores for type-changed neurons   ###
################################################

# Pull NBLAST feather from GCS
message("  Loading BANC-MANC NBLAST scores from GCS...")
gcs_nblast_path <- file.path(nblast.gcs.path, "banc_manc_v1.2.1_nblast.feather")
nblast_cache_dir <- file.path(tempdir(), "nblast_cache")
dir.create(nblast_cache_dir, showWarnings = FALSE, recursive = TRUE)
nblast_local <- file.path(nblast_cache_dir, basename(gcs_nblast_path))

if (!file.exists(nblast_local)) {
  message("  Downloading NBLAST feather from GCS...")
  system2("gsutil", c("cp", gcs_nblast_path, nblast_local),
          stdout = TRUE, stderr = TRUE)
}

nblast_available <- file.exists(nblast_local)
nb_all <- if (nblast_available) {
  arrow::read_feather(nblast_local)
} else NULL

if (nblast_available) {
  message(sprintf("  NBLAST feather: %d rows", nrow(nb_all)))

  # MANC IDs in same super_class (VNC intrinsic) for top-hit restriction
  vnc_match_ids <- fm %>%
    dplyr::filter(super_class == "ventral_nerve_cord_intrinsic",
                  !is.na(manc_id)) %>%
    dplyr::pull(manc_id) %>% as.character() %>% unique()

  # banc_id <-> root_id mapping (covers VNC intrinsic + effectors)
  id_map <- dplyr::bind_rows(
    if (exists("df") && is.data.frame(df))
      df %>% dplyr::select(banc_id, root_id) else NULL,
    if (exists("eff") && is.data.frame(eff))
      eff %>% dplyr::select(banc_id, root_id) else NULL
  ) %>%
    dplyr::filter(!is.na(root_id), root_id != "", root_id != banc_id) %>%
    dplyr::distinct(banc_id, root_id)
  rootid_to_bancid <- stats::setNames(id_map$banc_id, id_map$root_id)

  # ----- VNC intrinsic type-change scores -----
  if (exists("df") && is.data.frame(df)) {
    type_changed_df <- df %>% dplyr::filter(type_changed == "type changed")
    type_changed_ids <- unique(type_changed_df$banc_id)
    type_changed_rootids <- id_map$root_id[id_map$banc_id %in% type_changed_ids]

    nb_subset <- nb_all %>%
      dplyr::filter(root_626 %in% type_changed_ids |
                      pt_root_id %in% type_changed_rootids) %>%
      dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score) %>%
      dplyr::mutate(root_626 = dplyr::case_when(
        root_626 %in% type_changed_ids ~ root_626,
        pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
        TRUE ~ root_626
      )) %>%
      dplyr::select(-pt_root_id)

    top_scores <- nb_subset %>%
      dplyr::mutate(match_id = as.character(match_id)) %>%
      dplyr::filter(match_id %in% vnc_match_ids) %>%
      dplyr::group_by(root_626) %>%
      dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(banc_id = root_626,
                       nblast_score = score,
                       nblast_type = match_cell_type)

    st_scores <- type_changed_df %>%
      dplyr::transmute(banc_id, lookup_type = seatable_type) %>%
      dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
      dplyr::inner_join(nb_subset, by = c("banc_id" = "root_626")) %>%
      dplyr::filter(match_cell_type == lookup_type) %>%
      dplyr::group_by(banc_id) %>%
      dplyr::summarise(seatable_nblast = max(score, na.rm = TRUE), .groups = "drop")

    new_scores <- type_changed_df %>%
      dplyr::transmute(banc_id, lookup_type = effective_type) %>%
      dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
      dplyr::inner_join(nb_subset, by = c("banc_id" = "root_626")) %>%
      dplyr::filter(match_cell_type == lookup_type) %>%
      dplyr::group_by(banc_id) %>%
      dplyr::summarise(new_nblast = max(score, na.rm = TRUE), .groups = "drop")

    type_changed_df <- type_changed_df %>%
      dplyr::left_join(top_scores, by = "banc_id") %>%
      dplyr::left_join(st_scores, by = "banc_id") %>%
      dplyr::left_join(new_scores, by = "banc_id")
    message(sprintf("  Found NBLAST scores for %d/%d type-changed neurons",
                    sum(!is.na(type_changed_df$nblast_score)), nrow(type_changed_df)))
  }
} else {
  message("  NBLAST feather not available; skipping scores")
}

########################################
### MANC body ID lookup for bancsee  ###
########################################

if (nblast_available) {
  manc_lookup <- nb_all %>%
    dplyr::filter(!is.na(match_cell_type), match_cell_type != "",
                  !is.na(match_id), match_id != "") %>%
    dplyr::distinct(match_id, match_cell_type) %>%
    dplyr::rename(manc_id = match_id, cell_type = match_cell_type)
  message(sprintf("  MANC lookup: %d unique (manc_id, cell_type) pairs from NBLAST",
                  nrow(manc_lookup)))
} else {
  manc_lookup <- fm %>%
    dplyr::filter(!is.na(manc_id), !is.na(cell_type), cell_type != "") %>%
    dplyr::distinct(manc_id, cell_type)
}

# Fallback: use manc_nblast_match from SeaTable to fill nblast_type where the feather had no hit
if (exists("type_changed_df") && "manc_nblast_match" %in% names(type_changed_df)) {
  nblast_match_types <- manc_lookup %>%
    dplyr::distinct(manc_id, .keep_all = TRUE) %>%
    dplyr::rename(nblast_type_fb = cell_type)
  type_changed_df <- type_changed_df %>%
    dplyr::left_join(nblast_match_types,
                     by = c("manc_nblast_match" = "manc_id")) %>%
    dplyr::mutate(nblast_type = dplyr::coalesce(nblast_type, nblast_type_fb)) %>%
    dplyr::select(-nblast_type_fb)
}

###########################################
### Neuroglancer URLs for type changes  ###
###########################################

build_ngl_links <- function(df_rows, ngl_url, lookup, type_cols = c("seatable_type", "effective_type", "nblast_type")) {
  ngl_url2 <- sub("#!middleauth+", "?", ngl_url, fixed = TRUE)
  ngl_parts <- unlist(strsplit(ngl_url2, "?", fixed = TRUE))
  ngl_json <- fafbseg::flywire_fetch(ngl_parts[2], token = bancr:::banc_token(),
                                      return = "text", cache = TRUE)
  ngl_base <- fafbseg::ngl_decode_scene(
    fafbseg::ngl_encode_url(ngl_json, baseurl = ngl_parts[1]))

  ngl_ls <- fafbseg:::ngl_layer_summary(fafbseg::ngl_layers(ngl_base))
  banc_layer_idx <- match("v626 neurons", ngl_ls$name)
  if (is.na(banc_layer_idx)) banc_layer_idx <- match("segmentation proofreading", ngl_ls$name)

  manc_layer_idxs <- grep("manc", ngl_ls$name, ignore.case = TRUE)
  manc_old_idx <- manc_layer_idxs[1]
  manc_new_idx <- if (length(manc_layer_idxs) >= 2) manc_layer_idxs[2] else manc_layer_idxs[1]
  manc_nblast_idx <- grep("nblast", ngl_ls$name, ignore.case = TRUE)
  manc_nblast_idx <- if (length(manc_nblast_idx) > 0) manc_nblast_idx[1] else NA_integer_

  vapply(seq_len(nrow(df_rows)), function(i) {
    row <- df_rows[i, ]
    banc_rid <- row$banc_id
    type_to_ids <- function(t) {
      if (is.na(t) || t == "") return(character(0))
      lookup %>% dplyr::filter(cell_type == t) %>% dplyr::pull(manc_id)
    }
    st_manc <- type_to_ids(row[[type_cols[1]]])
    new_manc <- type_to_ids(row[[type_cols[2]]])
    nblast_manc <- type_to_ids(row[[type_cols[3]]])

    tryCatch({
      sc <- ngl_base
      sc[["layers"]][[banc_layer_idx]][["segments"]] <- as.character(banc_rid)
      sc[["layers"]][[banc_layer_idx]][["hiddenSegments"]] <- NULL
      sc[["layers"]][[manc_old_idx]][["segments"]] <- as.character(st_manc)
      sc[["layers"]][[manc_old_idx]][["hiddenSegments"]] <- NULL
      sc[["layers"]][[manc_new_idx]][["segments"]] <- as.character(new_manc)
      sc[["layers"]][[manc_new_idx]][["hiddenSegments"]] <- NULL
      if (!is.na(manc_nblast_idx) && length(nblast_manc) > 0) {
        sc[["layers"]][[manc_nblast_idx]][["segments"]] <- as.character(nblast_manc)
        sc[["layers"]][[manc_nblast_idx]][["hiddenSegments"]] <- NULL
      }
      as.character(sc)
    }, error = function(e) NA_character_)
  }, character(1))
}

ngl_url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6384965498437632"

if (exists("type_changed_df") && nrow(type_changed_df) > 0) {
  message(sprintf("  Building neuroglancer links for %d type-changed neurons...",
                  nrow(type_changed_df)))
  type_changed_df$neuroglancer_url <- tryCatch(
    build_ngl_links(type_changed_df, ngl_url, manc_lookup),
    error = function(e) { warning("Could not build neuroglancer links: ", e$message)
                          rep(NA_character_, nrow(type_changed_df)) })
  message(sprintf("  Generated %d/%d neuroglancer links",
                  sum(!is.na(type_changed_df$neuroglancer_url)), nrow(type_changed_df)))
}

#############################################
### Save CSVs: all changes + hemilineage  ###
#############################################

if (exists("type_changed_df") && nrow(type_changed_df) > 0) {
  all_changes_csv <- type_changed_df %>%
    dplyr::transmute(
      banc_id, root_id, seatable_type, effective_type, nblast_type,
      old_hemilineage, new_hemilineage, hemi_changed, dimorphism,
      nblast_score, seatable_nblast, new_nblast, neuroglancer_url
    ) %>%
    dplyr::arrange(dplyr::desc(seatable_nblast - new_nblast))

  csv_path_all <- file.path(data.dir, "vnc_type_changes.csv")
  readr::write_csv(all_changes_csv, csv_path_all)
  message(sprintf("  Saved all type changes CSV: %s (%d rows)",
                  csv_path_all, nrow(all_changes_csv)))

  hemi_csv <- all_changes_csv %>%
    dplyr::filter(hemi_changed == "hemilineage changed") %>%
    dplyr::select(-hemi_changed) %>%
    dplyr::arrange(dplyr::desc(new_nblast))
  csv_path_hemi <- file.path(data.dir, "vnc_hemilineage_changes.csv")
  readr::write_csv(hemi_csv, csv_path_hemi)
  message(sprintf("  Saved hemilineage changes CSV: %s (%d rows)",
                  csv_path_hemi, nrow(hemi_csv)))
}

# ---- Effector type changes CSV ----
if (exists("eff") && is.data.frame(eff)) {
  eff_type_changed <- eff %>% dplyr::filter(type_changed == "type changed")

  if (nblast_available && nrow(eff_type_changed) > 0) {
    eff_tc_ids <- unique(eff_type_changed$banc_id)
    eff_tc_rootids <- id_map$root_id[id_map$banc_id %in% eff_tc_ids]

    nb_eff <- nb_all %>%
      dplyr::filter(root_626 %in% eff_tc_ids |
                      pt_root_id %in% eff_tc_rootids) %>%
      dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score) %>%
      dplyr::mutate(root_626 = dplyr::case_when(
        root_626 %in% eff_tc_ids ~ root_626,
        pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
        TRUE ~ root_626
      )) %>%
      dplyr::select(-pt_root_id)

    eff_top <- nb_eff %>%
      dplyr::group_by(root_626) %>%
      dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(banc_id = root_626,
                       nblast_score = score,
                       nblast_type = match_cell_type)
    eff_st <- eff_type_changed %>%
      dplyr::transmute(banc_id, lookup_type = seatable_type) %>%
      dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
      dplyr::inner_join(nb_eff, by = c("banc_id" = "root_626")) %>%
      dplyr::filter(match_cell_type == lookup_type) %>%
      dplyr::group_by(banc_id) %>%
      dplyr::summarise(seatable_nblast = max(score, na.rm = TRUE), .groups = "drop")
    eff_new <- eff_type_changed %>%
      dplyr::transmute(banc_id, lookup_type = new_type) %>%
      dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
      dplyr::inner_join(nb_eff, by = c("banc_id" = "root_626")) %>%
      dplyr::filter(match_cell_type == lookup_type) %>%
      dplyr::group_by(banc_id) %>%
      dplyr::summarise(new_nblast = max(score, na.rm = TRUE), .groups = "drop")
    eff_type_changed <- eff_type_changed %>%
      dplyr::left_join(eff_top, by = "banc_id") %>%
      dplyr::left_join(eff_st, by = "banc_id") %>%
      dplyr::left_join(eff_new, by = "banc_id")
  } else {
    eff_type_changed$nblast_score    <- NA_real_
    eff_type_changed$nblast_type     <- NA_character_
    eff_type_changed$seatable_nblast <- NA_real_
    eff_type_changed$new_nblast      <- NA_real_
  }

  if ("manc_nblast_match" %in% names(eff_type_changed)) {
    eff_nblast_fb <- manc_lookup %>%
      dplyr::distinct(manc_id, .keep_all = TRUE) %>%
      dplyr::rename(nblast_type_fb = cell_type)
    eff_type_changed <- eff_type_changed %>%
      dplyr::left_join(eff_nblast_fb,
                       by = c("manc_nblast_match" = "manc_id")) %>%
      dplyr::mutate(nblast_type = dplyr::coalesce(nblast_type, nblast_type_fb)) %>%
      dplyr::select(-nblast_type_fb)
  }

  if (nrow(eff_type_changed) > 0) {
    eff_type_changed$neuroglancer_url <- tryCatch(
      build_ngl_links(eff_type_changed, ngl_url, manc_lookup),
      error = function(e) rep(NA_character_, nrow(eff_type_changed)))
    eff_changes_csv <- eff_type_changed %>%
      dplyr::transmute(
        banc_id, root_id, seatable_type, effective_type, nblast_type,
        manc_match_id, manc_match_score, dimorphism, dimorph_score,
        nblast_score, seatable_nblast, new_nblast, neuroglancer_url
      ) %>%
      dplyr::arrange(dplyr::desc(new_nblast))
    csv_path_eff <- file.path(data.dir, "vnc_effector_type_changes.csv")
    readr::write_csv(eff_changes_csv, csv_path_eff)
    message(sprintf("  Saved effector type changes CSV: %s (%d rows)",
                    csv_path_eff, nrow(eff_changes_csv)))
  }
}

#########################################################
### CSV: maleCNS effector type changes assessment     ###
#########################################################

mcns_eff_csv_path <- file.path(codex.dir, "id_id_matches", "banc_to_mcns_effectors.csv")
mcns_eff <- read_csv_if(mcns_eff_csv_path,
                         col_types = readr::cols(`banc ID` = "c",
                                                 `twin ID` = "c",
                                                 `mcns match ID` = "c"))
if (is.null(mcns_eff)) {
  message(sprintf("  maleCNS effector CSV not found: %s — skipping", mcns_eff_csv_path))
} else {
  message("  === maleCNS effector type changes ===")
  mcns_eff <- mcns_eff %>%
    dplyr::rename(banc_id = `banc ID`,
                  csv_type = type,
                  dimorph_score = `dimorph score`,
                  twin_id = `twin ID`,
                  twin_score = `twin score`,
                  mcns_match_id = `mcns match ID`,
                  mcns_match_score = `mcns match score`)

  mcns_meta_local <- malecns.meta %>%
    dplyr::mutate(malecns_09_id = as.character(malecns_09_id)) %>%
    dplyr::distinct(malecns_09_id, .keep_all = TRUE) %>%
    dplyr::select(malecns_09_id, cell_type)

  mcns_eff <- mcns_eff %>%
    dplyr::left_join(mcns_meta_local,
                     by = c("mcns_match_id" = "malecns_09_id")) %>%
    dplyr::rename(new_malecns_type = cell_type)

  mcns_eff_dimorph <- mcns_eff %>%
    dplyr::filter(!is.na(csv_type), csv_type != "") %>%
    dplyr::group_by(csv_type) %>%
    dplyr::summarise(median_dimorph = stats::median(dimorph_score, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(dimorphism = dplyr::case_when(
      median_dimorph > 4 ~ "sex-specific",
      median_dimorph > 2 ~ "dimorphic",
      TRUE ~ "isomorphic"
    ))

  mcns_eff <- mcns_eff %>%
    dplyr::left_join(mcns_eff_dimorph %>% dplyr::select(csv_type, dimorphism),
                     by = "csv_type") %>%
    dplyr::mutate(dimorphism = dplyr::if_else(is.na(dimorphism), "isomorphic", dimorphism)) %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, root_id, supervoxel_id, super_class,
                                           current_malecns_type = malecns_cell_type),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::mutate(type_changed = dplyr::case_when(
      is.na(current_malecns_type) & !is.na(new_malecns_type) ~ "type changed",
      !is.na(current_malecns_type) & !is.na(new_malecns_type) &
        current_malecns_type != new_malecns_type ~ "type changed",
      TRUE ~ "type did not change"
    ))

  mcns_type_changed <- mcns_eff %>% dplyr::filter(type_changed == "type changed")
  message(sprintf("  maleCNS effectors: %d neurons, %d type changes",
                  nrow(mcns_eff), nrow(mcns_type_changed)))

  # NBLAST scores from BANC-maleCNS feather
  mcns_nblast_path <- file.path(nblast.gcs.path, "banc_malecns_v0.9_nblast.feather")
  mcns_nblast_local <- file.path(nblast_cache_dir, basename(mcns_nblast_path))
  if (!file.exists(mcns_nblast_local)) {
    message("  Downloading maleCNS NBLAST feather from GCS...")
    system2("gsutil", c("cp", mcns_nblast_path, mcns_nblast_local),
            stdout = TRUE, stderr = TRUE)
  }

  mcns_nblast_available <- file.exists(mcns_nblast_local)
  if (mcns_nblast_available && nrow(mcns_type_changed) > 0) {
    nb_mcns <- arrow::read_feather(mcns_nblast_local)
    message(sprintf("  maleCNS NBLAST feather: %d rows", nrow(nb_mcns)))

    mcns_tc_ids <- unique(mcns_type_changed$banc_id)
    mcns_tc_rootids <- if (exists("id_map"))
      id_map$root_id[id_map$banc_id %in% mcns_tc_ids] else character(0)

    nb_mcns_sub <- nb_mcns %>%
      dplyr::filter(root_626 %in% mcns_tc_ids |
                      pt_root_id %in% mcns_tc_rootids) %>%
      dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score)
    if (exists("rootid_to_bancid")) {
      nb_mcns_sub <- nb_mcns_sub %>%
        dplyr::mutate(root_626 = dplyr::case_when(
          root_626 %in% mcns_tc_ids ~ root_626,
          pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
          TRUE ~ root_626
        ))
    }
    nb_mcns_sub <- nb_mcns_sub %>% dplyr::select(-pt_root_id)

    mcns_top <- nb_mcns_sub %>%
      dplyr::group_by(root_626) %>%
      dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(banc_id = root_626,
                       nblast_score = score,
                       nblast_type = match_cell_type)
    mcns_st <- mcns_type_changed %>%
      dplyr::transmute(banc_id, lookup_type = current_malecns_type) %>%
      dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
      dplyr::inner_join(nb_mcns_sub, by = c("banc_id" = "root_626")) %>%
      dplyr::filter(match_cell_type == lookup_type) %>%
      dplyr::group_by(banc_id) %>%
      dplyr::summarise(seatable_nblast = max(score, na.rm = TRUE), .groups = "drop")
    mcns_new <- mcns_type_changed %>%
      dplyr::transmute(banc_id, lookup_type = new_malecns_type) %>%
      dplyr::filter(!is.na(lookup_type), lookup_type != "") %>%
      dplyr::inner_join(nb_mcns_sub, by = c("banc_id" = "root_626")) %>%
      dplyr::filter(match_cell_type == lookup_type) %>%
      dplyr::group_by(banc_id) %>%
      dplyr::summarise(new_nblast = max(score, na.rm = TRUE), .groups = "drop")

    mcns_type_changed <- mcns_type_changed %>%
      dplyr::left_join(mcns_top, by = "banc_id") %>%
      dplyr::left_join(mcns_st, by = "banc_id") %>%
      dplyr::left_join(mcns_new, by = "banc_id")
    rm(nb_mcns, nb_mcns_sub); gc()
  } else {
    mcns_type_changed$nblast_score    <- NA_real_
    mcns_type_changed$nblast_type     <- NA_character_
    mcns_type_changed$seatable_nblast <- NA_real_
    mcns_type_changed$new_nblast      <- NA_real_
  }

  # maleCNS body lookup (for ngl links)
  mcns_lookup <- malecns.meta %>%
    dplyr::filter(!is.na(cell_type), cell_type != "") %>%
    dplyr::mutate(mcns_id = as.character(malecns_09_id)) %>%
    dplyr::select(mcns_id, cell_type) %>%
    dplyr::rename(manc_id = mcns_id)  # reuse builder column name

  mcns_ngl_url <- "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6343573275410432"
  if (nrow(mcns_type_changed) > 0) {
    mcns_type_changed$neuroglancer_url <- tryCatch({
      # Reuse builder; rename input columns to match
      tmp <- mcns_type_changed %>%
        dplyr::rename(seatable_type = current_malecns_type,
                      effective_type = new_malecns_type)
      build_ngl_links(tmp, mcns_ngl_url, mcns_lookup)
    }, error = function(e) {
      warning("Could not build maleCNS neuroglancer links: ", e$message)
      rep(NA_character_, nrow(mcns_type_changed))
    })

    mcns_changes_csv <- mcns_type_changed %>%
      dplyr::transmute(
        banc_id, root_id,
        current_malecns_type, new_malecns_type, nblast_type,
        mcns_match_id, mcns_match_score, dimorphism, dimorph_score,
        nblast_score, seatable_nblast, new_nblast, neuroglancer_url
      ) %>%
      dplyr::arrange(dplyr::desc(new_nblast))
    csv_path_mcns <- file.path(data.dir, "mcns_effector_type_changes.csv")
    readr::write_csv(mcns_changes_csv, csv_path_mcns)
    message(sprintf("  Saved maleCNS effector type changes CSV: %s (%d rows)",
                    csv_path_mcns, nrow(mcns_changes_csv)))
  }
}

###############################################################
### SeaTable update blocks (commented out — uncomment to push)
###############################################################
# The original bancpipeline script contains several large blocks that push
# updates back to SeaTable for reviewed accept_new == "T" / "F" / "A".
# Those have been omitted from this standalone analysis script. If you need
# to re-push reviewed assignments, copy the blocks from bancpipeline and
# uncomment the banctable_update_rows() calls.

###########################################################
### Reviewed type changes (reads in-repo CSVs)           ###
###########################################################

reviewed_csv_path <- file.path(data.dir, "vnc_type_changes_reveiwed.csv")
if (file.exists(reviewed_csv_path)) {
  message("  Reading reviewed type changes...")
  reviewed <- readr::read_csv(reviewed_csv_path,
                               col_types = readr::cols(banc_id = "c", root_id = "c",
                                                        .default = "c"),
                               show_col_types = FALSE) %>%
    dplyr::mutate(
      accept_new = trimws(accept_new),
      seatable_nblast = as.numeric(seatable_nblast),
      new_nblast = as.numeric(new_nblast),
      nblast_score = as.numeric(nblast_score)
    ) %>%
    dplyr::filter(accept_new %in% c("T", "F", "A")) %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, super_class, sexually_dimorphic),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::mutate(super_class = dplyr::if_else(is.na(super_class), "unknown", super_class))

  effector_classes <- c("motor", "visceral_circulatory", "endocrine")
  n_eff_override <- sum(reviewed$super_class %in% effector_classes &
                          reviewed$sexually_dimorphic == "female-specific" &
                          reviewed$accept_new != "A", na.rm = TRUE)
  reviewed <- reviewed %>%
    dplyr::mutate(accept_new = dplyr::if_else(
      super_class %in% effector_classes & sexually_dimorphic == "female-specific",
      "A", accept_new))
  reviewed$accept_new <- factor(reviewed$accept_new,
    levels = c("T", "F", "A"),
    labels = c("Connectivity match better", "Morphology match better", "All wrong"))

  message(sprintf("  %d reviewed neurons: %s",
                  nrow(reviewed),
                  paste(table(reviewed$accept_new), names(table(reviewed$accept_new)),
                        collapse = ", ")))
  if (n_eff_override > 0) {
    message(sprintf("  %d effector neurons overridden to 'All wrong' (female-specific in BANC)",
                    n_eff_override))
  }
} else {
  message("  Reviewed CSV not found at ", reviewed_csv_path)
}

###########################################################
### Left-right mirror matches (banc_match)               ###
###########################################################

lr_csv_path <- file.path(codex.dir, "banc_vnc_intrinsic_left_right_matches.csv")
lr_matches <- read_csv_if(lr_csv_path,
                           col_types = readr::cols(`Ego ID` = "c", `Twin ID` = "c", Conf = "d"))
if (is.null(lr_matches)) {
  message(sprintf("  LR matches CSV not found: %s — skipping", lr_csv_path))
} else {
  lr_matches <- lr_matches %>%
    dplyr::rename(ego_id = `Ego ID`, twin_id = `Twin ID`, conf = Conf)
  message(sprintf("  LR matches CSV: %d pairs", nrow(lr_matches)))

  ego_lookup <- bc %>%
    dplyr::select(`_id`, root_626, root_id) %>%
    dplyr::distinct(root_626, .keep_all = TRUE)
  twin_lookup <- bc %>%
    dplyr::select(root_626, root_id, supervoxel_id) %>%
    dplyr::distinct(root_626, .keep_all = TRUE) %>%
    dplyr::rename(twin_root_id = root_id, twin_svid = supervoxel_id)

  lr_update <- lr_matches %>%
    dplyr::inner_join(ego_lookup, by = c("ego_id" = "root_626")) %>%
    dplyr::inner_join(twin_lookup, by = c("twin_id" = "root_626")) %>%
    dplyr::select(`_id`, root_id, banc_match = twin_root_id,
                  banc_match_supervoxel_id = twin_svid) %>%
    dplyr::filter(!is.na(banc_match), banc_match != "",
                  !is.na(banc_match_supervoxel_id), banc_match_supervoxel_id != "") %>%
    as.data.frame()

  message(sprintf("  %d/%d pairs mapped to seatable rows with valid banc_match + supervoxel_id",
                  nrow(lr_update), nrow(lr_matches)))
  # banctable_update_rows(base = 'banc_meta', table = "banc_meta",
  #                       df = lr_update, append_allowed = FALSE, chunksize = 1000)
}

###########################################################
### Summary statistics                                   ###
###########################################################

if (exists("df") && is.data.frame(df)) {
  message("\n=== Summary: VNC intrinsic neurons ===")
  bc_vnc <- bc %>% dplyr::filter(super_class == "ventral_nerve_cord_intrinsic")
  n_vnc_st <- nrow(bc_vnc)
  n_Q <- bc_vnc %>% dplyr::filter(!is.na(manc_png_match), manc_png_match != "") %>% nrow()
  p_df <- bc_vnc %>%
    dplyr::filter(!is.na(seatable_type), seatable_type != "",
                  !is.na(manc_nblast_match), manc_nblast_match != "") %>%
    dplyr::mutate(manc_nblast_match = as.character(manc_nblast_match)) %>%
    dplyr::left_join(fm_types %>% dplyr::rename(nblast_match_type = cell_type),
                     by = c("manc_nblast_match" = "manc_id"))
  n_P <- sum(!is.na(p_df$nblast_match_type) &
             p_df$seatable_type != p_df$nblast_match_type, na.rm = TRUE)
  n_J <- nrow(df)
  n_K <- df %>%
    dplyr::filter(type_changed == "type did not change",
                  !is.na(seatable_type), seatable_type != "") %>% nrow()
  n_L <- sum(df$type_changed == "type changed")
  n_N <- df %>%
    dplyr::filter(!is.na(new_type), new_type != "",
                  is.na(seatable_type) | seatable_type == "") %>% nrow()

  vnc_changed_ids <- df %>% dplyr::filter(type_changed == "type changed") %>% dplyr::pull(banc_id)
  if (exists("reviewed") && is.data.frame(reviewed)) {
    vnc_rev <- reviewed %>% dplyr::filter(banc_id %in% vnc_changed_ids)
    n_X <- sum(vnc_rev$accept_new == "Connectivity match better")
    n_Y <- sum(vnc_rev$accept_new == "Morphology match better")
    n_R <- sum(vnc_rev$accept_new == "All wrong")
  } else {
    n_X <- n_Y <- n_R <- NA_integer_
  }
  H <- round(100 * (n_K + n_X) / (n_K + n_L), 1)
  G <- round(100 * (n_X + n_N) / n_J, 1)

  message(sprintf(paste0(
    "\n%d/%d (%.0f%%) ventral nerve cord intrinsic neurons were matched to a MANC ",
    "cell type. Of these, %d differed from their top NBLAST hit. Connectivity matching ",
    "suggested a match for %d VNC intrinsic neurons (%d already had a reviewed morphology ",
    "match, %d were type-change suggestions, %d were new matches). Of the %d type-change ",
    "suggestions, %d confirmed connectivity, %d preferred morphology, %d were poor in both. ",
    "Total agreement rate: %.1f%% [(K+X)/(K+L)], improvement: %.1f%% [(X+N)/J]."),
    n_Q, n_vnc_st, 100 * n_Q / n_vnc_st, n_P,
    n_J, n_K, n_L, n_N, n_L, n_X, n_Y, n_R, H, G))
}

message("### VNC type change analysis complete ###")
