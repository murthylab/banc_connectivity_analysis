###########################################################
### Figure panels: VNC morphology vs connectivity types
###
### Standalone plotting script — all data comes from
### SeaTable (already pushed) and GCS NBLAST feathers.
### Reviewed CSVs from ../data/.
###
### Panels:
###   A) vnc_reviewed_accept_counts
###   B) vnc_type_agreement_by_dimorphism
###   C) vnc_nblast_violin
###   D) vnc_reviewed_nblast_scatter
###   E) vnc_nblast_rank_manc
###   F) vnc_nblast_rank_malecns
###
### Additional plots:
###   vnc_type_changes_by_dimorphism
###   vnc_nblast_density
###   vnc_type_agreement (unfaceted)
###   lr_mirror_match_type_agreement
###
### Usage: Rscript panels-vnc-morphology-connectivity-types.R
###########################################################
source(file.path(dirname(sys.frame(1)$ofile), "startup.R"))

library(ggplot2)
library(dplyr)

message("### Figure: VNC morphology vs connectivity types ###")

.this_repo <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = FALSE)
if (!nzchar(.this_repo) || !dir.exists(.this_repo)) .this_repo <- getwd()
plot.dir <- file.path(.this_repo, "figs/figure_typing/links")
data.dir <- file.path(.this_repo, "data")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

###########################################################
### 1. Load data from SeaTable
###########################################################

message("\n=== Loading data from SeaTable ===")

bc <- banctable_query(paste0(
  "SELECT _id, root_id, root_626, supervoxel_id, super_class, cell_type, side, ",
  "manc_cell_type, malecns_cell_type, manc_png_match, manc_nblast_match, ",
  "sexually_dimorphic FROM banc_meta")) %>%
  dplyr::filter(!is.na(root_626), root_626 != "") %>%
  dplyr::mutate(root_626 = as.character(root_626),
                root_id = as.character(root_id),
                supervoxel_id = as.character(supervoxel_id),
                manc_png_match = as.character(manc_png_match)) %>%
  dplyr::distinct(root_626, .keep_all = TRUE)
message(sprintf("  BANC: %d neurons", nrow(bc)))

# franken_meta for MANC type lookups and side matching
fm <- franken_meta()
message(sprintf("  franken_meta: %d rows", nrow(fm)))

# Recover morphology-based type: manc_png_match -> franken_meta cell_type
fm_types <- fm %>%
  dplyr::filter(!is.na(manc_id), !is.na(cell_type), cell_type != "") %>%
  dplyr::mutate(manc_id = as.character(manc_id)) %>%
  dplyr::distinct(manc_id, .keep_all = TRUE) %>%
  dplyr::select(manc_id, morphology_type = cell_type)

bc <- bc %>%
  dplyr::left_join(fm_types, by = c("manc_png_match" = "manc_id"))

# VNC intrinsic neurons: connectivity type = manc_cell_type (pushed from CSV)
#                         morphology type  = morphology_type (from PNG match)
df <- bc %>%
  dplyr::filter(super_class == "ventral_nerve_cord_intrinsic",
                !is.na(manc_cell_type), manc_cell_type != "") %>%
  dplyr::mutate(
    connectivity_type = manc_cell_type,
    dimorphism = dplyr::if_else(
      is.na(sexually_dimorphic) | sexually_dimorphic == "",
      "isomorphic", sexually_dimorphic),
    type_changed = dplyr::if_else(
      !is.na(morphology_type) & morphology_type != "" &
        morphology_type != connectivity_type,
      "type changed", "type did not change")
  )
message(sprintf("  VNC intrinsic with connectivity type: %d neurons", nrow(df)))
message(sprintf("  Type changed: %d, unchanged: %d",
                sum(df$type_changed == "type changed"),
                sum(df$type_changed == "type did not change")))

###########################################################
### 2. Load NBLAST data from GCS
###########################################################

message("\n=== Loading NBLAST data ===")

nblast_cache_dir <- file.path(tempdir(), "nblast_cache")
dir.create(nblast_cache_dir, showWarnings = FALSE, recursive = TRUE)

# BANC-MANC NBLAST
gcs_nblast_path <- file.path(nblast.gcs.path, "banc_manc_v1.2.1_nblast.feather")
nblast_local <- file.path(nblast_cache_dir, basename(gcs_nblast_path))
if (!file.exists(nblast_local)) {
  message("  Downloading MANC NBLAST feather from GCS...")
  system2("gsutil", c("cp", gcs_nblast_path, nblast_local),
          stdout = TRUE, stderr = TRUE)
}

nblast_available <- file.exists(nblast_local)
if (nblast_available) {
  nb_all <- arrow::read_feather(nblast_local)
  message(sprintf("  MANC NBLAST: %d rows", nrow(nb_all)))
} else {
  message("  MANC NBLAST feather not available")
}

###########################################################
### 3. Prepare NBLAST score data
###########################################################

if (nblast_available) {
  message("\n=== Preparing NBLAST scores ===")

  # MANC IDs from same super_class (VNC intrinsic)
  vnc_match_ids <- fm %>%
    dplyr::filter(super_class == "ventral_nerve_cord_intrinsic",
                  !is.na(manc_id)) %>%
    dplyr::pull(manc_id) %>% as.character() %>% unique()

  # Map root_626 ↔ root_id
  id_map <- df %>%
    dplyr::select(banc_id = root_626, root_id) %>%
    dplyr::filter(!is.na(root_id), root_id != "", root_id != banc_id) %>%
    dplyr::distinct(banc_id, root_id)
  rootid_to_bancid <- stats::setNames(id_map$banc_id, id_map$root_id)

  vnc_ids <- unique(df$root_626)
  vnc_rootids <- id_map$root_id[id_map$banc_id %in% vnc_ids]

  # Subset NBLAST to VNC neurons
  nb_vnc <- nb_all %>%
    dplyr::filter(root_626 %in% vnc_ids |
                    pt_root_id %in% vnc_rootids) %>%
    dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score) %>%
    dplyr::mutate(root_626 = dplyr::case_when(
      root_626 %in% vnc_ids ~ root_626,
      pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
      TRUE ~ root_626
    )) %>%
    dplyr::select(-pt_root_id)
  message(sprintf("  %d NBLAST rows for VNC neurons", nrow(nb_vnc)))

  # (a) Top overall score per neuron (same super_class)
  top_scores <- nb_vnc %>%
    dplyr::mutate(match_id = as.character(match_id)) %>%
    dplyr::filter(match_id %in% vnc_match_ids) %>%
    dplyr::group_by(root_626) %>%
    dplyr::summarise(top_nblast = max(score, na.rm = TRUE), .groups = "drop")

  # (b) Top score for morphology type
  st_scores <- df %>%
    dplyr::filter(!is.na(morphology_type), morphology_type != "") %>%
    dplyr::transmute(banc_id = root_626, st_type = morphology_type) %>%
    dplyr::inner_join(nb_vnc, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(match_cell_type == st_type) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::summarise(st_nblast = max(score, na.rm = TRUE), .groups = "drop")

  # (c) Top score for connectivity type
  new_scores <- df %>%
    dplyr::transmute(banc_id = root_626, eff_type = connectivity_type) %>%
    dplyr::filter(!is.na(eff_type), eff_type != "") %>%
    dplyr::inner_join(nb_vnc, by = c("banc_id" = "root_626")) %>%
    dplyr::filter(match_cell_type == eff_type) %>%
    dplyr::group_by(banc_id) %>%
    dplyr::summarise(new_nblast = max(score, na.rm = TRUE), .groups = "drop")

  nblast_all_df <- top_scores %>%
    dplyr::left_join(st_scores, by = c("root_626" = "banc_id")) %>%
    dplyr::left_join(new_scores, by = c("root_626" = "banc_id")) %>%
    dplyr::left_join(df %>% dplyr::distinct(root_626, dimorphism),
                     by = "root_626")

  message(sprintf("  %d with top score, %d with morphology score, %d with connectivity score",
                  sum(!is.na(nblast_all_df$top_nblast)),
                  sum(!is.na(nblast_all_df$st_nblast)),
                  sum(!is.na(nblast_all_df$new_nblast))))

  # Top NBLAST hit cell type per neuron
  top_hit_type <- nb_vnc %>%
    dplyr::group_by(root_626) %>%
    dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(root_626, top_nblast_type = match_cell_type)

  # Agreement table
  agree_df <- df %>%
    dplyr::select(banc_id = root_626, connectivity_type, morphology_type, dimorphism) %>%
    dplyr::left_join(top_hit_type, by = c("banc_id" = "root_626"))
}

###########################################################
### 4. Read reviewed CSVs
###########################################################

message("\n=== Loading reviewed data ===")

reviewed_csv_path <- file.path(data.dir, "vnc_type_changes_reveiwed.csv")
has_reviewed <- file.exists(reviewed_csv_path)
if (has_reviewed) {
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
    dplyr::filter(accept_new %in% c("T", "F", "A"))

  # Join super_class and sexually_dimorphic from SeaTable
  reviewed <- reviewed %>%
    dplyr::left_join(bc %>% dplyr::select(root_626, super_class, sexually_dimorphic),
                     by = c("banc_id" = "root_626")) %>%
    dplyr::mutate(super_class = dplyr::if_else(is.na(super_class), "unknown", super_class))

  # Effector female-specific override
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
} else {
  message("  Reviewed CSV not found: ", reviewed_csv_path)
}

###########################################################
### Panel A: Reviewed accept counts
###########################################################

if (has_reviewed) {
  message("\n=== Panel A: Reviewed accept counts ===")

  clean_label <- function(x) gsub("_", " ", x)
  reviewed_prop <- reviewed %>%
    dplyr::mutate(super_class = clean_label(super_class)) %>%
    dplyr::count(super_class, accept_new) %>%
    dplyr::group_by(super_class) %>%
    dplyr::mutate(prop = n / sum(n), total = sum(n)) %>%
    dplyr::ungroup()

  sc_order <- reviewed_prop %>%
    dplyr::filter(accept_new == "Connectivity match better") %>%
    dplyr::arrange(dplyr::desc(prop)) %>%
    dplyr::pull(super_class)
  sc_order <- c(sc_order, setdiff(unique(reviewed_prop$super_class), sc_order))
  reviewed_prop$super_class <- factor(reviewed_prop$super_class, levels = rev(sc_order))

  sc_totals <- reviewed_prop %>% dplyr::distinct(super_class, total)

  plot_reviewed_counts <- ggplot(reviewed_prop,
                          aes(x = prop, y = super_class, fill = accept_new)) +
    geom_col(color = NA) +
    geom_text(data = sc_totals,
              aes(x = 1.02, y = super_class, label = paste0("n=", total), fill = NULL),
              hjust = 0, size = 3, inherit.aes = FALSE) +
    scale_fill_manual(values = paper.cols) +
    scale_x_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Proportion", y = NULL, fill = NULL) +
    theme_minimal() +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())

  ggsave(file.path(plot.dir, "vnc_reviewed_accept_counts.pdf"),
         plot = plot_reviewed_counts, width = 4, height = 3, bg = "white")
  message("  Saved: vnc_reviewed_accept_counts.pdf")
}

###########################################################
### Panel B: Type agreement by dimorphism
###########################################################

if (nblast_available) {
  message("\n=== Panel B: Type agreement by dimorphism ===")

  compute_agreement <- function(adf) {
    n_nb <- sum(!is.na(adf$top_nblast_type))
    n_st_nb <- sum(!is.na(adf$morphology_type) & !is.na(adf$top_nblast_type))
    n_st_new <- sum(!is.na(adf$morphology_type) & !is.na(adf$connectivity_type))

    cnt_new_ne_nblast <- sum(adf$connectivity_type != adf$top_nblast_type, na.rm = TRUE)
    cnt_st_ne_nblast <- sum(!is.na(adf$morphology_type) &
                              adf$morphology_type != adf$top_nblast_type, na.rm = TRUE)
    cnt_st_ne_new <- sum(!is.na(adf$morphology_type) &
                           adf$morphology_type != adf$connectivity_type, na.rm = TRUE)

    data.frame(
      category = factor(c("connectivity != nblast", "morphology != nblast",
                           "morphology != connectivity"),
                         levels = c("connectivity != nblast", "morphology != nblast",
                                    "morphology != connectivity")),
      pct = c(cnt_new_ne_nblast / n_nb * 100,
              cnt_st_ne_nblast / n_st_nb * 100,
              cnt_st_ne_new / n_st_new * 100),
      n_label = c(sprintf("%d / %d", cnt_new_ne_nblast, n_nb),
                  sprintf("%d / %d", cnt_st_ne_nblast, n_st_nb),
                  sprintf("%d / %d", cnt_st_ne_new, n_st_new)),
      stringsAsFactors = FALSE
    )
  }

  plot7_data <- agree_df %>%
    dplyr::group_split(dimorphism) %>%
    lapply(function(sub) {
      stats <- compute_agreement(sub)
      stats$dimorphism <- sub$dimorphism[1]
      stats$n_dimorphism <- nrow(sub)
      stats
    }) %>%
    dplyr::bind_rows()
  plot7_data$dimorphism_label <- sprintf("%s (n=%d)", plot7_data$dimorphism, plot7_data$n_dimorphism)

  plot_type_agreement_dimorphism <- ggplot(plot7_data,
                          aes(x = category, y = pct, fill = category)) +
    geom_col(width = 0.6, color = NA, show.legend = FALSE) +
    geom_text(aes(label = n_label), vjust = -0.5, size = 2.5) +
    scale_fill_manual(values = paper.cols) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.2))) +
    facet_wrap(~ dimorphism_label, scales = "free_x",
               labeller = labeller(.default = function(x) gsub("_", " ", x))) +
    labs(x = NULL, y = "Percentage") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  ggsave(file.path(plot.dir, "vnc_type_agreement_by_dimorphism.pdf"),
         plot = plot_type_agreement_dimorphism, width = 4, height = 3, bg = "white")
  message("  Saved: vnc_type_agreement_by_dimorphism.pdf")
}

###########################################################
### Panel C: NBLAST violin
###########################################################

if (nblast_available) {
  message("\n=== Panel C: NBLAST violin ===")

  plot_data <- dplyr::bind_rows(
    nblast_all_df %>% dplyr::filter(!is.na(top_nblast)) %>%
      dplyr::transmute(score = top_nblast, group = "nblast", dimorphism),
    nblast_all_df %>% dplyr::filter(!is.na(st_nblast)) %>%
      dplyr::transmute(score = st_nblast, group = "by morphology", dimorphism),
    nblast_all_df %>% dplyr::filter(!is.na(new_nblast)) %>%
      dplyr::transmute(score = new_nblast, group = "by connectivity", dimorphism)
  )
  plot_data$group <- factor(plot_data$group,
    levels = c("nblast", "by morphology", "by connectivity"))

  violin_data <- nblast_all_df %>%
    dplyr::filter(!is.na(top_nblast), !is.na(new_nblast)) %>%
    tidyr::pivot_longer(
      cols = c(top_nblast, st_nblast, new_nblast),
      names_to = "comparison",
      values_to = "score"
    ) %>%
    dplyr::filter(!is.na(score)) %>%
    dplyr::mutate(comparison = factor(
      dplyr::case_when(
        comparison == "top_nblast" ~ "nblast",
        comparison == "st_nblast" ~ "by morphology",
        comparison == "new_nblast" ~ "by connectivity"
      ),
      levels = c("nblast", "by morphology", "by connectivity")
    ))

  comparisons_list <- list(
    c("nblast", "by connectivity"),
    c("nblast", "by morphology"),
    c("by morphology", "by connectivity")
  )

  plot_nblast_violin <- ggplot(violin_data,
                                aes(x = comparison, y = score, fill = comparison)) +
    geom_violin(color = NA, alpha = 0.8, trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "grey40") +
    scale_fill_manual(values = c(
      "nblast" = unname(paper.cols["nblast"]),
      "by morphology" = unname(paper.cols["matched by morphology"]),
      "by connectivity" = unname(paper.cols["matched by connectivity"]))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    facet_wrap(~ dimorphism,
               labeller = labeller(.default = function(x) gsub("_", " ", x))) +
    ggpubr::stat_compare_means(method = "wilcox.test",
                                label = "p.signif",
                                comparisons = comparisons_list) +
    labs(x = NULL, y = "NBLAST score") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(plot.dir, "vnc_nblast_violin.pdf"),
         plot = plot_nblast_violin, width = 4, height = 3, bg = "white")
  message("  Saved: vnc_nblast_violin.pdf")
}

###########################################################
### Panel D: Reviewed NBLAST scatter
###########################################################

if (has_reviewed) {
  message("\n=== Panel D: Reviewed NBLAST scatter ===")

  scatter_df <- reviewed %>%
    dplyr::filter(!is.na(seatable_nblast), !is.na(new_nblast)) %>%
    dplyr::mutate(dimorphism = gsub("_", " ", dimorphism))

  plot_reviewed_scatter <- ggplot(scatter_df,
                           aes(x = seatable_nblast, y = new_nblast,
                               color = accept_new)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 1, stroke = 0.8, alpha = 0.8) +
    scale_color_manual(values = paper.cols) +
    labs(x = "Morphology match NBLAST score",
         y = "Connectivity match NBLAST score",
         color = "") +
    theme_minimal() +
    theme(legend.position = "top",
          panel.grid.minor = element_blank())

  ggsave(file.path(plot.dir, "vnc_reviewed_nblast_scatter.pdf"),
         plot = plot_reviewed_scatter, width = 4, height = 3, bg = "white")
  message("  Saved: vnc_reviewed_nblast_scatter.pdf")
}

###########################################################
### Panels E & F: NBLAST rank histograms
###########################################################

if (nblast_available) {
  message("\n=== Panels E & F: NBLAST rank histograms ===")

  # Non-sex-specific neurons with NBLAST data
  rank_neurons <- df %>%
    dplyr::filter(!dimorphism %in% c("female-specific", "male-specific")) %>%
    dplyr::filter(!is.na(connectivity_type), connectivity_type != "") %>%
    dplyr::filter(root_626 %in% unique(nb_vnc$root_626)) %>%
    dplyr::mutate(banc_id = root_626, effective_type = connectivity_type)

  # Rank computation helper
  compute_new_type_rank <- function(nb_data, neuron_df, id_col = "root_626") {
    rank_results <- lapply(seq_len(nrow(neuron_df)), function(i) {
      row <- neuron_df[i, ]
      hits <- nb_data %>%
        dplyr::filter(!!rlang::sym(id_col) == row$banc_id) %>%
        dplyr::arrange(dplyr::desc(score))
      if (nrow(hits) == 0) return(NULL)
      new_type_rows <- which(hits$match_cell_type == row$effective_type)
      if (length(new_type_rows) == 0) return(NULL)
      type_ranks <- hits %>%
        dplyr::group_by(match_cell_type) %>%
        dplyr::summarise(best_score = max(score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(best_score)) %>%
        dplyr::mutate(rank = dplyr::row_number())
      new_rank <- type_ranks$rank[type_ranks$match_cell_type == row$effective_type]
      data.frame(banc_id = row$banc_id, dimorphism = row$dimorphism,
                 rank = new_rank[1], stringsAsFactors = FALSE)
    })
    dplyr::bind_rows(rank_results)
  }

  # Rank plot helper
  make_rank_plot <- function(ranks_df) {
    if (nrow(ranks_df) == 0) return(NULL)
    binned <- ranks_df %>%
      dplyr::filter(!is.na(rank)) %>%
      dplyr::mutate(
        rank_bin = dplyr::if_else(rank >= 10, "10+", as.character(rank)),
        rank_bin = factor(rank_bin, levels = c(as.character(1:9), "10+")),
        dimorphism = factor(dimorphism, levels = c("isomorphic", "dimorphic"))
      )
    rank_summary <- binned %>%
      dplyr::count(rank_bin, dimorphism, .drop = FALSE) %>%
      dplyr::filter(!is.na(rank_bin))
    ggplot(rank_summary, aes(x = rank_bin, y = n, fill = dimorphism)) +
      geom_col(color = NA) +
      scale_fill_manual(values = paper.cols, drop = FALSE) +
      scale_x_discrete(drop = TRUE) +
      labs(x = "Rank of connectivity-matched type", y = "Count", fill = NULL) +
      theme_minimal() +
      theme(legend.position = "top",
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank())
  }

  # Side lookups
  rank_ids <- unique(rank_neurons$banc_id)
  rank_rootids <- id_map$root_id[id_map$banc_id %in% rank_ids]

  query_side_map <- rank_neurons %>%
    dplyr::distinct(banc_id, side)

  manc_side_map <- fm %>%
    dplyr::filter(!is.na(manc_id), !is.na(side), side != "") %>%
    dplyr::mutate(manc_id = as.character(manc_id)) %>%
    dplyr::distinct(manc_id, .keep_all = TRUE) %>%
    dplyr::select(manc_id, match_side = side)

  # --- Panel E: MANC ranks ---
  nb_manc_vnc <- nb_all %>%
    dplyr::filter((root_626 %in% rank_ids | pt_root_id %in% rank_rootids) &
                    score > 0) %>%
    dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score) %>%
    dplyr::mutate(
      root_626 = dplyr::case_when(
        root_626 %in% rank_ids ~ root_626,
        pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
        TRUE ~ root_626
      ),
      match_id = as.character(match_id)
    ) %>%
    dplyr::select(-pt_root_id) %>%
    dplyr::left_join(manc_side_map, by = c("match_id" = "manc_id")) %>%
    dplyr::left_join(query_side_map, by = c("root_626" = "banc_id")) %>%
    dplyr::filter(!is.na(match_side), !is.na(side), match_side == side) %>%
    dplyr::select(-match_side, -side, -match_id)

  manc_types <- unique(nb_manc_vnc$match_cell_type)
  manc_rank_neurons <- rank_neurons %>%
    dplyr::filter(effective_type %in% manc_types)
  message(sprintf("  MANC: %d/%d neurons with effective_type in MANC types",
                  nrow(manc_rank_neurons), nrow(rank_neurons)))

  manc_ranks <- compute_new_type_rank(nb_manc_vnc, manc_rank_neurons, "root_626")
  message(sprintf("  MANC ranks: %d neurons", nrow(manc_ranks)))

  plot_nblast_rank_manc <- make_rank_plot(manc_ranks)
  if (!is.null(plot_nblast_rank_manc)) {
    ggsave(file.path(plot.dir, "vnc_nblast_rank_manc.pdf"),
           plot = plot_nblast_rank_manc, width = 4, height = 3, bg = "white")
    message("  Saved: vnc_nblast_rank_manc.pdf")
  }

  # --- Panel F: maleCNS ranks ---
  malecns_nblast_path <- file.path(nblast.gcs.path, "banc_malecns_v0.9_nblast.feather")
  malecns_nblast_local <- file.path(nblast_cache_dir, basename(malecns_nblast_path))
  if (!file.exists(malecns_nblast_local)) {
    message("  Downloading maleCNS NBLAST feather from GCS...")
    system2("gsutil", c("cp", malecns_nblast_path, malecns_nblast_local),
            stdout = TRUE, stderr = TRUE)
  }

  if (file.exists(malecns_nblast_local)) {
    nb_malecns <- arrow::read_feather(malecns_nblast_local)
    message(sprintf("  maleCNS NBLAST: %d rows", nrow(nb_malecns)))

    # maleCNS metadata already loaded by startup.R
    malecns_meta <- malecns.meta %>%
      dplyr::select(dplyr::any_of(c("malecns_09_id", "cell_type", "manc_cell_type", "side")))
    malecns_type_map <- malecns_meta %>%
      dplyr::mutate(malecns_09_id = as.character(malecns_09_id)) %>%
      dplyr::distinct(malecns_09_id, .keep_all = TRUE) %>%
      dplyr::select(malecns_09_id, manc_cell_type, match_side = side)

    nb_malecns_vnc <- nb_malecns %>%
      dplyr::filter((root_626 %in% rank_ids | pt_root_id %in% rank_rootids) &
                      score > 0) %>%
      dplyr::select(root_626, pt_root_id, match_id, match_cell_type, score) %>%
      dplyr::mutate(
        root_626 = dplyr::case_when(
          root_626 %in% rank_ids ~ root_626,
          pt_root_id %in% names(rootid_to_bancid) ~ rootid_to_bancid[pt_root_id],
          TRUE ~ root_626
        ),
        match_id = as.character(match_id)
      ) %>%
      dplyr::select(-pt_root_id) %>%
      dplyr::left_join(malecns_type_map, by = c("match_id" = "malecns_09_id")) %>%
      dplyr::mutate(match_cell_type = dplyr::coalesce(manc_cell_type, match_cell_type)) %>%
      dplyr::left_join(query_side_map, by = c("root_626" = "banc_id")) %>%
      dplyr::filter(!is.na(match_side), !is.na(side), match_side == side) %>%
      dplyr::select(-manc_cell_type, -match_id, -match_side, -side)

    malecns_types <- unique(nb_malecns_vnc$match_cell_type)
    malecns_rank_neurons <- rank_neurons %>%
      dplyr::filter(effective_type %in% malecns_types)
    message(sprintf("  maleCNS: %d/%d neurons with effective_type in maleCNS types",
                    nrow(malecns_rank_neurons), nrow(rank_neurons)))

    malecns_ranks <- compute_new_type_rank(nb_malecns_vnc, malecns_rank_neurons, "root_626")
    message(sprintf("  maleCNS ranks: %d neurons", nrow(malecns_ranks)))

    plot_nblast_rank_malecns <- make_rank_plot(malecns_ranks)
    if (!is.null(plot_nblast_rank_malecns)) {
      ggsave(file.path(plot.dir, "vnc_nblast_rank_malecns.pdf"),
             plot = plot_nblast_rank_malecns, width = 4, height = 3, bg = "white")
      message("  Saved: vnc_nblast_rank_malecns.pdf")
    }
    rm(nb_malecns, nb_malecns_vnc, malecns_meta); gc()
  } else {
    message("  maleCNS NBLAST feather not available; skipping Panel F")
  }

  rm(nb_manc_vnc, manc_ranks); gc()
}

###########################################################
### Additional plots (not figure panels)
###########################################################

message("\n=== Additional plots ===")

# --- NBLAST density ---
if (nblast_available) {
  plot_data <- dplyr::bind_rows(
    nblast_all_df %>% dplyr::filter(!is.na(top_nblast)) %>%
      dplyr::transmute(score = top_nblast, group = "nblast", dimorphism),
    nblast_all_df %>% dplyr::filter(!is.na(st_nblast)) %>%
      dplyr::transmute(score = st_nblast, group = "by morphology", dimorphism),
    nblast_all_df %>% dplyr::filter(!is.na(new_nblast)) %>%
      dplyr::transmute(score = new_nblast, group = "by connectivity", dimorphism)
  )
  plot_data$group <- factor(plot_data$group,
    levels = c("nblast", "by morphology", "by connectivity"))

  plot_nblast_density <- ggplot(plot_data, aes(x = score, fill = group)) +
    geom_density(alpha = 0.4, color = NA) +
    scale_fill_manual(values = paper.cols) +
    facet_wrap(~ dimorphism) +
    labs(x = "NBLAST score", y = "Density", fill = NULL) +
    theme_minimal() +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  ggsave(file.path(plot.dir, "vnc_nblast_density.pdf"),
         plot = plot_nblast_density, width = 12, height = 5, dpi = 300, bg = "white")
  message("  Saved: vnc_nblast_density.pdf")
}

# --- Type agreement (unfaceted) ---
if (nblast_available) {
  n_total <- nrow(agree_df)
  n_with_nblast <- sum(!is.na(agree_df$top_nblast_type))
  plot6_data <- compute_agreement(agree_df)

  plot_type_agreement <- ggplot(plot6_data,
                          aes(x = category, y = pct, fill = category)) +
    geom_col(width = 0.6, color = NA, show.legend = FALSE) +
    geom_text(aes(label = n_label), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = paper.cols) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(subtitle = sprintf("n=%d neurons, %d with MANC NBLAST data",
                            n_total, n_with_nblast),
         x = NULL, y = "Percentage") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(plot.dir, "vnc_type_agreement.pdf"),
         plot = plot_type_agreement, width = 7, height = 5, dpi = 300, bg = "white")
  message("  Saved: vnc_type_agreement.pdf")
}

# --- Type changes by dimorphism ---
df1_summary <- df %>%
  dplyr::count(dimorphism, type_changed) %>%
  dplyr::group_by(dimorphism) %>%
  dplyr::mutate(prop = n / sum(n), total = sum(n)) %>%
  dplyr::ungroup()
df1_totals <- df1_summary %>% dplyr::distinct(dimorphism, total)

plot_type_changes_dimorphism <- ggplot(df1_summary,
                     aes(x = dimorphism, y = prop, fill = type_changed)) +
  geom_col(color = NA) +
  geom_text(data = df1_totals,
            aes(x = dimorphism, y = 1.03, label = total, fill = NULL),
            size = 3.5) +
  scale_fill_manual(values = paper.cols) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Proportion", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

ggsave(file.path(plot.dir, "vnc_type_changes_by_dimorphism.pdf"),
       plot = plot_type_changes_dimorphism, width = 7, height = 5, dpi = 300, bg = "white")
message("  Saved: vnc_type_changes_by_dimorphism.pdf")

message("\n### Figure panels complete ###")
