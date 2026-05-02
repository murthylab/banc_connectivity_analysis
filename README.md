# BANC connectivity analysis paper resources

Repository for data, code, figures, and manuscript assets for
["Effective Connectome Alignment For Cell Typing and Identifying Sex
Differences."](https://docs.google.com/document/d/1MUOX8YmrFuuWjsmmGci5Kq9HxQTjkg0KoB2Pt4m_6Xo/edit?tab=t.0)

## Repository structure

```
R/
  startup.R                                  — data loading from GCS + SeaTable
  panels-vnc-dimorphic-influence.R           — sensory → dimorphic influence heatmaps
  panels-vnc-morphology-connectivity-types.R — VNC type agreement & NBLAST panels
  banc-dimorphic-density.R                   — synapse density maps (VNC + brain)
  banc-vnc-type-changes.R                    — VNC type change CSVs + summary
  banc-vnc-type-review-images.R              — comparison PNGs for type-change review
settings/
  paper_colours_lacroix.csv                  — shared colour palette
data/
  banc/                                      — BANC-specific data & influence cache
  manc/, mcns/                               — cell inventory CSVs per dataset
  codex/                                     — drop-in inputs for type-changes script
  cache/                                     — GCS download cache (auto-created)
figs/
  figure_typing/                             — cell typing panels (.ai + linked PDFs)
  figure_dimorphic/                          — sexual dimorphism panels + density maps
    links/                                   — generated PDFs/PNGs
```

## Analysis scripts

### `panels-vnc-dimorphic-influence.R`

Computes sensory → dimorphic/sex-specific neuron influence for BANC (female)
and maleCNS (male), producing paired heatmaps across three super_class groups
(VNC intrinsic, ascending, descending). Influence is recalculated from scratch
each run using the `influencer` package (PETSc solver). Outputs go to
`figs/figure_dimorphic/links/`.

### `panels-vnc-morphology-connectivity-types.R`

Compares VNC cell type assignments between morphology-based (NBLAST) and
connectivity-based approaches, stratified by sexual dimorphism. Uses NBLAST
scores from GCS and reviewed type-change CSVs from `data/`. Outputs go to
`figs/figure_typing/links/`.

### `banc-dimorphic-density.R`

Density plots of dimorphic / sex-specific synapses across BANC, MANC, FAFB,
and maleCNS, registered to JRCVNC2018U (VNC) and JRC2018U (brain) reference
spaces. Produces faceted density grids and BANC-vs-other subtraction maps for
both pre- and post-synapses. Requires `gsutil`, `nat.flybrains`,
`nat.templatebrains`, `nat.jrcbrains`, and the `navis` Python package
(via reticulate) for maleCNS↔BANC coordinate transforms. Outputs go to
`figs/figure_dimorphic/links/`. Adapted from
`bancpipeline/banc/update/banc-dimorphic-density.R`; runs standalone.

### `banc-vnc-type-changes.R`

Reads BANC seatable + reviewed CSVs, joins to `franken_meta`, identifies VNC
intrinsic + effector type changes, computes hemilineage changes, looks up
NBLAST scores per type-change category, and writes:

- `data/vnc_type_changes.csv`
- `data/vnc_hemilineage_changes.csv`
- `data/vnc_effector_type_changes.csv`
- `data/mcns_effector_type_changes.csv`

Optional input CSVs (BANC dimorphism, effector matches, LR mirror matches,
alignment scores) live under `data/codex/`; sections that depend on them are
skipped if files are missing. SeaTable update blocks from the original
bancpipeline script have been intentionally omitted — re-add and uncomment
`banctable_update_rows()` calls if you need to push back. Adapted from
`bancpipeline/banc/update/banc-vnc-type-changes.R`; runs standalone.

### `banc-vnc-type-review-images.R`

For each entry in `data/vnc_type_changes_reveiwed.csv` with `accept_new = NA`,
generates a 3-panel PNG comparing the BANC neuron against the best MANC
NBLAST match for both `seatable_type` and `effective_type`. Uses the
BANC↔MANC NBLAST feather (GCS) for best-match selection with `franken_meta`
fallback. Output PNGs go to `figs/figure_dimorphic/links/`. Adapted from
`bancpipeline/banc/update/banc-vnc-type-review-images.R`; runs standalone.

## Data sources

All connectome data is loaded from Google Cloud Storage at runtime:

- **BANC metadata & edgelist**: `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/`
- **maleCNS metadata & edgelist**: `gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/malecns_09/`
- **NBLAST scores**: `gs://lee-lab_brain-and-nerve-cord-fly-connectome/nblast/`
- **SeaTable** (BANC + maleCNS): queried via `bancr::banctable_query()` for
  columns not in the GCS feathers (e.g. `sexually_dimorphic`, `malecns_cell_type`)

Downloaded files are cached locally in `data/cache/`. Delete this directory to
force a fresh download.

## Dependencies

- **R packages**: `bancr`, `influencer`, `malevnc`, `dplyr`, `arrow`,
  `reticulate`, `ComplexHeatmap`, `circlize`, `ggplot2`, `ggpubr`,
  `nat.flybrains`, `nat.templatebrains`, `nat.jrcbrains`, `nat.ggplot`,
  `Rvcg`, `rgl`, `fafbseg`
- **Python** (via reticulate): `petsc4py`, `gcsfs`, `pyarrow`,
  `navis` + `flybrains` (for maleCNS↔BANC coordinate transforms)
- **CLI**: `gsutil` (for NBLAST + synapse feather/parquet downloads)

## Running

```bash
# From the repository root:
Rscript R/panels-vnc-dimorphic-influence.R
Rscript R/panels-vnc-morphology-connectivity-types.R
Rscript R/banc-dimorphic-density.R
Rscript R/banc-vnc-type-changes.R
Rscript R/banc-vnc-type-review-images.R
```

Set `BANC_VERSION` to override the default data version (e.g. `BANC_VERSION=banc_900`).
