# BANC connectivity analysis paper resources

Repository for data, code, figures, and manuscript assets for
["Effective Connectome Alignment For Cell Typing and Identifying Sex
Differences."](https://docs.google.com/document/d/1MUOX8YmrFuuWjsmmGci5Kq9HxQTjkg0KoB2Pt4m_6Xo/edit?tab=t.0)

## Repository structure

```
R/
  startup.R                                 — data loading from GCS + SeaTable
  panels-vnc-dimorphic-influence.R          — sensory → dimorphic influence heatmaps
  panels-vnc-morphology-connectivity-types.R — VNC type agreement & NBLAST panels
settings/
  paper_colours_lacroix.csv                 — shared colour palette
data/
  banc/                                     — BANC-specific data & influence cache
  manc/, mcns/                              — cell inventory CSVs per dataset
  cache/                                    — GCS download cache (auto-created)
figs/
  figure_typing/                            — cell typing panels (.ai + linked PDFs)
  figure_dimorphic/                         — sexual dimorphism panels + density maps
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

- **R packages**: `bancr`, `influencer`, `dplyr`, `arrow`, `reticulate`,
  `ComplexHeatmap`, `circlize`, `ggplot2`, `ggpubr`
- **Python** (via reticulate): `petsc4py`, `gcsfs`, `pyarrow`
- **CLI**: `gsutil` (for NBLAST feather downloads)

## Running

```bash
# From the repository root:
Rscript R/panels-vnc-dimorphic-influence.R
Rscript R/panels-vnc-morphology-connectivity-types.R
```

Set `BANC_VERSION` to override the default data version (e.g. `BANC_VERSION=banc_900`).
