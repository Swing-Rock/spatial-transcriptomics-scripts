# CRC Spatial Transcriptomics Pipeline

An end-to-end pipeline for processing spatial transcriptomics data (10x Visium / CytAssist) from colorectal cancer (CRC) and liver metastasis (CRLM) samples — spanning ingestion of heterogeneous GEO data formats, QC, signature scoring, and cell type deconvolution via `cell2location`.

The pipeline is split across **R** (data ingestion, QC, format standardization, signature scoring) and **Python** (deconvolution), connected by a common on-disk data layout modeled on 10x Space Ranger's output structure.

---

## Pipeline Overview

```
                     ┌─────────────────────────┐
  Raw GEO downloads  │   reorganize samples    │  
  (.mtx/.h5/.h5ad)    →                        |
                     └─────────────────────────┘
                                 │
                                 ▼
                     ┌─────────────────────────┐
                     │     st_data_qc.R        │  Auto-detect format (MTX /
                     │  (run_qc, format        →   h5 / h5ad), load, QC filter,
                     │   detection, loaders)    │   SCTransform normalize
                     └─────────────────────────┘
                                 │
                                 ▼
                     ┌─────────────────────────┐
                     │  CRC_analysis.R         │  Batch QC, pseudobulking,
                     │                         │  EdgeR CPM, MAP signature
                     │                         │  scoring (MSS/MSI etc.)
                     └─────────────────────────┘
                                 │
                                 ▼
                     ┌─────────────────────────┐
                     │    Cell2location.py     │  Load Visium (h5/h5ad),
                     │  (Python / scanpy /      │  filter MT genes, match
                     │   cell2location)         │  reference signatures,
                     │                          │  train, export abundances
                     └─────────────────────────┘
                                 │
                                 ▼
                 Per-sample cell type abundance CSVs
                 + deconvolved .h5ad + trained model
```

---

## Repository Structure

```
CRC_project/
├── CRC_analysis.R  # Batch pipeline: QC → pseudobulk → EdgeR → MAP
├── Cell2location.py              # Python deconvolution stage
├── inf_aver.csv                  # Reference cell type signature matrix (for cell2location)
├── <GSE_ID>/                     # One folder per GEO series
│   ├── <GSE_ID>_features.tsv.gz  # Shared gene annotation (Ensembl ID + symbol) (if applicable)
│   ├── <GSE_ID>_barcodes.tsv.gz  # Shared barcode list (if applicable)
│   └── <sample_name>/            # One folder per sample (e.g. GSM8655157)
│       ├── filtered_feature_bc_matrix.h5             # or matrix.mtx
│       └── spatial/
│           ├── tissue_positions_list.csv.gz
│           ├── scalefactors_json.json.gz
│           ├── tissue_hires_image.png.gz
└──         └── tissue_lowres_image.png.gz

visium_scripts/
├── st_data_qc.R              # Core QC + multi-format loader (MTX/h5/h5ad)
├── st_helper.R               # Plotting helpers (make_v_plot, make_sf_plot, etc.)
├── split_dual_tissue.R       # Splits dual-tissue capture areas into 2 samples
└── export_rda_to_spaceranger.R  # Converts legacy .rda Seurat objects to SR layout
```

---

## Supported Input Formats

The pipeline was originally built around standard Space Ranger `.h5` output, then extended to handle two additional formats commonly found in GEO deposits:

| Format | How it's detected | Loader |
|---|---|---|
| Space Ranger `.h5` | `filtered_feature_bc_matrix.h5` present | `Load10X_Spatial()` |
| `.h5ad` | `*.h5ad` present | `load_h5ad_visium()` |
| Raw MTX + coordinates | `*.mtx(.gz)` present | `load_mtx_spatial()` |

All formats converge on the same `run_qc()` output (a QC'd, normalized Seurat object), so everything downstream is format-agnostic. See `MTX_Spatial_Data_Challenges.md` for the detailed engineering notes on the MTX-specific loader (gene symbol resolution, duplicate aggregation, image registration caveats, etc.).

---

## How to Run the Full Pipeline

### Step 1 — Organize raw downloads

### Step 2 — Run QC, pseudobulking, and MAP signature scoring

```r
source("visium_scripts/CRC_analysis_mtx_fixed.R")
```

This is the main R entry point. For each sample folder found in `data_folder`:
1. `run_qc()` auto-detects the input format and loads the sample
2. QC filtering is applied (`nFeature_Spatial`, `nCount_Spatial`, `percent.mt`, `percent.ribo` thresholds)
3. `SCTransform` normalization (if `run_qc()` parameter _normalization_ is true, but for this project we don't want to normalize cuz MAP takes raw counts)
4. QC diagnostic plots saved to `output/<sample>_qc_plots.png`
5. Pseudobulk aggregation (`AggregateExpression`) across all spots per sample
6. `edgeR` CPM normalization with `filterByExpr` (min.count = 50)
7. MAP signature scoring (patched to use the current `GSVA::ssgseaParam` API) — computes consensus molecular subtype-like scores; missing marker genes are imputed via per-sample median as a fallback

Outputs land in `output/cache/<sample_name>/`:
- `MAP_input.txt` — CPM matrix used as MAP input
- `ss_ssgsea.txt`, `ssgsea_score1-10.txt` — signature scores

Set `batch_processing <- TRUE` at the top of the script to process every sample folder in `data_folder`, or `FALSE` to run a single sample.

### Step 3 — Cell type deconvolution (Python)

Once samples are QC'd and in a consistent Space Ranger-style layout, run the deconvolution stage:

```bash
python Cell2location.py
```

**Before running**, set these paths at the top of the script:
```python
INPUT_FOLDER  = "CRC_project/GSE283052"      # parent folder of sample subdirectories
INF_AVER_PATH = "CRC_project/inf_aver.csv"   # reference cell type signature matrix
```

For each valid sample folder (must contain either `filtered_feature_bc_matrix.h5` or a `.h5ad`), the script will:
1. Load the Visium data via `squidpy.read.visium()` (h5) or `load_h5ad_visium()` (h5ad, reattaching images/scalefactors from `spatial/`)
2. Filter out mitochondrial genes
3. Intersect genes with the reference signature matrix (`inf_aver.csv`)
4. Train a `cell2location.models.Cell2location` model (default: 30,000 epochs, GPU if available)
5. Export the posterior and save:
   - `<sample>_deconvolved.h5ad` — full AnnData with abundance estimates
   - `<sample>_abundance.csv` — mean cell type abundance per spot
   - `<sample>_conservative_abundance.csv` — 5th-percentile (conservative) abundance per spot
   - `model/` — saved trained model

All outputs are written to `CRC_project/GSE283052/output/deconvolved/<sample_name>/`.

**Note:** the `filtered_feature_bc_matrix.h5` used at this stage can either be a genuine Space Ranger export or one produced by reverse-engineering an H5 from an MTX-derived Seurat object (see `MTX_Spatial_Data_Challenges.md`, Challenge 9) — both are read identically by `squidpy.read.visium()`.

---

## Key Configuration Points

| Script | Variable | Purpose |
|---|---|---|
| `st_data_qc.R` | `data_folder` | Parent folder of QC-ready sample subfolders |
| `st_data_qc.R` (`run_qc`) | QC thresholds | `nFeature_Spatial`, `nCount_Spatial`, `percent.mt`, `percent.ribo` cutoffs |
| `Cell2location.py` | `INPUT_FOLDER`, `INF_AVER_PATH` | Sample parent folder and reference signature path |
| `Cell2location.py` | `N_CELLS_PER_LOCATION`, `DETECTION_ALPHA`, `MAX_EPOCHS`, `N_POSTERIOR_SAMPLES` | Model hyperparameters |

---

## Known Limitations / Outstanding Items

- **Image registration for MTX-derived samples** is currently disabled in `load_mtx_spatial()` — histology-overlay spatial feature plots (`SpatialFeaturePlot`) are unavailable for these samples until a proper `VisiumV1`/`CreateImage()` object is constructed. QC violin plots and all downstream numerical analysis are unaffected.
- **Duplicate gene symbols** (multiple Ensembl IDs mapping to one symbol) are resolved by summing expression via sparse matrix multiplication rather than dropping duplicates 
- **`recover_raw_counts()`** in `Cell2location.py` is currently commented out in the main loop; re-enable it if a given sample's `adata.X` may contain normalized (non-integer) values rather than raw counts, since `cell2location` requires raw counts.
