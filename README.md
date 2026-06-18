# spatial-transcriptomics-scripts
Ongkeko Lab spatial transcriptomics data processing pipeline

## flow chart of script dependencies
[Flowchart.pdf](https://github.com/user-attachments/files/26293063/Flowchart.pdf)
https://lucid.app/lucidchart/9ca23dbf-cc43-4b66-b452-ce79e9b11e3e/edit?viewport_loc=-547%2C-325%2C3078%2C1476%2C0_0&invitationId=inv_80ce0e20-4ee8-4184-8f4a-6fc41391dc04

---

## visium_scripts/

### st_master.R
- entry point for the Visium pipeline; set all path variables and operation flags here
- `batch_processing`: if TRUE, loops over all sample subfolders in `data_folder`; if FALSE, runs a single sample
- `deconvolution`: if TRUE, runs RCTD deconvolution (st_deconvolution.R) instead of SingleR annotation
- `single_cell_annotation`: only relevant when `deconvolution = FALSE`; toggles between spot-level vs. cluster-level SingleR labeling
- sources all other visium scripts and orchestrates the full pipeline: QC → deconvolution/annotation → merge → HPV analysis
- caches intermediate results (`.rds` files) to avoid re-running expensive steps

### st_data_qc.R
- loads Visium `.h5` data using Seurat's `Load10X_Spatial`
- filters spots by the following thresholds:
  - nFeature_Spatial < 200 or > 7,500
  - nCount_Spatial < 250 or > 50,000
  - percent mitochondrial genes > 15%
  - percent ribosomal genes > 40%
- normalizes data using SCTransform

### st_clustering.R
- runs dimensionality reduction (PCA → UMAP) and cell type annotation using SingleR with the Human Primary Cell Atlas as reference
- **single cell labeling** (`single_cell = TRUE`): labels each spot individually before clustering
- **cluster labeling** (`single_cell = FALSE`): clusters spots first (FindNeighbors dims = 1:30, FindClusters resolution = 0.5), then assigns a cell type label to each cluster
- saves annotated UMAP and spatial plots as PNGs

### st_merge_data.R
- merges multiple per-sample Seurat objects into one combined object
- runs standard workflow: FindVariableFeatures → PCA → FindNeighbors → FindClusters → UMAP
- runs Harmony batch correction (`group.by = orig.ident`, theta = 8) and re-clusters on the corrected embedding
- saves a side-by-side UMAP plot comparing pre- and post-Harmony clustering

### st_deconvolution.R
- performs cell-type proportion deconvolution of Visium spots using RCTD (spacexr package)
- reference: Puram et al. 2017 HNSCC scRNA-seq dataset (GEO: GSE103322); downloaded and cached automatically
- `get_hnscc_reference()`: downloads, parses, and caches the Puram 2017 reference as a spacexr `Reference` object
- `run_rctd()`: runs RCTD in `doublet_mode = "full"` to estimate proportions for all cell types per spot; adds per-cell-type proportion columns and a `RCTD_dominant` column to Seurat metadata; saves spatial plots of dominant cell type and per-cell-type proportion heatmaps

### st_hpv_pathway_analysis.R
- differential expression and pathway enrichment analysis comparing HPV+ vs HPV− malignant spots
- `add_hpv_status()`: tags each spot with HPV status based on sample ID
- `subset_malignant()`: isolates malignant spots using either RCTD proportion threshold (`RCTD_Malignant >= 0.5`) or SingleR label matching
- `run_de()`: finds differentially expressed genes (Wilcoxon test, FDR < 0.05, log2FC threshold = 0.5); saves results CSV and a volcano plot
- `run_go()`: GO biological process enrichment (clusterProfiler) for HPV+ and HPV− upregulated gene sets; saves dotplots and result CSVs
- `run_kegg()`: KEGG pathway enrichment for HPV+ and HPV− gene sets; saves dotplots and result CSVs
- `run_spatial_scoring()`: scores top GO gene sets across all spots using `AddModuleScore` and saves spatial heatmap PNGs

### st_helper.R
- shared helper functions for plot formatting used across other scripts
- includes `title_plot()` for generating title panels and `make_v_plot()` for violin plots

---

## scripts/

### xenium_pipeline.R
- early-stage pipeline for processing 10x Xenium spatial transcriptomics data (work in progress)
- loads Xenium cell feature matrix (`.h5`), creates a Seurat object, and runs normalization → PCA → clustering → UMAP
- annotates clusters using SingleR with the Human Primary Cell Atlas reference
- contains commented-out code for loading spatial cell coordinates from `cells.parquet` and adding them as a spatial embedding

---

## data_files_info.md
- reference document describing the files produced by 10x Genomics Visium and Xenium platforms
- covers the purpose and format of key files: `filtered_feature_bc_matrix.h5`, `tissue_positions_list.csv`, `scalefactors_json.json`, `tissue_hires/lowres_image.png` (Visium) and `cell_feature_matrix.h5`, `cells.parquet`, `cell_boundaries.parquet`, `transcripts.parquet`, `morphology.ome.tif` (Xenium)

---

## file organization
Make sure each sample's file structure is organized as follows:
- (Folder) Sample name 
  - (File) (SampleName)_filtered_feature_bc_matrix.h5
  - (Folder) spatial
    - (File) tissue_lowres_image.png
    - (File) scalefactors_json.json
    - (File) tissue_positions_list.csv


lots of thanks and credits to Alfred and Riya for providing the code for workflow and many helps
