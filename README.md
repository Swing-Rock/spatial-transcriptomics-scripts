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
- **reference dataset:** Puram et al. 2017 HNSCC scRNA-seq dataset (GEO: [GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322)); downloaded and cached automatically
  - **rationale:** this is a landmark single-cell RNA-seq atlas of head and neck squamous cell carcinoma (HNSCC) profiling ~6,000 cells from 18 patients, with expert-annotated cell types including Malignant, T cell, B cell, Fibroblast, Macrophage, Endothelial, Dendritic, and Mast cells. Using a disease-matched HNSCC reference is critical because RCTD learns cell-type expression signatures from the reference — a generic or non-HNSCC reference would fail to accurately distinguish tumor cells from the stromal and immune compartments that are characteristic of this cancer type
- `get_hnscc_reference()`: downloads, parses, and caches the Puram 2017 reference as a spacexr `Reference` object; note that GSE103322 stores log2(TPM+1) values which are back-transformed to integer pseudo-counts for compatibility with spacexr
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
- entry point for the Xenium pipeline; set all path variables and operation flags here
- sources `xenium_qc.R`, `xenium_annotation.R`, and `xenium_merge.R`
- `batch_processing`: if TRUE, loops over all sample subfolders in `data_folder` (excludes the output folder automatically); if FALSE, runs `data_folder` as a single sample
- `single_cell_annotation`: if TRUE, labels each cell individually with SingleR; if FALSE, clusters first then labels each cluster
- orchestrates the full pipeline: QC → annotation → merge, with `.rds` caching at each stage
- also contains `hpv_analysis()`, which runs HPV pathway analysis on the merged Xenium object (same steps as the Visium HPV pipeline: DE → GO → KEGG → spatial scoring)
- includes a run-time timer that prints total elapsed time on completion

### xenium_qc.R
- `run_xenium_qc(data_path, data_name, graph_output, out_path)`: QC and normalization for a single Xenium sample
- loads Xenium data using Seurat's `LoadXenium`; the raw assay is named `"Xenium"` (not `"Spatial"`)
- computes QC metrics with special handling for targeted panels that may lack MT or ribosomal genes:
  - `percent.mt`: mitochondrial content (MT- prefix); set to 0 if no MT features in the panel
  - `percent.ribo`: ribosomal content (RPS/RPL prefix); set to 0 if absent
  - `blank_rate`: percentage of blank/negative-control codeword counts (BLANK_ or NegControlProbe_ prefix); set to 0 if absent
- filters cells using thresholds tuned for targeted Xenium panels (typically 100–500 genes):
  - nFeature_Xenium > 5
  - nCount_Xenium > 10
  - percent.mt < 20%
  - blank_rate < 5%
- normalizes with SCTransform (stored in the `"SCT"` assay)
- saves a combined violin plot comparing raw / filtered / normalized distributions as a PNG

### xenium_annotation.R
- `run_xenium_annotation(data_obj, data_name, out_path, single_cell)`: cluster and annotate a single Xenium Seurat object
- because Xenium is single-cell resolution (not spot-level), annotation is performed at the individual cell level rather than the spot level
- runs PCA on the SCT assay, then UMAP (dims 1:30)
- annotates using SingleR with the **Human Primary Cell Atlas** reference (`celldex::HumanPrimaryCellAtlasData()`), mirroring the annotation strategy in `st_clustering.R`
- **single cell mode** (`single_cell = TRUE`): labels each cell individually; saves annotated UMAP and `ImageDimPlot` spatial plots
- **cluster mode** (`single_cell = FALSE`): runs FindNeighbors/FindClusters (resolution = 0.5) first, then assigns one SingleR label per cluster; saves cluster UMAP and spatial plots
- all labels stored in `SingleR_label` metadata column

### xenium_merge.R
- `merge_xenium(data_obj, out_path)`: merges a named list of annotated Xenium Seurat objects and integrates across samples
- re-runs SCTransform on the merged object to produce a shared feature space across samples (necessary because per-sample SCT residuals are not directly comparable after merging)
- runs the standard pre-integration workflow: FindVariableFeatures (3,000 features, VST) → PCA → FindNeighbors → FindClusters → UMAP; embeddings stored as `pca.SCT` and `umap.SCT`
- applies Harmony batch correction (`group.by = orig.ident`, theta = 8) on the SCT PCA; corrected embedding stored as `harmony.SCT`
- re-clusters on the Harmony embedding; final UMAP stored as `umap.harmony.SCT`
- saves a three-panel PNG: pre-Harmony UMAP by sample, post-Harmony UMAP by sample, and post-Harmony UMAP by cell type

### tcga_survival.R
- associates HPV pathway gene sets (derived from spatial GO enrichment in `st_hpv_pathway_analysis.R`) with overall survival outcomes in TCGA bulk RNA-seq data
- **reference datasets:** two clinical data sources are supported (see below); expression data is always from TCGA-HNSC via `TCGAbiolinks`
  - **TCGA-HNSC expression (TCGAbiolinks):** the largest publicly available, clinically annotated bulk RNA-seq cohort for HNSCC (~500 patients). Using TCGA allows us to ask whether spatially-derived HPV pathway signatures have prognostic relevance at the population level — bridging single-sample spatial findings to a generalizable clinical context
  - **Pre-processed clinical CSV — `survival_data_HNSC.csv` (recommended):** a pre-processed TCGA-HNSC clinical table (`GSE281978/output/survival_data_HNSC.csv`) with 489 patients, pre-computed overall survival times (`overall_survival`), a clean boolean event indicator (`deceased`), and HPV status directly encoded in the `Subtype` column as `HNSC_HPV+` / `HNSC_HPV-`. This is preferred over the raw TCGAbiolinks clinical download because it avoids the ambiguous multi-column HPV parsing in `tidy_clinical()` and is validated to match all 489 expression matrix samples. Pass its path via the `clinical_csv` argument of `run_tcga_survival()`
  - **TCGAbiolinks clinical fallback:** if `clinical_csv` is not provided (or the path is invalid), the pipeline falls back to downloading clinical data via `TCGAbiolinks` and parsing it with `tidy_clinical()`; HPV status is sourced from `TCGAquery_subtype("HNSC")` (`paper_HPV.status` field)
- `load_external_clinical()`: loads and standardizes the pre-processed CSV; handles the one known NA in `overall_survival` (patient `TCGA-H7-A6C4`, where `days_to_death` was stored as `'--` in the source) by falling back to `diagnoses.days_to_last_follow_up`; outputs the same `submitter_id / os_time / os_event / hpv_status` format as `tidy_clinical()`
- `get_tcga_hnsc()`: downloads TCGA-HNSC STAR counts and clinical data via `TCGAbiolinks`; also downloads HPV subtype annotations; caches everything as a single `.rds`
- `tidy_clinical()`: extracts and cleans overall survival time, event status, and HPV status from the raw TCGAbiolinks clinical table; standardizes HPV labels to `"positive"` / `"negative"`; falls back to the subtype table for HPV if not present in the clinical table
- `build_expr_matrix()`: extracts log2(TPM+1) expression (or log2(CPM+1) if TPM unavailable); deduplicates to one sample per patient using 12-character TCGA barcodes
- `extract_go_genesets()`: pulls the top N significant GO terms (default 5) from each HPV direction (HPV+ upregulated, HPV− upregulated) from the spatial `run_go()` output
- `run_gsva()`: scores extracted gene sets in TCGA bulk expression using GSVA; filters sets to genes present in the expression matrix
- `run_km_survival()`: dichotomizes each pathway score at a configurable quantile (default: median split) and produces Kaplan-Meier survival curves; optionally stratifies by HPV status; saves PNGs via `ggsurvplot`
- `run_cox_survival()`: runs univariate Cox proportional hazards regression for each pathway score and one joint multivariate Cox model (adjusting for HPV status and age where available); saves a forest plot and results CSVs
- `run_tcga_survival()`: main entry point; call after `hpv_analysis()` in `st_master.R`, passing the `go_results` object; accepts optional `clinical_csv` argument — pass the path to `survival_data_HNSC.csv` to use the pre-processed data

---

## data_files_info.md
- reference document describing the files produced by 10x Genomics Visium and Xenium platforms
- covers the purpose and format of key files: `filtered_feature_bc_matrix.h5`, `tissue_positions_list.csv`, `scalefactors_json.json`, `tissue_hires/lowres_image.png` (Visium) and `cell_feature_matrix.h5`, `cells.parquet`, `cell_boundaries.parquet`, `transcripts.parquet`, `morphology.ome.tif` (Xenium)

---

## file organization

### Visium
Make sure each sample's file structure is organized as follows:
- (Folder) Sample name
  - (File) `(SampleName)_filtered_feature_bc_matrix.h5`
  - (Folder) spatial
    - (File) `tissue_lowres_image.png`
    - (File) `scalefactors_json.json`
    - (File) `tissue_positions_list.csv`

### Xenium
Each sample folder should contain the standard Xenium output bundle:
- (Folder) Sample name
  - (File) `cell_feature_matrix.h5`
  - (File) `cells.csv.gz` (or `cells.parquet`)
  - (File) `experiment.xenium`
  - (File) `morphology.ome.tif` (optional)

In batch mode, all sample folders should live inside a shared `data_folder`; an `output/` subfolder is automatically excluded from processing.

---

lots of thanks and credits to Alfred and Riya for providing the code for workflow and many helps
