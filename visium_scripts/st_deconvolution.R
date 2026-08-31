# st_deconvolution.R
# Cell-type proportion deconvolution for Visium spots using RCTD (spacexr).
# Reference: Puram et al. 2017 HNSCC scRNA-seq dataset (GEO: GSE103322).
#
# Exported functions:
#   get_hnscc_reference(cache_path)           — download/cache the scRNA-seq reference
#   run_rctd(seurat_obj, reference, ...)      — deconvolve one sample, returns updated Seurat obj

library(Matrix)
library(Seurat)
library(spacexr)
library(ggplot2)

# ── Reference: Puram et al. 2017 (GSE103322) ──────────────────────────────────
# The data matrix has 5 metadata rows before gene counts begin:
#   Row 1: "processed by Maxima enzyme"   — 0/1 flag
#   Row 2: "Lymph node"                   — 0/1 flag
#   Row 3: "classified  as cancer cell"   — 1 = malignant, 0 = non-malignant
#   Row 4: "classified as non-cancer cells" — 1 = non-cancer, 0 = cancer (inverse of row 3)
#   Row 5: "non-cancer cell type"         — string label (e.g. "T cell", "B cell",
#                                           "Fibroblast", "Macrophage", "Endothelial",
#                                           "Dendritic", "Mast"); 0 for malignant cells
N_META_ROWS <- 5


#' Download, parse, and cache the Puram 2017 HNSCC scRNA-seq reference.
#'
#' @param cache_path  Directory to store the downloaded file and .rds cache.
#' @return            A spacexr::Reference object.
get_hnscc_reference <- function(cache_path) {
  ref_rds <- file.path(cache_path, "puram2017_reference.rds")

  if (file.exists(ref_rds)) {
    cat("Loading cached Puram 2017 HNSCC reference...\n")
    return(readRDS(ref_rds))
  }

  cat("Downloading Puram et al. 2017 HNSCC scRNA-seq reference (GSE103322)...\n")
  gz_file <- file.path(cache_path, "GSE103322_HNSCC_data_matrix.txt.gz")

  if (!file.exists(gz_file)) {
    download.file(
      url      = paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz"),
      destfile = gz_file,
      mode     = "wb",
      timeout  = 600
    )
    cat("Download complete.\n")
  }

  cat("Parsing reference matrix (may take a minute)...\n")
  raw <- read.table(
    gzfile(gz_file), sep = "\t", header = TRUE,
    row.names = 1, check.names = FALSE, comment.char = ""
  )

  # Print metadata rows so the user can sanity-check
  cat("\n--- Metadata rows (first", N_META_ROWS, "rows) ---\n")
  print(rownames(raw)[seq_len(N_META_ROWS)])
  cat("---\n\n")

  meta <- raw[seq_len(N_META_ROWS), ]

  # ── Identify malignancy row ────────────────────────────────────────────────
  # Row 3: "classified  as cancer cell"; encoding: 1 = malignant, 0 = non-malignant
  malignant_row <- grep("classified  as cancer cell", rownames(meta), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(malignant_row)) stop("Could not find malignancy row in metadata. Check N_META_ROWS.")
  malignant_vals <- as.character(meta[malignant_row, ])

  # ── Identify non-malignant cell type row ──────────────────────────────────
  # Row 5: "non-cancer cell type"; already contains plain-text labels
  nm_row <- grep("non-cancer cell type", rownames(meta), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(nm_row)) stop("Could not find non-malignant cell type row in metadata. Check N_META_ROWS.")
  nm_vals <- as.character(meta[nm_row, ])
  # Trim whitespace and strip any leading dashes (e.g. "-Fibroblast")
  nm_vals <- trimws(sub("^-+", "", nm_vals))

  # ── Build unified cell type label vector ──────────────────────────────────
  cell_labels <- dplyr::case_when(
    malignant_vals == "1"          ~ "Malignant",
    nm_vals != "" & nm_vals != "0" ~ nm_vals,
    TRUE                           ~ NA_character_
  )

  # ── Counts matrix ──────────────────────────────────────────────────────────
  # GSE103322 stores log2(TPM+1) values, not raw counts. Back-transform to
  # TPM-scale pseudo-counts (round(2^x - 1)) so spacexr::Reference gets integers.
  counts_raw <- as.matrix(raw[(N_META_ROWS + 1):nrow(raw), ])
  storage.mode(counts_raw) <- "numeric"
  counts_raw <- round(2^counts_raw - 1)
  storage.mode(counts_raw) <- "integer"

  # Keep only cells with a resolved cell type
  keep <- !is.na(cell_labels)
  counts_raw  <- counts_raw[, keep, drop = FALSE]
  cell_labels <- cell_labels[keep]
  names(cell_labels) <- colnames(counts_raw)

  # Drop cell types with fewer than 25 cells (spacexr minimum)
  type_counts <- table(cell_labels)
  valid_types <- names(type_counts[type_counts >= 25])
  keep2 <- cell_labels %in% valid_types
  if (any(!keep2)) {
    dropped <- names(type_counts[type_counts < 25])
    cat("Dropping low-count cell types (<25 cells):", paste(dropped, collapse = ", "), "\n")
    counts_raw  <- counts_raw[, keep2, drop = FALSE]
    cell_labels <- cell_labels[keep2]
  }

  cat("Cell type distribution in reference:\n")
  print(sort(table(cell_labels), decreasing = TRUE))

  counts_sp <- Matrix::Matrix(counts_raw, sparse = TRUE)
  reference <- spacexr::Reference(
    counts     = counts_sp,
    cell_types = factor(cell_labels),
    nUMI       = colSums(counts_sp)
  )

  saveRDS(reference, ref_rds)
  cat("\nPuram 2017 reference saved to:", ref_rds, "\n")
  return(reference)
}


#' Run RCTD deconvolution on one Visium Seurat object.
#'
#' Adds per-cell-type proportion columns and a dominant-cell-type column to
#' the Seurat object's metadata. Saves a spatial plot of the dominant cell type.
#'
#' @param seurat_obj   Seurat object (post-QC Visium sample).
#' @param reference    spacexr Reference object from get_hnscc_reference().
#' @param out_path     Directory for PNG output.
#' @param sample_name  Sample identifier used in output file names.
#' @param n_cores      Parallel cores for RCTD (default 4).
#' @return             Seurat object with RCTD proportion metadata added.
run_rctd <- function(seurat_obj, reference, out_path, sample_name, n_cores = 4) {
  cat("\nRunning RCTD deconvolution for:", sample_name, "\n")

  # ── Build SpatialRNA object ────────────────────────────────────────────────
  counts <- GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")

  coords <- GetTissueCoordinates(seurat_obj)
  # Keep only the x/y position columns (Seurat v5 returns "x", "y", "cell")
  coords <- as.data.frame(coords[, c("x", "y")])

  # Align barcodes between counts and coords
  common <- intersect(colnames(counts), rownames(coords))
  counts <- counts[, common, drop = FALSE]
  coords <- coords[common, , drop = FALSE]

  spatialRNA <- spacexr::SpatialRNA(
    coords = coords,
    counts = counts,
    nUMI   = colSums(counts)
  )

  # ── Run RCTD ──────────────────────────────────────────────────────────────
  # doublet_mode = "full" estimates proportions for all cell types per spot
  rctd <- spacexr::create.RCTD(spatialRNA, reference, max_cores = n_cores)
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = "full")

  # ── Extract proportions ───────────────────────────────────────────────────
  weights <- rctd@results$weights
  props   <- sweep(weights, 1, rowSums(weights), "/")   # normalize to sum to 1
  props   <- as.data.frame(as.matrix(props))
  # Prefix columns to avoid conflicts with existing metadata
  colnames(props) <- paste0("RCTD_", make.names(colnames(props)))

  # Add dominant (highest-proportion) cell type per spot
  type_cols <- colnames(props)
  props$RCTD_dominant <- sub("^RCTD_", "", type_cols[apply(props[, type_cols], 1, which.max)])

  seurat_obj <- AddMetaData(seurat_obj, metadata = props)

  cat("Cell type proportions added to metadata. Dominant types:\n")
  print(sort(table(seurat_obj$RCTD_dominant), decreasing = TRUE))

  # ── Spatial plot: dominant cell type ──────────────────────────────────────
  p <- SpatialDimPlot(seurat_obj, group.by = "RCTD_dominant", pt.size.factor = 2.5) +
    ggplot2::ggtitle(paste0(sample_name, " — RCTD dominant cell type"))

  png(file.path(out_path, paste0(sample_name, "_RCTD_dominant_cell_type.png")),
      width = 900, height = 800, res = 150)
  print(p)
  dev.off()

  # ── Spatial plots: per-cell-type proportion heatmaps ──────────────────────
  ct_features <- colnames(props)[grepl("^RCTD_", colnames(props)) & colnames(props) != "RCTD_dominant"]
  n_ct <- length(ct_features)
  cols <- min(3, n_ct)
  rows <- ceiling(n_ct / cols)

  prop_plots <- SpatialFeaturePlot(seurat_obj, features = ct_features,
                                   pt.size.factor = 2.5, ncol = cols)
  png(file.path(out_path, paste0(sample_name, "_RCTD_proportions.png")),
      width = cols * 700, height = rows * 600, res = 150)
  print(prop_plots)
  dev.off()

  cat("RCTD plots saved for:", sample_name, "\n")
  return(seurat_obj)
}
