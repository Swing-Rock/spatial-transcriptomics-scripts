source("visium_scripts/st_helper.R")
library(Seurat)
library(patchwork)

# INPUT:
#   (str) data_path  — directory containing the Xenium output bundle
#                      (cell_feature_matrix.h5, cells.csv.gz, experiment.xenium, etc.)
#   (str) data_name  — sample identifier used for output file names
#   (bool) graph_output — whether to save QC plots (default TRUE)
#   (str) out_path   — directory for PNG output (default = data_path)
# OUTPUT: normalized Seurat object (SCTransform) after QC filtering

run_xenium_qc <- function(data_path, data_name, graph_output = TRUE, out_path = data_path, verbose_output = FALSE) {

  cat("=============== Starting QC for Xenium sample:", data_name, "===============\n")

  # ── Load ──────────────────────────────────────────────────────────────────
  cat("Loading Xenium data...\n")
  xen_obj <- LoadXenium(data_path, fov = "fov")

  # ── QC metrics ────────────────────────────────────────────────────────────
  # Mitochondrial content — Xenium panels often lack MT genes, which causes
  # PercentageFeatureSet to return NA for all cells and silently drops every
  # cell in the subset step. Fall back to 0 when no MT features are present.
  mt_features <- grep("^MT-", rownames(xen_obj), value = TRUE)
  if (length(mt_features) > 0) {
    xen_obj[["percent.mt"]] <- PercentageFeatureSet(xen_obj, features = mt_features)
  } else {
    xen_obj[["percent.mt"]] <- 0
    cat("Note: no mitochondrial features (MT-) in panel; percent.mt set to 0.\n")
  }

  # Ribosomal content — same guard
  ribo_features <- grep("^RP[SL]", rownames(xen_obj), value = TRUE)
  if (length(ribo_features) > 0) {
    xen_obj[["percent.ribo"]] <- PercentageFeatureSet(xen_obj, features = ribo_features)
  } else {
    xen_obj[["percent.ribo"]] <- 0
    cat("Note: no ribosomal features (RPS/RPL) in panel; percent.ribo set to 0.\n")
  }

  # Blank / negative-control codeword rate.
  # Xenium stores these as features prefixed "BLANK_" or "NegControlProbe_".
  blank_features <- grep("^(BLANK_|NegControlProbe_)", rownames(xen_obj), value = TRUE)
  if (length(blank_features) > 0) {
    xen_obj[["blank_rate"]] <- PercentageFeatureSet(xen_obj, features = blank_features)
  } else {
    xen_obj[["blank_rate"]] <- 0
    cat("Note: no blank/NegControl features detected in panel.\n")
  }

  cat("Raw cell count:", ncol(xen_obj), "\n")
  cat("QC summary (raw):\n")
  print(summary(xen_obj@meta.data[, c("nFeature_Xenium", "nCount_Xenium", "percent.mt", "blank_rate")]))

  # ── Filtering ─────────────────────────────────────────────────────────────
  # Thresholds tuned for targeted Xenium panels (typically 100-500 gene panels).
  # Adjust nFeature / nCount limits if using a larger panel.
  xen_sub <- subset(
    xen_obj,
    subset = nFeature_Xenium > 5       &   # at least 5 unique transcripts
             nCount_Xenium  > 10       &   # at least 10 total transcripts
             percent.mt     < 20       &   # <20% mitochondrial
             blank_rate     < 5            # <5% blank codewords
  )
  cat("Cells after QC:", ncol(xen_sub), "(removed", ncol(xen_obj) - ncol(xen_sub), ")\n")

  # ── Normalization ─────────────────────────────────────────────────────────
  cat("Normalizing with SCTransform...\n")
  # SCTransform residuals are stored in the "SCT" assay
  xen_norm <- SCTransform(xen_sub, assay = "Xenium", verbose = verbose_output)

  # ── Plots ─────────────────────────────────────────────────────────────────
  if (graph_output) {
    cat("Generating QC plots...\n")
    features  <- c("nFeature_Xenium", "nCount_Xenium", "percent.mt", "blank_rate")
    row_len   <- 4

    p_raw  <- make_v_plot(xen_obj,  features, row_len)
    p_filt <- make_v_plot(xen_sub,  features, row_len)
    p_norm <- make_v_plot(xen_norm, features, row_len)

    t_raw  <- title_plot("Raw Data")
    t_filt <- title_plot("After QC Filtering")
    t_norm <- title_plot("After SCTransform")

    combined <- t_raw / p_raw / t_filt / p_filt / t_norm / p_norm

    png(
      file.path(out_path, paste0(data_name, "_xenium_qc_plots.png")),
      width = 2200, height = 2800, res = 150
    )
    print(combined)
    dev.off()
    cat("QC plots saved.\n")
  }

  cat("=============== QC done for sample:", data_name, "===============\n")
  return(xen_norm)
}
