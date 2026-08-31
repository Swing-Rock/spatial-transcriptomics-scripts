source("scripts/xenium_qc.R")
source("scripts/xenium_annotation.R")
source("scripts/xenium_merge.R")
source("visium_scripts/st_hpv_pathway_analysis.R")
options(future.globals.maxSize = 2000 * 1024^2)  # 2GB

# ── Paths ─────────────────────────────────────────────────────────────────────
# data_folder  : root folder containing one sub-folder per sample,
#                OR the single sample directory if batch_processing = FALSE
# picture_folder: where all PNG outputs are written
# cache_folder : where .rds checkpoints are saved

# Point to the PARENT folder that contains all Xenium sample sub-folders
data_folder    <- "GSE300147"
picture_folder <- file.path(data_folder, "output")
cache_folder   <- file.path(picture_folder, "cache")
dir.create(cache_folder, showWarnings = FALSE, recursive = TRUE)

# ── Parameters ────────────────────────────────────────────────────────────────
# batch_processing  : TRUE  → every sub-folder in data_folder is a sample
#                     FALSE → data_folder IS the sample directory
# single_cell_annotation:
#   TRUE  → label each cell individually with SingleR (more labels, more noise)
#   FALSE → cluster first, then label each cluster (fewer but cleaner cell types)
#   Tip: run both and compare results

batch_processing       <- TRUE
single_cell_annotation <- FALSE

# ── HPV Sample Labels ─────────────────────────────────────────────────────────
# Known HPV status from GEO metadata (p16 IHC).
# Samples listed as "N/A" in GEO will be classified by infer_hpv_status()
# after merging, using proxy gene expression.
# HPV_plus_samples will be set after infer_hpv_status() runs (see Run section)
KNOWN_HPV_PLUS  <- c("GSM9054471", "GSM9054472", "GSM9054473", 
                     "GSM9054474", "GSM9054475", "GSM9054476", 
                     "GSM9054477","GSM9054486", "GSM9054487")
KNOWN_HPV_MINUS <- c("GSM9054483", "GSM9054484")
HPV_STATUS_UNKNOWN <- c("GSM9054478", "GSM9054479", "GSM9054480", 
                        "GSM9054481", "GSM9054482", "GSM9054485", "GSM9054488")

start_time <- Sys.time()
MALIGNANT_LABELS <- c("Epithelial_cells", "Keratinocytes")

# ── DATA FOLDER FORMAT ────────────────────────────────────────────────────────
# BATCH MODE (batch_processing = TRUE):
#   (Folder) data_folder/
#     -> (Folder) SampleName1/
#          cell_feature_matrix.h5
#          cells.csv.gz  (or cells.parquet)
#          experiment.xenium
#          morphology.ome.tif  (optional)
#     -> (Folder) SampleName2/
#     -> (Folder) output/   <- excluded automatically
#
# SINGLE SAMPLE (batch_processing = FALSE):
#   data_folder/
#     cell_feature_matrix.h5
#     cells.csv.gz
#     experiment.xenium
# ─────────────────────────────────────────────────────────────────────────────


# ── Helper: QC one sample (with caching) ─────────────────────────────────────
run_sample_qc <- function(sample_path, sample_name) {
  qc_cache <- file.path(cache_folder, paste0(sample_name, "_qc.rds"))
  if (file.exists(qc_cache)) {
    cat("\n::::::::::LOADING QC FROM CACHE:", sample_name, "::::::::::\n")
    return(readRDS(qc_cache))
  }
  data <- run_xenium_qc(sample_path, sample_name, graph_output = TRUE, out_path = picture_folder)
  saveRDS(data, qc_cache)
  cat("\nQC cache saved:", qc_cache, "\n")
  return(data)
}


# ── Batch QC ─────────────────────────────────────────────────────────────────
batch_qc <- function() {
  batch_cache <- file.path(cache_folder, "data_list_qc.rds")
  if (file.exists(batch_cache)) {
    cat("\n\n::::::::::LOADING BATCH QC FROM CACHE::::::::::\n")
    return(readRDS(batch_cache))
  }

  all_entries <- list.files(data_folder, full.names = TRUE)
  # Keep only directories, skip the output folder
  sample_dirs <- all_entries[
    file.info(all_entries)$isdir &
    basename(all_entries) != basename(picture_folder)
  ]

  cat(length(sample_dirs), "SAMPLES WILL BE PROCESSED\n")

  data_list <- list()
  for (sample_path in sample_dirs) {
    sample_name <- basename(sample_path)
    cat("\n\n::::::::::::::: Processing sample:", sample_name, ":::::::::::::::\n")
    data <- run_xenium_qc(sample_path, sample_name, graph_output = TRUE, out_path = picture_folder)
    data$orig.ident <- sample_name
    data_list[[sample_name]] <- data
  }
  rm(data)

  saveRDS(data_list, batch_cache)
  cat("\n\nBatch QC cache saved:", batch_cache, "\n")
  cat("\n\n::::::::::DONE WITH QC, STARTING ANNOTATION::::::::::\n")
  return(data_list)
}


# ── Per-sample annotation (with caching) ─────────────────────────────────────
per_sample_annotation <- function(data_list) {
  annot_cache <- file.path(cache_folder, "data_list_annotated.rds")
  if (file.exists(annot_cache)) {
    cat("\n\n::::::::::LOADING ANNOTATION FROM CACHE::::::::::\n")
    return(readRDS(annot_cache))
  }

  for (name in names(data_list)) {
    cat("\n\n::::::::::::::: Annotating sample:", name, ":::::::::::::::\n")
    data_list[[name]] <- run_xenium_annotation(
      data_obj    = data_list[[name]],
      data_name   = name,
      out_path    = picture_folder,
      single_cell = single_cell_annotation
    )
  }

  saveRDS(data_list, annot_cache)
  cat("\n\nAnnotation cache saved:", annot_cache, "\n")
  cat("\n\n::::::::::DONE WITH ANNOTATION, MERGING SAMPLES::::::::::\n")
  return(data_list)
}


# ── Merge (with caching) ──────────────────────────────────────────────────────
merge_samples <- function(data_list) {
  merge_cache <- file.path(cache_folder, "merged_data.rds")
  if (file.exists(merge_cache)) {
    cat("\n\n::::::::::LOADING MERGED DATA FROM CACHE::::::::::\n")
    return(readRDS(merge_cache))
  }

  merged_data <- merge_xenium(data_list, picture_folder)
  saveRDS(merged_data, merge_cache)
  cat("\n\nMerged data cache saved:", merge_cache, "\n")
  return(merged_data)
}


# ── Single-sample workflow ────────────────────────────────────────────────────
single_sample_workflow <- function() {
  sample_name <- basename(data_folder)

  data <- run_sample_qc(data_folder, sample_name)
  data$orig.ident <- sample_name

  data <- run_xenium_annotation(
    data_obj    = data,
    data_name   = sample_name,
    out_path    = picture_folder,
    single_cell = single_cell_annotation
  )
  return(data)
}

# ── Infer HPV status for N/A samples ─────────────────────────────────────────
# Uses proxy gene expression to classify samples with unknown p16 status.
# Strategy:
#   1. Score each sample using genes known to be elevated in HPV+ HNSCC
#      (CXCL10, CD274, MKI67, PCNA — immune response and E7-driven proliferation)
#   2. Compute HPV+ and HPV- centroids from the known-status samples
#   3. Classify unknown samples by proximity to each centroid
#   4. Validate: HPV+ centroid must score higher than HPV- centroid
infer_hpv_status <- function(seurat_obj, out_path) {
  cat("\n\n::::::::::INFERRING HPV STATUS FOR UNKNOWN SAMPLES::::::::::\n")

  # Proxy genes: elevated in HPV+ HNSCC
  #   CXCL10, CD274 — interferon/immune checkpoint response, consistently higher in HPV+
  #   MKI67, PCNA   — proliferation driven by HPV E7 (inactivates pRb → cell cycle entry)
  proxy_genes <- c("CXCL10", "CD274", "MKI67", "PCNA")
  available   <- proxy_genes[proxy_genes %in% rownames(seurat_obj)]

  if (length(available) < 2) {
    stop("Fewer than 2 proxy genes found in panel. Cannot infer HPV status.")
  }
  cat("Proxy genes available in panel:", paste(available, collapse = ", "), "\n")

  # Per-cell module score
  seurat_obj <- AddModuleScore(seurat_obj,
                               features = list(available),
                               name     = "HPV_proxy",
                               assay    = "SCT",
                               ctrl     = 5)
  score_col <- "HPV_proxy1"   # Seurat appends "1"

  # Aggregate to per-sample mean score
  library(dplyr)
  known_lookup <- c(
    setNames(rep("Positive", length(KNOWN_HPV_PLUS)),  KNOWN_HPV_PLUS),
    setNames(rep("Negative", length(KNOWN_HPV_MINUS)), KNOWN_HPV_MINUS)
  )

  sample_scores <- seurat_obj@meta.data |>
    as.data.frame() |>
    dplyr::select(orig.ident, score = dplyr::all_of(score_col)) |>
    dplyr::group_by(orig.ident) |>
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(known_p16 = known_lookup[orig.ident],
                  known_p16 = ifelse(is.na(known_p16), "Unknown", known_p16))

  # Centroids from known samples
  centroid_pos <- mean(sample_scores$mean_score[sample_scores$known_p16 == "Positive"])
  centroid_neg <- mean(sample_scores$mean_score[sample_scores$known_p16 == "Negative"])
  threshold    <- (centroid_pos + centroid_neg) / 2

  cat(sprintf("HPV+ centroid: %.4f | HPV- centroid: %.4f | threshold: %.4f\n",
              centroid_pos, centroid_neg, threshold))

  if (centroid_pos <= centroid_neg) {
    warning(paste(
      "HPV+ centroid is NOT higher than HPV- centroid.",
      "Proxy gene direction may not hold for this dataset.",
      "Review hpv_status_classification.csv and set HPV_plus_samples manually if needed."
    ))
  }

  # Classify all samples
  sample_scores <- sample_scores |>
    mutate(
      hpv_status = case_when(
        known_p16 == "Positive"              ~ "HPV+ (known)",
        known_p16 == "Negative"              ~ "HPV- (known)",
        mean_score >= threshold              ~ "HPV+ (inferred)",
        TRUE                                 ~ "HPV- (inferred)"
      )
    )

  cat("\nFull HPV classification:\n")
  print(as.data.frame(sample_scores[, c("orig.ident", "mean_score", "known_p16", "hpv_status")]))

  # Bar plot ordered by score
  status_colors <- c(
    "HPV+ (known)"    = "#E64B35",
    "HPV- (known)"    = "#4DBBD5",
    "HPV+ (inferred)" = "#FF9D9A",
    "HPV- (inferred)" = "#9DCFDE"
  )
  p <- ggplot(sample_scores,
              aes(x = reorder(orig.ident, mean_score), y = mean_score, fill = hpv_status)) +
    geom_col() +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "grey30", linewidth = 0.6) +
    annotate("text", x = 0.6, y = threshold + 0.002, label = "threshold", hjust = 0,
             size = 3, color = "grey30") +
    scale_fill_manual(values = status_colors) +
    labs(
      title    = "Per-sample HPV proxy signature score",
      subtitle = paste("Proxy genes:", paste(available, collapse = ", ")),
      x = "Sample", y = "Mean module score", fill = "HPV status"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  ggsave(file.path(out_path, "hpv_inferred_status.png"), p, width = 10, height = 5, dpi = 150)
  cat("Classification plot saved:", file.path(out_path, "hpv_inferred_status.png"), "\n")

  # Save table for manual review
  write.csv(sample_scores, file.path(out_path, "hpv_status_classification.csv"), row.names = FALSE)
  cat("Classification table saved:", file.path(out_path, "hpv_status_classification.csv"), "\n")

  # Return all HPV+ sample IDs (confirmed + inferred)
  hpv_plus <- sample_scores$orig.ident[grepl("HPV\\+", sample_scores$hpv_status)]
  cat("\nFinal HPV_plus_samples:\n")
  print(hpv_plus)
  return(hpv_plus)
}


# -- HPV analysis --------------------------------------------------------------
hpv_analysis <- function(merged_data) {
  hpv_cache <- file.path(cache_folder, "hpv_analysis.rds")
  if (file.exists(hpv_cache)) {
    cat("\n\n::::::::::LOADING HPV ANALYSIS FROM CACHE::::::::::\n")
    return(readRDS(hpv_cache))
  }

  cat("\n\n::::::::::STARTING HPV PATHWAY ANALYSIS::::::::::\n")

  # 1. Tag each cell with HPV status based on sample ID (orig.ident)
  #    HPV_plus_samples is defined in the config section above
  merged_data <- add_hpv_status(merged_data, HPV_plus_samples)

  # 2. Isolate malignant cells
  #    !! Adjust MALIGNANT_LABELS in visium_scripts/st_hpv_pathway_analysis.R
  #    !! to match your SingleR labels (e.g. "Epithelial_cells", "Keratinocytes")
  malignant <- subset_malignant(merged_data, deconv = FALSE, malignant_labels = MALIGNANT_LABELS)

  # 3. DE: HPV+ vs HPV- within malignant compartment
  de_cache <- file.path(cache_folder, "de_results.rds")
  if (file.exists(de_cache)) {
    cat("\n\n::::::::::LOADING DE RESULTS FROM CACHE::::::::::\n")
    de_results <- readRDS(de_cache)
  } else {
    cat("\nRunning differential expression...\n")
    #malignant <- JoinLayers(malignant, assay = "SCT")
    de_results <- run_de(malignant, picture_folder, platform = "xenium")
    saveRDS(de_results, de_cache)
    cat("\n\nDE cache saved to:", de_cache, "\n")
  }

  # 4. GO biological process enrichment
  cat("\nRunning GO enrichment...\n")
  go_results <- run_go(de_results, picture_folder)

  # 5. KEGG pathway enrichment
  cat("\nRunning KEGG enrichment...\n")
  kegg_results <- run_kegg(de_results, picture_folder)

  # 6. Score top GO pathways spatially across all cells
  cat("\nScoring pathways spatially...\n")
  merged_data <- run_spatial_scoring(merged_data, go_results, picture_folder, platform = "xenium")

  cat("\n\n::::::::::HPV PATHWAY ANALYSIS COMPLETE::::::::::\n")

  result <- list(
    merged_data  = merged_data,
    go_results   = go_results,
    kegg_results = kegg_results
  )
  saveRDS(result, hpv_cache)
  cat("\n\nHPV analysis cache saved to:", hpv_cache, "\n")
  return(result)
}

# ── Run ───────────────────────────────────────────────────────────────────────
if (batch_processing) {
  data_list   <- batch_qc()
  data_list   <- per_sample_annotation(data_list)
  merged_data <- merge_samples(data_list)
  # Release per-sample list from global env after merge
  rm(data_list); gc()
} else {
  merged_data <- single_sample_workflow()
}

# Infer HPV status for the N/A samples, then combine with known labels
HPV_plus_samples <- infer_hpv_status(merged_data, picture_folder)

hpv_results <- hpv_analysis(merged_data)


# ── Timer ─────────────────────────────────────────────────────────────────────
end_time   <- Sys.time()
total_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
hours      <- floor(total_secs / 3600)
minutes    <- floor((total_secs %% 3600) / 60)
seconds    <- round(total_secs %% 60, 2)

cat("\nPROGRAM START TIME:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("PROGRAM END TIME:  ", format(end_time,   "%Y-%m-%d %H:%M:%S"), "\n")
cat("TOTAL RUN TIME:    ", hours, "h", minutes, "m", seconds, "s\n")
cat("\nALL DONE :D\n")
