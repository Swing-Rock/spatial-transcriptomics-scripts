source("scripts/xenium_qc.R")
source("scripts/xenium_annotation.R")
source("scripts/xenium_merge.R")

# ── Paths ─────────────────────────────────────────────────────────────────────
# data_folder  : root folder containing one sub-folder per sample,
#                OR the single sample directory if batch_processing = FALSE
# picture_folder: where all PNG outputs are written
# cache_folder : where .rds checkpoints are saved

data_folder    <- "C:/ongkeko lab works/spatial transcriptomics/xenium dataset/GSM9054471"
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

batch_processing       <- FALSE
single_cell_annotation <- FALSE


start_time <- Sys.time()

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

# -- HPV analysis --------------------------------------------------------------
hpv_analysis <- function(merged_data, deconvolution){
  cat("\n\n::::::::::STARTING HPV PATHWAY ANALYSIS::::::::::\n")
  
  # 1. Tag each spot with HPV status based on sample ID
  merged_data <- add_hpv_status(merged_data, HPV_plus_samples)
  
  # 2. Isolate malignant spots
  #    !! Adjust MALIGNANT_LABELS in st_hpv_pathway_analysis.R to match
  #    !! whatever labels your annotation used (e.g. "Tumor", "Epithelial")
  malignant <- subset_malignant(merged_data, FALSE, MALIGNANT_LABELS)
  
  # 3. DE: HPV+ vs HPV- within malignant compartment
  de_cache <- file.path(cache_folder, "de_results.rds")
  if (file.exists(de_cache)) {
    cat("\n\n::::::::::LOADING DE RESULTS FROM CACHE::::::::::\n")
    de_results <- readRDS(de_cache)
  } else {
    cat("\nRunning differential expression...\n")
    DefaultAssay(merged_data) <- "Spatial"
    malignant <- NormalizeData(merged_data, assay = "Spatial")
    malignant <- JoinLayers(malignant, assay = "Spatial")
    de_results <- run_de(malignant, picture_folder)
    de_results$gene <- rownames(de_results)
    saveRDS(de_results, de_cache)
    cat("\n\nDE cache saved to:", de_cache, "\n")
  }
  
  # 4. GO biological process enrichment
  cat("\nRunning GO enrichment...\n")
  go_results <- run_go(de_results, picture_folder)
  
  # 5. KEGG pathway enrichment
  cat("\nRunning KEGG enrichment...\n")
  kegg_results <- run_kegg(de_results, picture_folder)
  
  # 6. Score top GO pathways spatially across all spots
  cat("\nScoring pathways spatially...\n")
  merged_data <- run_spatial_scoring(merged_data, go_results, picture_folder)
  
  cat("\n\n::::::::::HPV PATHWAY ANALYSIS COMPLETE::::::::::\n")
}

# ── Run ───────────────────────────────────────────────────────────────────────
if (batch_processing) {
  data_list   <- batch_qc()
  data_list   <- per_sample_annotation(data_list)
  merged_data <- merge_samples(data_list)
} else {
  merged_data <- single_sample_workflow()
}


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
