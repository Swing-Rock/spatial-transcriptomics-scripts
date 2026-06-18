source("visium_scripts/st_data_qc.R")
source("visium_scripts/st_clustering.R")
source("visium_scripts/st_merge_data.R")
source("visium_scripts/st_hpv_pathway_analysis.R")
source("visium_scripts/st_deconvolution.R")

# define path for everything 
# data folder
data_folder = "C:/ongkeko lab works/spatial transcriptomics/data set 1" 
# graphical output folder
picture_folder = "C:/ongkeko lab works/spatial transcriptomics/data set 1/output"
# data checkpoint folder
cache_folder  = file.path(picture_folder, "cache") 
dir.create(cache_folder, showWarnings = FALSE, recursive = TRUE)

#define operation parameters
#if there's multiple samples in the folder that need to be merged
batch_processing <- TRUE 
# RCTD deconvolution will be ran if true, if false, SingleR will be used to annotate
deconvolution <- TRUE 
#this only matters if we are not doing deconvolution. if we are, we can ignore this. 
#TRUE makes it so each spot is labeled a cell type before clustering. it can produce more cell types but it might be noisy
#FALSE makes it cluster before giving each cluster a cell type label. there might be less cell types but each type would be clearly defined
#i would run both and compare output if needed
single_cell_annotation <- TRUE 


start_time <- Sys.time()
HPV_plus_samples <- c("GSM8633893", "GSM8633894")

#----------all the functions that actually runs the code are here----------
batch_qc <- function() {
  qc_cache <- file.path(cache_folder, "data_list_qc.rds")
  if (file.exists(qc_cache)) {
    cat("\n\n::::::::::LOADING QC FROM CACHE::::::::::\n")
    data_list <- readRDS(qc_cache)
    return(data_list)
    
  } else {
    # Get all files in the folder (not recursive)
    folders <- list.files(data_folder, full.names = TRUE)
    #initiate variables
    n <- 0
    file_list <- list()
    
    # Loop through folders
    for (folder in folders) {
      # Get folder name only
      folder_name <- basename(folder)
      # Skip the "output" folder
      if (folder_name == basename(picture_folder)) {
        next
      }
      file_list[[folder_name]] <- folder   # key = file name, value = absolute path
      n = n + 1
    }
    cat(n, "SAMPLES WILL BE PROCESSED\n")
    data_list <- list()
    cat("no qc data found, running qc")
    for (name in names(file_list)) {
      cat("\n\n::::::::::::::: Processing sample:", name, " :::::::::::::::\n")
      data <- run_qc(file_list[[name]], name, TRUE, picture_folder)
      data_list[[name]] <- data
      data_list[[name]]$orig.ident <- name
    }
    rm(data)
    saveRDS(data_list, qc_cache)
    cat("\n\nQC cache saved to:", qc_cache, "\n")
    cat("\n\n::::::::::DONE WITH QC, STARTING DECONVOLUTION::::::::::\n")
    return(data_list)
    
  }
  
}

single_sample_qc <- function(data_name){
  qc_cache <- file.path(cache_folder, paste0(data_name, "_qc.rds"))
  if (file.exists(qc_cache)) {
    cat("\n\n::::::::::LOADING QC FROM CACHE::::::::::\n")
    return (data <- readRDS(qc_cache))
  } else {
    data <- run_qc(data_folder, data_name, TRUE, picture_folder)
    saveRDS(data, qc_cache)
    cat("\n\nQC cache saved to:", qc_cache, "\n")
    return (data)
  }
}


per_sample_deconv <- function(data_list){
  deconv_cache <- file.path(cache_folder, "data_list_deconv.rds")
  if (file.exists(deconv_cache)) {
    cat("\n\n::::::::::LOADING DECONVOLUTION FROM CACHE::::::::::\n")
    data_list <- readRDS(deconv_cache)
    return (data_list)
  } else {
    hnscc_ref <- get_hnscc_reference(cache_folder)
    for (name in names(data_list)) {
      cat("\n\n::::::::::::::: Deconvolving sample:", name, " :::::::::::::::\n")
      data_list[[name]] <- run_rctd(data_list[[name]], hnscc_ref, picture_folder, name)
    }
    saveRDS(data_list, deconv_cache)
    cat("\n\nDeconvolution cache saved to:", deconv_cache, "\n")
  }
  cat("\n\n::::::::::DONE WITH DECONVOLUTION, MERGING SAMPLES::::::::::\n")
  return(data_list)
}

merge_data <- function(data_list){
  merge_cache <- file.path(cache_folder, "merged_data.rds")
  if (file.exists(merge_cache)) {
    cat("\n\n::::::::::LOADING MERGED DATA FROM CACHE::::::::::\n")
    merged_data <- readRDS(merge_cache)
    return (merged_data)
  } else {
    cat("no merged data found, merging samples")
    merged_data <- merge_dataset(data_list, picture_folder)
    merged_data <- FindVariableFeatures(merged_data, assay = "SCT", selection.method = "vst", nfeatures = 3000)
    saveRDS(merged_data, merge_cache)
    cat("\n\nMerged data cache saved to:", merge_cache, "\n")
  }
  return (merged_data)
}

cluster_n_annotate <- function(merged_data){
  cluster_cache <- file.path(cache_folder, "clustered_data.rds")
  if (file.exists(cluster_cache)) {
    cat("\n\n::::::::::LOADING CLUSTERED DATA FROM CACHE::::::::::\n")
    merged_data <- readRDS(cluster_cache)
    return(merged_data)
  } else {
    cat("no annotated data found, annotating samples")
    merged_data <- clustering(data_obj = merged_data, data_name = "merged", out_path = picture_folder, single_cell = FALSE)
    saveRDS(merged_data, cluster_cache)
    cat("\n\nClustered data cache saved to:", cluster_cache, "\n")
    return (merged_data)
  }
}


single_sample_cluster <- function(data, data_name){
  cluster_cache <- file.path(cache_folder, paste0(data_name, "_clustered.rds"))
  if (file.exists(cluster_cache)) {
    cat("\n\n::::::::::LOADING CLUSTERED DATA FROM CACHE::::::::::\n")
    annotated_data <- readRDS(cluster_cache)
    return (annotated_data)
  } else {
    annotated_data <- clustering(data, data_name, picture_folder)
    saveRDS(merged_data, cluster_cache)
    cat("\n\nClustered data cache saved to:", cluster_cache, "\n")
    return (annotated_data)
  }
}


hpv_analysis <- function(merged_data, deconvolution){
  cat("\n\n::::::::::STARTING HPV PATHWAY ANALYSIS::::::::::\n")
  
  # 1. Tag each spot with HPV status based on sample ID
  merged_data <- add_hpv_status(merged_data, HPV_plus_samples)
  
  # 2. Isolate malignant spots
  #    !! Adjust MALIGNANT_LABELS in st_hpv_pathway_analysis.R to match
  #    !! whatever labels your annotation used (e.g. "Tumor", "Epithelial")
  malignant <- subset_malignant(merged_data, deconvolution, MALIGNANT_LABELS)
  
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
#--------------------end of functions--------------------


# ---------- FOR RUNNING SCRIPT ON A SINGLE SAMPLE: ----------
# DATA FOLDER FORMAT:
#   (Folder) Sample name
#   -> (File) (SampleName)_filtered_feature_bc_matrix.h5
#   -> (Folder) (SampleName)
#   ->-> (File) tissue_lowres_image.png
#   ->-> (File) scalefactors_json.json
#   ->-> (File) tissue_positions_list.csv
# data_folder IS THE FOLDER CONTAINING THE .h5
# ---------- FOR RUNNING SCRIPT ON A FOLDER CONTAINING MULTIPLE SAMPLES: ----------
# MAKE SURE FOLDER IS FORMATTED AS FOLLOW: 
#  (Folder) whatever_name
#   -> (Folder) (SampleName)
#   -> (Folder) (SampleName)
#   -> (Folder) (SampleName)
#   -> (Folder) (output folder): not required
# EACH CHILD FOLDER SHOULD HAVE THE SAME FORMAT AS DESCRIBED UNDER 'RUNNING SCRIPT ON A SINGLE SAMPLE' 
# data_folder IS THE whatever_name FOLDER
# ---------


if (batch_processing){
  data_list = batch_qc()

  if (deconvolution)  {
    data_list = per_sample_deconv(data_list)
    merged_data = merge_data(data_list)
  } else {
    merged_data = merge_data(data_list)
    merged_data = cluster_n_annotate(merged_data)
  }
  
} else { #NOTE TO SELF: REMEMBER TO ADD DECONVOLUTION CODE TO THIS SOMETIMES... HOPEFULLY this code actually works idk when's the last time ive used it tbh
  data_name = basename(data_folder)
  dataObj = single_sample_qc(data_name)
  merged_data = single_sample_cluster(dataObj, data_name)
}

hpv_analysis(merged_data, deconvolution)



end_time <- Sys.time()
cat("PROGRAM START TIME:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("PROGRAM END TIME:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

total_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
# Convert to hours, minutes, seconds
hours <- floor(total_secs / 3600)
minutes <- floor((total_secs %% 3600) / 60)
seconds <- round(total_secs %% 60, 2)
# Print nicely
cat("TOTAL RUN TIME: ",
    hours, "h", 
    minutes, "m", 
    seconds, "s\n")

cat("\nALL DONE :D\n")

