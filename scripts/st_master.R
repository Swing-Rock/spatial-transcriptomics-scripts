source("scripts/st_data_qc.R")
source("scripts/st_clustering.R")
source("scripts/st_merge_data.R")


data_folder = "C:/ongkeko lab works/spatial transcriptomics/data set 1"
picture_folder = "C:/ongkeko lab works/spatial transcriptomics/data set 1/output"
start_time <- Sys.time()
batch_processing <- TRUE


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



if (batch_processing){
  
  # Get all files in the folder (not recursive)
  folders <- list.files(data_folder, full.names = TRUE)
  #initiate variables
  n <- 0
  file_list <- list()
  data_list <- list()
  
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
  
  
  # run QC on each sample
  for (name in names(file_list)) {
    cat("\n\n::::::::::::::: Processing sample:", name, " :::::::::::::::\n")
    data <- run_qc(file_list[[name]], name, TRUE, picture_folder)
    data_list[[name]] <- data
    data_list[[name]]$orig.ident <- name
  }
  rm(data)
  cat("\n\n::::::::::DONE WITH QC, MERGING SAMPLES::::::::::\n")
  
  merged_data <- merge_dataset(data_list, picture_folder)
  merged_data <- FindVariableFeatures(merged_data, assay = "SCT", selection.method = "vst", nfeatures = 3000)
  clustering(data_obj = merged_data, data_name = "merged", out_path = picture_folder, single_cell = TRUE)
    
} else {
  data_name = basename(data_folder)
  data <- run_qc(data_folder, data_name, TRUE, picture_folder)
  clustering(data, data_name, picture_folder)
}


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