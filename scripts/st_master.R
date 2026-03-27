source("scripts/st_data_qc.R")
source("scripts/st_clustering.R")

data_folder = "C:/ongkeko lab works/spatial transcriptomics/data set 1/GSM8633893"
picture_folder = "C:/ongkeko lab works/spatial transcriptomics/data set 1/output"
start_time <- Sys.time()
batch_processing <- FALSE

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
  # Get all subdirectories (full paths)
  folders <- list.dirs(data_folder, full.names = TRUE, recursive = FALSE)
  n <- 0
  
  # Loop through folders
  for (folder in folders) {
    
    # Get folder name only
    folder_name <- basename(folder)
    
    # Skip the "output" folder
    if (folder_name == basename(picture_folder)) {
      next
    }
    cat("\n\n::::::::::::::: Processing sample:", folder_name, " :::::::::::::::\n")
    n = n + 1
    data <- run_qc(folder, folder_name, TRUE, picture_folder)
    clustering(data_obj = data, data_name = folder_name, out_path = picture_folder, single_cell = FALSE)
    
  }
  cat(n, " SAMPLE PROCESSED\n")
} else {
  data_name = basename(data_folder)
  data3 <- run_qc(data_folder, data_name, TRUE, picture_folder)
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