# ============================================================
# Reorganize GEO-downloaded sample folders into Space Ranger layout
#
# For each sample folder:
#   1. Strips the "GSMxxxxxxx_PXX" prefix from every filename
#   2. Moves spatial-related files into a "spatial/" subfolder
#   3. Leaves matrix/barcodes/features files at the sample root
#
# Example:
#   GSM8655158_P04tissue_positions_list.csv.gz
#     -> spatial/tissue_positions_list.csv.gz
#   GSM8655158_P04matrix.mtx.gz
#     -> matrix.mtx.gz   (stays at root)
# ============================================================

library(stringr)

# ---- CONFIG ----
data_folder <- "CRC_project/GSE283052"   # parent folder containing sample subfolders

# Files that belong in the "spatial/" subfolder (Space Ranger convention)
spatial_patterns <- c(
  "tissue_positions",
  "scalefactors",
  "tissue_hires_image",
  "tissue_lowres_image",
  "aligned_fiducials",
  "detected_tissue_image"
)

# ---- FUNCTIONS ----

# Strip a GSM-style prefix like "GSM8655158_P04" from a filename
strip_gsm_prefix <- function(filename) {
  str_replace(filename, "^GSM\\d+_P\\d+", "")
}

# Decide if a file belongs in spatial/ based on its (already-stripped) name
is_spatial_file <- function(filename) {
  any(sapply(spatial_patterns, function(p) grepl(p, filename)))
}

reorganize_sample <- function(sample_path) {
  sample_name <- basename(sample_path)
  cat("\n::::::::::::::: Processing sample:", sample_name, ":::::::::::::::\n")
  
  files <- list.files(sample_path, full.names = FALSE)
  # Skip if already reorganized (avoid double-processing)
  files <- files[files != "spatial"]
  
  if (length(files) == 0) {
    cat("  No files found, skipping.\n")
    return(invisible(NULL))
  }
  
  spatial_dir <- file.path(sample_path, "spatial")
  dir.create(spatial_dir, showWarnings = FALSE)
  
  for (f in files) {
    old_path <- file.path(sample_path, f)
    
    # Skip directories (e.g. if spatial/ already has content)
    if (dir.exists(old_path)) next
    
    new_name <- strip_gsm_prefix(f)
    
    if (is_spatial_file(new_name)) {
      new_path <- file.path(spatial_dir, new_name)
    } else {
      new_path <- file.path(sample_path, new_name)
    }
    
    if (old_path == new_path) next  # nothing to do
    
    file.rename(old_path, new_path)
    cat("  ", f, " -> ", sub(paste0(sample_path, "/"), "", new_path), "\n")
  }
}

# ---- MAIN ----

folders <- list.files(data_folder, full.names = TRUE)
folders <- folders[dir.exists(folders)]  # only actual sample directories

cat(length(folders), "sample folder(s) found.\n")

for (folder in folders) {
  reorganize_sample(folder)
}

cat("\nAll samples reorganized.\n")