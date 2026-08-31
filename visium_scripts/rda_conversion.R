# ============================================================
# export_rda_to_spaceranger.R
#
# Converts .rda spatial Seurat objects (e.g. ST-colon1.rda) into a
# Space-Ranger-style folder structure:
#
#   <output_root>/<sample_name>/filtered_feature_bc_matrix.h5
#   <output_root>/<sample_name>/spatial/tissue_positions.csv
#   <output_root>/<sample_name>/spatial/scalefactors_json.json
#   <output_root>/<sample_name>/spatial/tissue_hires_image.png
#   <output_root>/<sample_name>/spatial/tissue_lowres_image.png
#
# This is the FIRST format your cell2location Python script already
# checks for (sq.read.visium()), so no Python changes are needed -
# just point INPUT_FOLDER at a directory containing both your real
# Space Ranger folders and these exported ones side by side.
#
# Reuses load_spatial_object() from st_data_qc.R, so all the version-
# compatibility patching we already debugged (VisiumV1 "misc" slot,
# UpdateSeuratObject, image/assay alignment) applies here too.
# ============================================================

library(Seurat)
library(Matrix)
library(jsonlite)
library(png)

# DropletUtils gives us a proper, well-tested 10x-format HDF5 writer.
# It's a Bioconductor package - install once with:
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("DropletUtils")
if (!requireNamespace("DropletUtils", quietly = TRUE)) {
  stop("Package 'DropletUtils' is required to write the 10x-format .h5 file.\n",
       "Install it with:\n",
       '  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")\n',
       '  BiocManager::install("DropletUtils")')
}
load_spatial_object <- function(data_path, data_name) {
  
  if (tolower(tools::file_ext(data_path)) == "rda") {
    if (!file.exists(data_path)) {
      stop("File not found: ", data_path)
    }
    
    # .rda restores variables under whatever name they were saved as,
    # not the filename - load into a scratch env and search for the object.
    env <- new.env()
    loaded_names <- load(data_path, envir = env)
    cat("Loaded", data_path, "-> variable(s):", paste(loaded_names, collapse = ", "), "\n")
    
    obj <- NULL
    if (length(loaded_names) == 1) {
      obj <- get(loaded_names[1], envir = env)
    } else {
      for (nm in loaded_names) {
        candidate <- get(nm, envir = env)
        if (inherits(candidate, "Seurat")) {
          obj <- candidate
          break
        }
      }
    }
    
    # unwrap if it's a list containing a Seurat object rather than the object itself
    if (!is.null(obj) && !inherits(obj, "Seurat") && is.list(obj)) {
      inner <- Filter(function(x) inherits(x, "Seurat"), obj)
      if (length(inner) >= 1) obj <- inner[[1]]
    }
    
    if (is.null(obj) || !inherits(obj, "Seurat")) {
      stop("Could not find a Seurat object inside ", data_path,
           ". Variable(s) found: ", paste(loaded_names, collapse = ", "),
           ". Class of first object: ",
           if (length(loaded_names) >= 1) class(get(loaded_names[1], envir = env))[1] else "none")
    }
    
    # --- patch old-format image slots BEFORE calling any Seurat accessor ---
    # .rda files saved with an older SeuratObject version may carry image
    # objects (e.g. class "VisiumV1") that are missing slots the currently
    # installed package's class definition now expects (commonly "misc").
    # Calling Assays()/Images()/etc. on the raw object triggers R's S4
    # validObject() check and errors with:
    #   "invalid class VisiumV1 object: slots in class definition but not in object: misc"
    # Fix: rebuild each image under the CURRENT class definition, copying
    # over whatever slots the old object actually has and letting any new
    # slots (like "misc") take their default value. This uses raw slot()
    # access (not the Seurat generics), which does not trigger validation.
    if (length(obj@images) > 0) {
      for (img_name in names(obj@images)) {
        old_img <- obj@images[[img_name]]
        img_class <- class(old_img)[1]
        new_img <- tryCatch(methods::new(img_class), error = function(e) NULL)
        if (is.null(new_img)) {
          warning("Could not rebuild image slot '", img_name, "' (class ", img_class,
                  ") for sample '", data_name, "' - leaving as-is, downstream spatial plots may fail.")
          next
        }
        for (s in methods::slotNames(img_class)) {
          val <- tryCatch(methods::slot(old_img, s), error = function(e) NULL)
          if (!is.null(val)) {
            methods::slot(new_img, s) <- val
          }
        }
        obj@images[[img_name]] <- new_img
      }
      cat("  patched", length(obj@images), "image slot(s) for version compatibility\n")
    }
    
    # Seurat's own plotting/analysis functions (SpatialPlot, SCTransform, etc.)
    # check an internal version stamp on the object and refuse to run until
    # it's been brought up to date - separate from the image class patch above.
    obj <- suppressWarnings(UpdateSeuratObject(obj))
    cat("  ran UpdateSeuratObject() for version compatibility\n")
    
    # Belt-and-braces: UpdateSeuratObject doesn't always catch every legacy
    # naming mismatch. If the image's spot coordinates use a different cell
    # naming convention than the assay (e.g. old "_1" barcode suffix vs.
    # current "-1"), SpatialPlot fails with a generic "run UpdateSeuratObject"
    # error even after that's already been run. Detect and fix that directly:
    # if cell names match one-to-one once suffixes are stripped, realign and
    # rename the image's coordinates to the assay's actual cell names.
    assay_cells <- Cells(obj)
    strip_suffix <- function(x) sub("[-_][0-9]+$", "", x)
    if (length(obj@images) > 0) {
      for (img_name in names(obj@images)) {
        img <- obj@images[[img_name]]
        coord_cells <- rownames(img@coordinates)
        if (is.null(coord_cells)) next
        if (setequal(assay_cells, coord_cells)) next  # already aligned, nothing to do
        
        base_assay <- strip_suffix(assay_cells)
        base_coord <- strip_suffix(coord_cells)
        if (length(base_assay) == length(base_coord) &&
            setequal(base_assay, base_coord)) {
          match_idx <- match(base_assay, base_coord)
          img@coordinates <- img@coordinates[match_idx, , drop = FALSE]
          rownames(img@coordinates) <- assay_cells
          obj@images[[img_name]] <- img
          cat("  realigned image '", img_name,
              "' spot names to match assay cell names (suffix mismatch fixed)\n", sep = "")
        } else {
          warning("Image '", img_name, "' coordinates do not match assay cell names, ",
                  "even after stripping barcode suffixes, for sample '", data_name,
                  "'. Spatial plots will likely fail or misalign spots - this needs manual review.")
        }
      }
    }
    
    # normalize assay name to "Spatial" so downstream code (percent.mt calc,
    # SCTransform(assay="Spatial"), AggregateExpression()$Spatial etc.) works untouched
    assays_present <- Assays(obj)
    if (!"Spatial" %in% assays_present) {
      alt_names <- intersect(c("SCT", "RNA", "Spatial.008um", "Spatial.016um"), assays_present)
      if (length(alt_names) >= 1) {
        cat("  no 'Spatial' assay found; renaming '", alt_names[1], "' -> 'Spatial'\n", sep = "")
        obj[["Spatial"]] <- obj[[alt_names[1]]]
        DefaultAssay(obj) <- "Spatial"
      } else {
        warning("No 'Spatial' or recognizable alternate assay found in ", data_path,
                ". Assays present: ", paste(assays_present, collapse = ", "))
      }
    }
    
    if (length(Images(obj)) == 0) {
      warning("Sample '", data_name, "' has no @images slot populated - ",
              "spot coordinates/scale factors may be missing from this .rda.")
    }
    
    return(obj)
  }
  
  # ---- original folder-based loading path (unchanged) ----
  data.dir <- data_path
  h5_path  <- file.path(data.dir, "filtered_feature_bc_matrix.h5")
  h5ad_files <- list.files(data.dir, pattern = "\\.h5ad$")
  if (file.exists(h5_path)) {
    return(Load10X_Spatial(data.dir, filename = "filtered_feature_bc_matrix.h5"))
  } else if (length(h5ad_files) > 0) {
    cat("No filtered_feature_bc_matrix.h5 found; loading from h5ad file instead.\n")
    return(load_h5ad_visium(data.dir))
  } else {
    stop("No valid input file found in: ", data.dir,
         "\n  Expected: filtered_feature_bc_matrix.h5, a .h5ad file, or pass a .rda file path directly.")
  }
}

#' Export a single .rda sample to a Space-Ranger-style folder.
#'
#' @param rda_path path to the .rda file (e.g. "ST-colon1.rda")
#' @param output_root directory under which "<sample_name>/" will be created
#' @param sample_name optional override; defaults to the .rda filename minus extension
export_rda_to_spaceranger <- function(rda_path, output_root, sample_name = NULL) {
  
  if (is.null(sample_name)) {
    sample_name <- tools::file_path_sans_ext(basename(rda_path))
  }
  
  cat("=============== Exporting", sample_name, "===============\n")
  
  # --- load + patch (reuses everything we already fixed for the R QC pipeline) ---
  obj <- load_spatial_object(rda_path, sample_name)
  
  out_dir <- file.path(output_root, sample_name)
  spatial_dir <- file.path(out_dir, "spatial")
  dir.create(spatial_dir, showWarnings = FALSE, recursive = TRUE)
  
  # --- 1. raw counts -> filtered_feature_bc_matrix.h5 ---
  cat("  writing counts matrix...\n")
  counts <- tryCatch(
    GetAssayData(obj, assay = "Spatial", layer = "counts"),   # Seurat v5 style
    error = function(e) GetAssayData(obj, assay = "Spatial", slot = "counts")  # older style
  )
  
  h5_path <- file.path(out_dir, "filtered_feature_bc_matrix.h5")
  if (file.exists(h5_path)) file.remove(h5_path)  # write10xCounts refuses to overwrite silently
  
  DropletUtils::write10xCounts(
    path        = h5_path,
    x           = counts,
    barcodes    = colnames(counts),
    gene.id     = rownames(counts),
    gene.symbol = rownames(counts),
    type        = "HDF5",
    version     = "3",
    overwrite   = TRUE
  )
  
  # --- 2. spatial coordinates -> tissue_positions.csv ---
  cat("  writing tissue_positions.csv...\n")
  if (length(Images(obj)) == 0) {
    stop("Sample '", sample_name, "' has no image/coordinate slot - cannot export spatial ",
         "position data. (This should have been caught earlier as a warning during loading.)")
  }
  image_name <- Images(obj)[1]
  img_obj <- obj@images[[image_name]]
  
  coords <- img_obj@coordinates
  # Seurat's VisiumV1 coordinates columns: tissue, row, col, imagerow, imagecol
  # Space Ranger's tissue_positions.csv columns (with header, current format):
  #   barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
  tissue_positions <- data.frame(
    barcode              = rownames(coords),
    in_tissue            = coords$tissue,
    array_row            = coords$row,
    array_col            = coords$col,
    pxl_row_in_fullres   = coords$imagerow,
    pxl_col_in_fullres   = coords$imagecol,
    stringsAsFactors     = FALSE
  )
  write.csv(tissue_positions, file.path(spatial_dir, "tissue_positions.csv"),
            row.names = FALSE, quote = FALSE)
  
  # --- 3. scale factors -> scalefactors_json.json ---
  cat("  writing scalefactors_json.json...\n")
  sf <- img_obj@scale.factors
  # "scalefactors" is a plain S3 list in Seurat (not S4) - use $ accessors.
  # Guard with a fallback in case a given .rda was saved under a slightly
  # different naming/class convention for this slot.
  get_sf <- function(sf, name) {
    val <- tryCatch(sf[[name]], error = function(e) NULL)
    if (is.null(val)) {
      warning("scale.factors field '", name, "' not found for sample '", sample_name,
              "' - writing NA. Downstream plotting scale may be off; check the source object.")
      return(NA)
    }
    val
  }
  scalefactors <- list(
    tissue_hires_scalef       = get_sf(sf, "hires"),
    tissue_lowres_scalef      = get_sf(sf, "lowres"),
    fiducial_diameter_fullres = get_sf(sf, "fiducial"),
    spot_diameter_fullres     = get_sf(sf, "spot")
  )
  jsonlite::write_json(scalefactors, file.path(spatial_dir, "scalefactors_json.json"),
                       auto_unbox = TRUE, digits = NA)
  
  # --- 4. tissue image -> tissue_hires_image.png / tissue_lowres_image.png ---
  cat("  writing tissue image(s)...\n")
  img_arr <- img_obj@image
  if (is.null(img_arr)) {
    warning("Sample '", sample_name, "' has no image array in its image slot - ",
            "writing scalefactors/coordinates only. Squidpy's sq.read.visium() ",
            "will likely fail without a tissue image present; you may need to source ",
            "the original tissue image separately for this sample.")
  } else {
    # We only have ONE resolution baked into the Seurat object (Load10X_Spatial
    # keeps just one). We write it to both filenames rather than relying on
    # the Python script's hires-from-lowres fallback, since we can't reliably
    # tell which resolution this array actually is without the original files.
    png::writePNG(img_arr, file.path(spatial_dir, "tissue_hires_image.png"))
    png::writePNG(img_arr, file.path(spatial_dir, "tissue_lowres_image.png"))
  }
  
  cat("  done ->", out_dir, "\n")
  cat("=============== Export complete for", sample_name, "===============\n\n")
  
  invisible(out_dir)
}

#' Convenience wrapper: export every .rda file in a folder.
#'
#' @param rda_folder folder containing .rda files (e.g. "CRC_project/rda_raw")
#' @param output_root where "<sample_name>/" folders will be created
#'                     (point this at the same INPUT_FOLDER your Python script uses,
#'                      so exported samples sit alongside real Space Ranger folders)
export_all_rda <- function(rda_folder, output_root) {
  rda_files <- list.files(rda_folder, pattern = "\\.rda$", full.names = TRUE)
  if (length(rda_files) == 0) {
    stop("No .rda files found in: ", rda_folder)
  }
  cat("Found", length(rda_files), ".rda file(s) to export.\n\n")
  results <- lapply(rda_files, export_rda_to_spaceranger, output_root = output_root)
  cat("\nAll exports complete. ", length(results), " sample folder(s) written to: ", output_root, "\n", sep = "")
  invisible(results)
}

# ============================================================
# Example usage:
#
#   source("export_rda_to_spaceranger.R")
#
#   # single file:
#   export_rda_to_spaceranger("CRC_project/rda_raw/ST-colon1.rda",
#                              output_root = "CRC_project/GSE326101")
#
#   # everything in a folder:
#   export_all_rda("CRC_project/rda_raw", output_root = "CRC_project/GSE326101")
#
# After this, your existing cell2location.py script needs NO changes -
# it will find these exported samples as ordinary Space Ranger folders
# alongside any real ones already in INPUT_FOLDER.
# ============================================================