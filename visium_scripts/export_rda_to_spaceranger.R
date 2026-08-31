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

source("visium_scripts/rda_conversion.R")   # brings in load_spatial_object()
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

#' Get full-resolution pixel coordinates for an image slot, working across
#' Seurat's different spatial image classes (VisiumV1, VisiumV2, etc.)
#' without needing to know each class's internal slot structure.
#'
#' Why this matters: VisiumV1 (older Seurat) stores coordinates in an
#' @coordinates slot with columns imagerow/imagecol. VisiumV2 (Seurat >=5.1,
#' used for CytAssist/newer Visium data) has NO @coordinates slot at all -
#' you have to go through the GetTissueCoordinates() generic instead, which
#' returns x/y/cell columns. This helper normalizes both cases to a
#' consistent barcode/imagerow/imagecol data.frame so the rest of the
#' pipeline doesn't care which image class it's given.
#'
#' IMPORTANT: scale = NULL is required for BOTH classes to get true full-
#' resolution pixel coordinates (matching the raw TIFF and Space Ranger's
#' pxl_row_in_fullres/pxl_col_in_fullres) - VisiumV1's default is actually
#' scale="lowres", which would silently return coordinates scaled DOWN to
#' match the lowres image instead.
get_full_res_coords <- function(obj, image_name = NULL) {
  if (is.null(image_name)) image_name <- Images(obj)[1]
  
  raw <- GetTissueCoordinates(obj, image = image_name, scale = NULL)
  cn <- colnames(raw)
  
  if (all(c("imagerow", "imagecol") %in% cn)) {
    row_col <- "imagerow"; col_col <- "imagecol"
  } else if (all(c("x", "y") %in% cn)) {
    # VisiumV2 convention: x = imagecol (horizontal), y = imagerow (vertical)
    row_col <- "y"; col_col <- "x"
  } else {
    stop("Unrecognized coordinate column names from GetTissueCoordinates(): ",
         paste(cn, collapse = ", "), " - this image class may need explicit handling.")
  }
  
  barcodes <- if ("cell" %in% cn) raw[["cell"]] else rownames(raw)
  
  out <- data.frame(
    barcode   = barcodes,
    imagerow  = raw[[row_col]],
    imagecol  = raw[[col_col]],
    stringsAsFactors = FALSE
  )
  rownames(out) <- out$barcode
  out
}

#' Write an already-loaded, in-memory Seurat spatial object out to a
#' Space-Ranger-style folder. This is the shared "writer" used both by
#' export_rda_to_spaceranger() (loads from .rda first) and by
#' split_dual_tissue.R (subsets an object in-memory first) - factored out
#' so both paths produce byte-identical folder structures.
#'
#' @param obj a Seurat object with a Spatial assay + at least one image
#' @param sample_name name used for the output subfolder
#' @param output_root directory under which "<sample_name>/" will be created
write_seurat_to_spaceranger <- function(obj, sample_name, output_root) {
  
  if (!requireNamespace("DropletUtils", quietly = TRUE)) {
    stop("Package 'DropletUtils' is required to write the 10x-format .h5 file.\n",
         "Install it with:\n",
         '  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")\n',
         '  BiocManager::install("DropletUtils")')
  }
  
  cat("=============== Exporting", sample_name, "===============\n")
  
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
         "position data.")
  }
  image_name <- Images(obj)[1]
  img_obj <- obj@images[[image_name]]
  
  # Uses GetTissueCoordinates() generic (works across VisiumV1/VisiumV2/etc.)
  # rather than raw @coordinates slot access, which only exists on VisiumV1.
  coords <- get_full_res_coords(obj, image_name)
  # IMPORTANT for split samples: after subset(), coords may still contain
  # rows for barcodes that got dropped from the assay - restrict to exactly
  # the cells actually present in this object.
  coords <- coords[coords$barcode %in% colnames(obj), , drop = FALSE]
  
  # array_row/array_col/in_tissue live in VisiumV1's raw @coordinates slot,
  # but VisiumV2 (newer Seurat, CytAssist/HD-style data) doesn't expose them
  # via any generic - try the raw slot as a bonus, fall back gracefully if
  # it's not there rather than failing the whole export.
  extra <- tryCatch({
    raw_coords <- img_obj@coordinates
    raw_coords[coords$barcode, c("tissue", "row", "col")]
  }, error = function(e) NULL)
  
  if (!is.null(extra)) {
    in_tissue_vals <- extra$tissue
    array_row_vals <- extra$row
    array_col_vals <- extra$col
  } else {
    warning("Sample '", sample_name, "': could not retrieve array_row/array_col/in_tissue ",
            "(image class '", class(img_obj)[1], "' doesn't expose these the way VisiumV1 does). ",
            "Writing in_tissue=1 for all spots (objects normally only retain in-tissue spots ",
            "anyway) and array_row/array_col as NA. This does NOT affect pixel-coordinate-based ",
            "plotting or cell2location deconvolution - only tools that specifically need the ",
            "hex-grid array position would be affected.")
    in_tissue_vals <- 1
    array_row_vals <- NA
    array_col_vals <- NA
  }
  
  # Seurat's GetTissueCoordinates output columns: barcode, imagerow, imagecol (normalized
  # by get_full_res_coords() regardless of underlying image class)
  # Space Ranger's tissue_positions.csv columns (with header, current format):
  #   barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
  tissue_positions <- data.frame(
    barcode              = coords$barcode,
    in_tissue            = in_tissue_vals,
    array_row            = array_row_vals,
    array_col            = array_col_vals,
    pxl_row_in_fullres   = coords$imagerow,
    pxl_col_in_fullres   = coords$imagecol,
    stringsAsFactors     = FALSE
  )
  write.csv(tissue_positions, file.path(spatial_dir, "tissue_positions.csv"),
            row.names = FALSE, quote = FALSE)
  
  # --- 3. scale factors -> scalefactors_json.json ---
  cat("  writing scalefactors_json.json...\n")
  # ScaleFactors() is Seurat's generic accessor, working across image classes -
  # try that first, then fall back to raw @scale.factors (VisiumV1-only) for
  # older objects, so neither class silently breaks this step.
  sf <- tryCatch(ScaleFactors(img_obj), error = function(e) {
    tryCatch(img_obj@scale.factors, error = function(e2) NULL)
  })
  if (is.null(sf)) {
    warning("Sample '", sample_name, "': could not retrieve scale factors from image class '",
            class(img_obj)[1], "' via either ScaleFactors() or @scale.factors. ",
            "Writing all fields as NA - downstream plotting scale will likely be wrong ",
            "until this is fixed manually.")
    sf <- list()
  }
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
  # GetImage() is Seurat's generic for the raster itself, works across
  # VisiumV1/VisiumV2 - falls back to raw @image slot (VisiumV1-only) if
  # the generic isn't available for some reason.
  img_arr <- tryCatch(GetImage(img_obj, mode = "raster"), error = function(e) {
    tryCatch(img_obj@image, error = function(e2) NULL)
  })
  if (is.null(img_arr)) {
    warning("Sample '", sample_name, "' has no retrievable image array (tried GetImage() and ",
            "@image on class '", class(img_obj)[1], "') - writing scalefactors/coordinates only. ",
            "You'll need to source the tissue image separately for this sample.")
  } else {
    # GetImage(mode="raster") returns an R "raster" object (character matrix of
    # hex colors), not a numeric array - convert to a numeric array png::writePNG can use.
    if (inherits(img_arr, "raster")) {
      img_arr <- grDevices::as.raster(img_arr)
      img_dims <- dim(img_arr)
      rgb_vals <- grDevices::col2rgb(img_arr) / 255
      img_arr <- array(0, dim = c(img_dims[1], img_dims[2], 3))
      img_arr[,,1] <- matrix(rgb_vals[1, ], nrow = img_dims[1], byrow = TRUE)
      img_arr[,,2] <- matrix(rgb_vals[2, ], nrow = img_dims[1], byrow = TRUE)
      img_arr[,,3] <- matrix(rgb_vals[3, ], nrow = img_dims[1], byrow = TRUE)
    }
    # Note for split samples: this writes the FULL original image (both
    # tissue pieces still visible), just with tissue_positions.csv now only
    # listing spots for THIS piece. That's fine - Space Ranger/squidpy don't
    # require every pixel to be "used"; downstream tools only plot the listed
    # spots. If you want a visually cropped image instead, that's a separate,
    # optional step - not needed for cell2location/analysis correctness.
    png::writePNG(img_arr, file.path(spatial_dir, "tissue_hires_image.png"))
    png::writePNG(img_arr, file.path(spatial_dir, "tissue_lowres_image.png"))
  }
  
  cat("  done ->", out_dir, "\n")
  cat("=============== Export complete for", sample_name, "===============\n\n")
  
  invisible(out_dir)
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
  
  # --- load + patch (reuses everything we already fixed for the R QC pipeline) ---
  obj <- load_spatial_object(rda_path, sample_name)
  
  write_seurat_to_spaceranger(obj, sample_name = sample_name, output_root = output_root)
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